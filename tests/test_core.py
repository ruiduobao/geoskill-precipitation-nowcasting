"""Core algorithm tests for precipitation-nowcasting."""
import sys
import os

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from conftest import mod


def _pearson(a, b):
    am = a - a.mean(); bm = b - b.mean()
    den = np.sqrt(np.sum(am * am) * np.sum(bm * bm))
    return float(np.sum(am * bm) / den) if den > 0 else 0.0


def _blob(h=64, w=64, cy=30.0, cx=28.0, sigma=6.0, amp=30.0):
    yy, xx = np.mgrid[0:h, 0:w].astype(np.float32)
    return amp * np.exp(-((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * sigma ** 2))


class TestEstimateDisplacement:
    def test_recovers_injected_shift(self):
        """互相关应精确恢复注入的整数平移 (3, 2)。"""
        a = _blob().astype(np.float32)
        b = mod._shift_array(a, 3, 2)
        vy, vx, corr = mod.estimate_displacement(a, b, search=8)
        assert (vy, vx) == (3.0, 2.0)
        assert corr > 0.9

    def test_negative_shift(self):
        a = _blob(cy=32, cx=32).astype(np.float32)
        b = mod._shift_array(a, -4, -2)
        vy, vx, _ = mod.estimate_displacement(a, b, search=8)
        assert (vy, vx) == (-4.0, -2.0)

    def test_zero_shift_identity(self):
        a = _blob().astype(np.float32)
        vy, vx, corr = mod.estimate_displacement(a, a, search=6)
        assert (vy, vx) == (0.0, 0.0)
        assert corr == pytest.approx(1.0, abs=1e-6)

    def test_shape_mismatch_raises(self):
        with pytest.raises(mod.ValidationError):
            mod.estimate_displacement(np.zeros((10, 10)), np.zeros((8, 8)))


class TestEstimateMotion:
    def test_recovers_synthetic_velocity(self):
        """多帧序列平均速度应接近注入的 (vy=2, vx=3)。"""
        cube, info = mod.generate_synthetic_cube(
            [116, 39, 117, 40], n_frames=4, vy=2.0, vx=3.0)
        motion = mod.estimate_motion(cube, search=8)
        assert motion["vy"] == pytest.approx(2.0, abs=1.0)
        assert motion["vx"] == pytest.approx(3.0, abs=1.0)
        assert len(motion["pairs"]) == 3
        assert motion["mean_corr"] > 0.8

    def test_too_few_frames_raises(self):
        with pytest.raises(mod.ValidationError):
            mod.estimate_motion(np.zeros((1, 8, 8)))

    def test_wrong_ndim_raises(self):
        with pytest.raises(mod.ValidationError):
            mod.estimate_motion(np.zeros((8, 8)))


class TestExtrapolate:
    def test_output_shape(self):
        frame = _blob().astype(np.float32)
        fc = mod.extrapolate(frame, 2.0, 3.0, 4)
        assert fc.shape == (4, 64, 64)

    def test_invalid_n_steps_raises(self):
        with pytest.raises(mod.ValidationError):
            mod.extrapolate(_blob().astype(np.float32), 1.0, 1.0, 0)

    def test_shifts_field_in_velocity_direction(self):
        """外推 1 步后，场体质心应沿速度方向移动。"""
        frame = _blob(cy=28, cx=28).astype(np.float32)
        fc = mod.extrapolate(frame, 3.0, 3.0, 1)
        yy, xx = np.mgrid[0:64, 0:64].astype(np.float32)
        cy0 = (frame * yy).sum() / frame.sum()
        cx0 = (frame * xx).sum() / frame.sum()
        cy1 = (fc[0] * yy).sum() / fc[0].sum()
        cx1 = (fc[0] * xx).sum() / fc[0].sum()
        assert cy1 - cy0 == pytest.approx(3.0, abs=0.5)
        assert cx1 - cx0 == pytest.approx(3.0, abs=0.5)


class TestNowcastSkill:
    def test_forecast_correlates_with_truth(self):
        """外推场应与按真实速度生成的未来场高度相关（位移估计有效）。"""
        cube, info = mod.generate_synthetic_cube(
            [116, 39, 117, 40], n_frames=4, vy=2.0, vx=3.0)
        motion = mod.estimate_motion(cube, search=8)
        fc = mod.extrapolate(cube[-1], motion["vy"], motion["vx"], 2)
        truth = mod.make_truth_future(2, info)
        r1 = _pearson(fc[0], truth[0])
        assert r1 > 0.6, f"step-1 forecast-truth correlation too low: {r1}"

    def test_estimated_velocity_close_to_truth(self):
        cube, info = mod.generate_synthetic_cube(
            [116, 39, 117, 40], n_frames=5, vy=1.0, vx=2.0)
        motion = mod.estimate_motion(cube, search=8)
        assert abs(motion["vy"] - info["truth_vy"]) < 1.0
        assert abs(motion["vx"] - info["truth_vx"]) < 1.0


class TestSynthetic:
    def test_cube_shape(self):
        cube, info = mod.generate_synthetic_cube(
            [116, 39, 117, 40], n_frames=4)
        assert cube.shape == (4, 64, 64)
        assert cube.min() >= 0.0

    def test_field_moves_over_time(self):
        """随时间推移，场体质心应沿 (vy, vx) 移动。"""
        cube, _ = mod.generate_synthetic_cube(
            [116, 39, 117, 40], n_frames=3, vy=2.0, vx=3.0)
        yy, xx = np.mgrid[0:64, 0:64].astype(np.float32)
        cx = [(cube[k] * xx).sum() / cube[k].sum() for k in range(3)]
        # x 质心随帧递增，增量 ~3
        assert cx[2] - cx[0] == pytest.approx(6.0, abs=1.0)


class TestGeoTiffIO:
    def test_roundtrip(self, tmp_path):
        cube = np.random.uniform(0, 5, (3, 16, 16)).astype(np.float32)
        bbox = [116.0, 39.0, 117.0, 40.0]
        path = str(tmp_path / "t.tif")
        mod.write_geotiff(path, cube, bbox)
        back, rb = mod.read_geotiff(path)
        assert back.shape == cube.shape
        np.testing.assert_allclose(rb, bbox, atol=1e-6)

    def test_read_missing_raises(self):
        with pytest.raises(mod.UsageError):
            mod.read_geotiff("/nonexistent/x.tif")
