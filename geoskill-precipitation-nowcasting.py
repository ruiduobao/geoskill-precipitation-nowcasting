#!/usr/bin/env python3
"""precipitation-nowcasting — 降水临近预报

基于**光流法**的拉格朗日持久性（Lagrangian persistence）降水临近预报：

- **位移估计**（光流）：用归一化互相关在搜索窗内匹配相邻两帧降水场，
  估计场体的平移矢量 (vy, vx)；对多个相邻帧对取平均以提高稳健性
  （简化交叉相关光流，等价于全场 Lucas-Kanade 的平移假设）。
- **拉格朗日外推**：假设位移场在短时间内保持不变，将最新一帧沿估计速度
  平移，外推未来 0–6 小时的降水场（scipy.ndimage.shift 双线性重采样）。

数据源：本地多期 GeoTIFF（每波段 = 一个时间步雷达/卫星降水场），或
``--synthetic`` 生成以已知速度平移的模拟降水场用于离线验证。

隐私声明 / Privacy：
- 默认完全离线运行，``--synthetic`` 模式不读取任何外部数据。
- 所有处理在本地完成，不上传任何用户数据。

Usage:
    python precipitation-nowcasting.py --bbox 116 39 117 40 --lead-time 60 --n-frames 4
    python precipitation-nowcasting.py --input radar_stack.tif --lead-time 90

License: MIT
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.ndimage import shift as _nd_shift

VERSION = "1.0.0"
SKILL_NAME = "precipitation-nowcasting"

# ---- 复用共享核心库（本地 vendored，随脚本目录一起分发）----
try:
    from _geoskill_core.errors import (
        GeoSkillError, UsageError, ValidationError, ProcessError, to_exit_code,
    )
    from _geoskill_core.manifest import OutputManifest, OutputFile
except ImportError:  # pragma: no cover - fallback minimal definitions
    class GeoSkillError(Exception):
        def __init__(self, message: str, code: int = 7, kind: str = "EGeo", **kw):
            super().__init__(message)
            self.message, self.code, self.kind = message, code, kind

    class UsageError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=2, kind="EUsage", **k)

    class ValidationError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=6, kind="EValidate", **k)

    class ProcessError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=7, kind="EProcess", **k)

    def to_exit_code(exc):
        return getattr(exc, "code", 7)

    OutputManifest = None
    OutputFile = None


def _utc_now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat().replace("+00:00", "Z")


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------
def validate_bbox(bbox: List[float], source: str = "bbox") -> None:
    """Validate geographic bbox: W<=E, S<=N, lon/lat in range, min area.

    Cross-dateline (W>E) is a ValidationError with a hint to split.
    """
    if bbox is None or len(bbox) != 4:
        raise ValidationError(f"{source}: expected 4 floats [W S E N], got {bbox!r}")
    try:
        w, s, e, n = (float(v) for v in bbox)
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"{source}: non-numeric bbox values: {bbox!r}") from exc
    for v, name in ((w, "W"), (s, "S"), (e, "E"), (n, "N")):
        if not (v == v):
            raise ValidationError(f"{source}: bbox contains NaN at {name}")
    if not (-180.0 <= w <= 180.0 and -180.0 <= e <= 180.0):
        raise ValidationError(f"{source}: lon out of [-180,180]: W={w} E={e}")
    if not (-90.0 <= s <= 90.0 and -90.0 <= n <= 90.0):
        raise ValidationError(f"{source}: lat out of [-90,90]: S={s} N={n}")
    if w > e:
        raise ValidationError(
            f"{source}: W ({w}) > E ({e}); cross-dateline bboxes are not supported. "
            "Split into two bboxes on each side of the 180\u00b0 meridian and run separately."
        )
    if s > n:
        raise ValidationError(f"{source}: S ({s}) > N ({n})")
    if (e - w) < 1e-9 or (n - s) < 1e-9:
        raise ValidationError(
            f"{source}: bbox too small (dlon={e - w}, dlat={n - s}); need > 1e-9 degrees"
        )


# ---------------------------------------------------------------------------
# 核心算法
# ---------------------------------------------------------------------------
def _shift_array(arr: np.ndarray, dy: float, dx: float) -> np.ndarray:
    """将二维数组内容平移 (dy, dx) 个像元（整数），空出部分补零（无环绕）。"""
    out = np.zeros_like(arr)
    h, w = arr.shape
    dy_i, dx_i = int(round(dy)), int(round(dx))
    sy0, sy1 = max(0, -dy_i), min(h, h - dy_i)
    ty0, ty1 = max(0, dy_i), min(h, h + dy_i)
    sx0, sx1 = max(0, -dx_i), min(w, w - dx_i)
    tx0, tx1 = max(0, dx_i), min(w, w + dx_i)
    if ty1 > ty0 and tx1 > tx0:
        out[ty0:ty1, tx0:tx1] = arr[sy0:sy1, sx0:sx1]
    return out


def _normalized_corr(a: np.ndarray, b: np.ndarray) -> float:
    """零均值归一化互相关（Pearson r）。

    NaN-safe: only the pixels that are finite in *both* a and b contribute
    to the correlation. This keeps the optical-flow search well-defined
    when the input has NoData pixels (the correlation simply uses the
    common-mask sub-array).
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    mask = np.isfinite(a) & np.isfinite(b)
    if not mask.any():
        return 0.0
    am = a[mask] - a[mask].mean()
    bm = b[mask] - b[mask].mean()
    den = np.sqrt(np.sum(am * am) * np.sum(bm * bm))
    if den == 0:
        return 0.0
    return float(np.sum(am * bm) / den)


def estimate_displacement(
    frame_a: np.ndarray, frame_b: np.ndarray, search: int = 12
) -> Tuple[float, float, float]:
    """用互相关搜索估计 frame_a → frame_b 的平移矢量。

    返回 (vy, vx, corr)：frame_b ≈ shift(frame_a, vy, vx)，corr 为峰值相关系数。
    在 [-search, search] 的整数位移格点上搜索最大归一化互相关。
    """
    if frame_a.shape != frame_b.shape:
        raise ValidationError(
            f"frame shape mismatch: {frame_a.shape} vs {frame_b.shape}",
            shape_a=tuple(frame_a.shape), shape_b=tuple(frame_b.shape),
        )
    best_corr = -2.0
    best = (0.0, 0.0)
    for cy in range(-search, search + 1):
        for cx in range(-search, search + 1):
            sa = _shift_array(frame_a, cy, cx)
            c = _normalized_corr(sa, frame_b)
            if c > best_corr:
                best_corr = c
                best = (float(cy), float(cx))
    return best[0], best[1], best_corr


def estimate_motion(
    cube: np.ndarray, search: int = 12
) -> Dict[str, Any]:
    """对 (n_frames, H, W) 序列逐相邻帧对估计位移并取平均。

    返回 dict：vy, vx（平均速度，单位：像元/时间步）、pairs（逐对结果）、
    mean_corr（平均峰值相关）。
    """
    if cube.ndim != 3:
        raise ValidationError(
            f"cube must be 3-D (n_frames, H, W), got {cube.shape}",
            shape=tuple(cube.shape),
        )
    n = cube.shape[0]
    if n < 2:
        raise ValidationError(
            f"need at least 2 frames for motion estimation, got {n}", n=int(n),
        )
    pairs: List[Dict[str, Any]] = []
    vys: List[float] = []
    vxs: List[float] = []
    for k in range(n - 1):
        vy, vx, corr = estimate_displacement(cube[k], cube[k + 1], search=search)
        vys.append(vy); vxs.append(vx)
        pairs.append({"from_frame": k, "to_frame": k + 1,
                      "vy": vy, "vx": vx, "corr": corr})
    return {
        "vy": float(np.mean(vys)),
        "vx": float(np.mean(vxs)),
        "pairs": pairs,
        "mean_corr": float(np.mean([p["corr"] for p in pairs])),
    }


def extrapolate(
    frame: np.ndarray, vy: float, vx: float, n_steps: int
) -> np.ndarray:
    """拉格朗日持久性外推：将 frame 沿 (vy, vx) 平移 k 步，k=1..n_steps。

    返回 (n_steps, H, W) 预报序列（双线性重采样，边界外补零）。
    """
    if n_steps < 1:
        raise ValidationError(f"n_steps must be >= 1, got {n_steps}", n=int(n_steps))
    frame = np.asarray(frame, dtype=np.float32)
    out = np.zeros((n_steps, frame.shape[0], frame.shape[1]), dtype=np.float32)
    for k in range(1, n_steps + 1):
        out[k - 1] = _nd_shift(frame, (vy * k, vx * k), order=1,
                               mode="constant", cval=0.0).astype(np.float32)
    return out


# ---------------------------------------------------------------------------
# 合成数据：以已知速度平移的降水场
# ---------------------------------------------------------------------------
def generate_synthetic_cube(
    bbox: List[float],
    n_frames: int = 4,
    width: int = 64,
    height: int = 64,
    vy: float = 2.0,
    vx: float = 3.0,
    seed: int = 42,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """生成 (n_frames, H, W) 的平移降水场序列（高斯雨团 + 噪声）。

    雨团以 (vy, vx) 像元/帧 的速度平移；初始中心按总位移回退，
    保证整个序列中雨团基本位于画面内部。
    """
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:height, 0:width].astype(np.float32)
    cy0 = height / 2.0 - vy * (n_frames - 1) / 2.0
    cx0 = width / 2.0 - vx * (n_frames - 1) / 2.0
    sigma = 6.0
    amp = 30.0
    cube = np.zeros((n_frames, height, width), dtype=np.float32)
    for k in range(n_frames):
        cy = cy0 + vy * k
        cx = cx0 + vx * k
        g = amp * np.exp(-((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * sigma ** 2))
        noise = rng.normal(0, 0.5, (height, width)).astype(np.float32)
        cube[k] = np.clip(g + noise, 0.0, None)
    info = {
        "bbox": bbox, "width": width, "height": height,
        "n_frames": n_frames, "truth_vy": vy, "truth_vx": vx,
    }
    return cube, info


def make_truth_future(
    n_steps: int, cube_info: Dict[str, Any]
) -> np.ndarray:
    """供测试使用：按真实速度生成未来 n_steps 帧（与合成序列同一模型）。"""
    height = cube_info["height"]
    width = cube_info["width"]
    vy = cube_info["truth_vy"]
    vx = cube_info["truth_vx"]
    n_frames = cube_info["n_frames"]
    yy, xx = np.mgrid[0:height, 0:width].astype(np.float32)
    cy0 = height / 2.0 - vy * (n_frames - 1) / 2.0
    cx0 = width / 2.0 - vx * (n_frames - 1) / 2.0
    sigma = 6.0
    amp = 30.0
    out = np.zeros((n_steps, height, width), dtype=np.float32)
    for s in range(1, n_steps + 1):
        k = n_frames - 1 + s        # 最后一帧之后的第 s 步
        cy = cy0 + vy * k
        cx = cx0 + vx * k
        out[s - 1] = amp * np.exp(
            -((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * sigma ** 2))
    return out


# ---------------------------------------------------------------------------
# GeoTIFF I/O（直接用 rasterio，保证离线可用）
# ---------------------------------------------------------------------------
def write_geotiff(
    path: str,
    cube: np.ndarray,
    bbox: List[float],
    nodata: float = -9999.0,
) -> None:
    import rasterio
    from rasterio.transform import from_bounds

    if cube.ndim == 2:
        cube = cube[np.newaxis, ...]
    nb, h, w = cube.shape
    transform = from_bounds(bbox[0], bbox[1], bbox[2], bbox[3], w, h)
    profile = {
        "driver": "GTiff", "height": h, "width": w, "count": nb,
        "dtype": "float32", "crs": "EPSG:4326", "transform": transform,
        "nodata": nodata, "compress": "deflate",
    }
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with rasterio.open(path, "w", **profile) as dst:
        for b in range(nb):
            dst.write(cube[b].astype("float32"), b + 1)


def read_geotiff(path: str) -> Tuple[np.ndarray, List[float]]:
    import rasterio
    if not os.path.exists(path):
        raise UsageError(f"input raster not found: {path}", path=path)
    with rasterio.open(path) as src:
        cube = src.read().astype(np.float32)
        b = src.bounds
        bbox = [b.left, b.bottom, b.right, b.top]
    return cube, bbox


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------
def write_manifest(
    output_dir: str,
    args: argparse.Namespace,
    outputs: List[Dict[str, Any]],
    qa: Dict[str, Any],
    started_at: str,
    exit_code: int,
    bbox: List[float],
) -> Optional[str]:
    if OutputManifest is None:
        return None
    cmd = " ".join([SKILL_NAME] + sys.argv[1:])
    man = OutputManifest(
        skill=SKILL_NAME,
        skill_version=VERSION,
        command=cmd,
        started_at=started_at,
        finished_at=_utc_now(),
        exit_code=exit_code,
        inputs={
            "input": getattr(args, "input", None),
            "lead_time": getattr(args, "lead_time", None),
            "n_frames": getattr(args, "n_frames", None),
            "synthetic": bool(getattr(args, "synthetic", False)),
        },
        outputs=[OutputFile(**o) for o in outputs],
        qa=qa,
        software={"python": sys.version.split()[0], "skill": SKILL_NAME},
    )
    path = os.path.join(output_dir, "output-manifest.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(man.to_dict(), f, ensure_ascii=False, indent=2, default=str)
    return path


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def process(args: argparse.Namespace) -> int:
    started_at = _utc_now()
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    bbox = list(args.bbox) if args.bbox else None

    # CLI-layer validation: lead-time and dt-minutes must be strictly positive
    if not (args.lead_time == args.lead_time) or args.lead_time <= 0:
        raise ValidationError(
            f"--lead-time must be a positive number of minutes (got {args.lead_time!r}); "
            "a 0 lead time would produce a redundant 'now' forecast."
        )
    if not (args.dt_minutes == args.dt_minutes) or args.dt_minutes <= 0:
        raise ValidationError(
            f"--dt-minutes must be a positive number of minutes (got {args.dt_minutes!r})"
        )
    if not isinstance(args.search, int) or args.search < 1:
        raise ValidationError(
            f"--search must be a positive integer (got {args.search!r})"
        )

    synth_info: Optional[Dict[str, Any]] = None
    if args.input and not args.synthetic:
        cube, file_bbox = read_geotiff(args.input)
        bbox = bbox if bbox is not None else file_bbox
        # Replace NoData sentinel with NaN so the cross-correlation
        # displacement search doesn't see -9999 as a strong "feature".
        import rasterio as _rio
        with _rio.open(args.input) as _src:
            _nd = _src.nodata
        if _nd is not None:
            cube = np.where(cube == _nd, np.nan, cube).astype(np.float32)
        if not np.isfinite(cube).any():
            raise ValidationError(
                f"input raster '{args.input}' contains only NoData pixels; nothing to forecast"
            )
        source_note = args.input
    else:
        if bbox is None:
            raise UsageError("provide --bbox (synthetic mode) or --input <raster>")
        validate_bbox(bbox, source="--bbox")
        if not isinstance(args.n_frames, int) or args.n_frames < 2:
            raise ValidationError(
                f"--n-frames must be a positive integer >= 2 (got {args.n_frames!r}); "
                "at least 2 frames are needed for optical-flow displacement estimation"
            )
        cube, synth_info = generate_synthetic_cube(
            bbox, n_frames=args.n_frames,
        )
        source_note = "synthetic"

    if cube.size == 0:
        raise ValidationError("input raster is empty")
    if cube.ndim != 3:
        raise ValidationError(
            f"input must be a frame stack (n_frames, H, W), got {cube.shape}",
            shape=tuple(cube.shape),
        )
    if cube.shape[0] < 2:
        raise ValidationError(
            f"need at least 2 frames, got {cube.shape[0]}", n=int(cube.shape[0]),
        )

    # If --bbox is also given with --input, validate the user-supplied bbox
    if bbox is not None and args.bbox is not None:
        validate_bbox(bbox, source="--bbox")

    # 1) 光流位移估计
    motion = estimate_motion(cube, search=args.search)
    vy, vx = motion["vy"], motion["vx"]

    # 2) 预报步数
    n_forecast = max(1, int(round(args.lead_time / args.dt_minutes)))

    # 3) 拉格朗日外推
    last = cube[-1]
    forecast = extrapolate(last, vy, vx, n_forecast)

    # Mask NoData pixels in the forecast (forecast.tif is the latest frame
    # translated by (vy, vx); pixels that were NoData in the latest frame
    # stay NoData here so the user can identify them downstream).
    finite_last = np.isfinite(last)
    forecast_out = np.where(
        finite_last[np.newaxis, :, :], forecast, -9999.0,
    ).astype(np.float32)

    # 写出产物
    fc_tif = os.path.join(output_dir, "forecast.tif")
    write_geotiff(fc_tif, forecast_out, bbox, nodata=-9999.0)

    disp_path = os.path.join(output_dir, "displacement.json")
    disp_payload = {
        "source": source_note,
        "n_frames_input": int(cube.shape[0]),
        "dt_minutes": args.dt_minutes,
        "lead_time_minutes": args.lead_time,
        "n_forecast_steps": n_forecast,
        "displacement_per_step_px": {"vy": vy, "vx": vx},
        "mean_peak_correlation": motion["mean_corr"],
        "frame_pairs": motion["pairs"],
        "forecast_lead_times_minutes": [
            args.dt_minutes * (k + 1) for k in range(n_forecast)
        ],
    }
    if synth_info is not None:
        disp_payload["synthetic_truth_velocity_px"] = {
            "vy": synth_info["truth_vy"], "vx": synth_info["truth_vx"],
        }
    with open(disp_path, "w", encoding="utf-8") as f:
        json.dump(disp_payload, f, ensure_ascii=False, indent=2)

    qa: Dict[str, Any] = {
        "source": source_note,
        "n_frames_input": int(cube.shape[0]),
        "n_forecast_steps": n_forecast,
        "velocity_vy_px": vy,
        "velocity_vx_px": vx,
        "mean_peak_correlation": motion["mean_corr"],
        "lead_time_minutes": args.lead_time,
        "n_valid_pixels_last_frame": int(finite_last.sum()),
        "n_total_pixels_last_frame": int(finite_last.size),
    }
    if synth_info is not None:
        qa["synthetic_truth_velocity"] = {
            "vy": synth_info["truth_vy"], "vx": synth_info["truth_vx"],
        }

    outputs = [
        {"path": fc_tif, "kind": "raster", "crs_epsg": 4326,
         "bbox_wgs84": bbox, "band_count": int(forecast.shape[0])},
        {"path": disp_path, "kind": "json"},
    ]

    man_path = write_manifest(output_dir, args, outputs, qa, started_at, 0, bbox)

    if not args.quiet:
        print(f"[{SKILL_NAME}] source: {source_note}")
        print(f"[{SKILL_NAME}] estimated velocity: vy={vy:.2f} vx={vx:.2f} px/step "
              f"(corr={motion['mean_corr']:.3f})")
        print(f"[{SKILL_NAME}] lead time: {args.lead_time} min → "
              f"{n_forecast} forecast steps (dt={args.dt_minutes} min)")
        print(f"[{SKILL_NAME}] output: {fc_tif}")
        if man_path:
            print(f"[{SKILL_NAME}] manifest: {man_path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog=SKILL_NAME,
        description="Optical-flow Lagrangian-persistence precipitation nowcasting.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--bbox", nargs=4, type=float, metavar=("W", "S", "E", "N"),
                   help="geographic extent [minLon minLat maxLon maxLat]")
    p.add_argument("--input", help="input multi-band frame-stack GeoTIFF (band=frame)")
    p.add_argument("--lead-time", type=float, default=60.0,
                   help="total forecast lead time in minutes (default: 60)")
    p.add_argument("--dt-minutes", type=float, default=15.0,
                   help="time step between frames / forecasts in minutes (default: 15)")
    p.add_argument("--n-frames", type=int, default=4,
                   help="number of input frames for synthetic mode (default: 4)")
    p.add_argument("--search", type=int, default=12,
                   help="max displacement search radius in pixels (default: 12)")
    p.add_argument("--synthetic", action="store_true",
                   help="generate a synthetic translating precipitation field")
    p.add_argument("--output-dir", default="./output", help="output directory")
    p.add_argument("--quiet", action="store_true", help="suppress console output")
    p.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return process(args)
    except GeoSkillError as exc:
        print(f"[{SKILL_NAME}] ERROR [{exc.kind}] {exc.message}", file=sys.stderr)
        return exc.code
    except KeyboardInterrupt:
        return 130
    except Exception as exc:  # noqa: BLE001
        print(f"[{SKILL_NAME}] ERROR {exc}", file=sys.stderr)
        return to_exit_code(exc)


if __name__ == "__main__":
    sys.exit(main())
