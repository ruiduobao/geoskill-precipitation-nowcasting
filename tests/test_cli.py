"""CLI argument parsing tests for precipitation-nowcasting."""
import subprocess
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from conftest import mod, SCRIPT


def run_cli(args, timeout=120):
    return subprocess.run(
        [sys.executable, SCRIPT] + args,
        capture_output=True, text=True, timeout=timeout,
    )


class TestCLIBasics:
    def test_help(self):
        r = run_cli(["--help"])
        assert r.returncode == 0
        assert "--bbox" in r.stdout
        assert "--lead-time" in r.stdout
        assert "--synthetic" in r.stdout

    def test_version(self):
        r = run_cli(["--version"])
        assert r.returncode == 0
        assert "1.0.0" in r.stdout

    def test_no_args_exit_2(self):
        r = run_cli([])
        assert r.returncode == 2

    def test_invalid_lead_time_type_exit_2(self):
        """非法数值参数（非浮点）→ argparse error → exit 2。"""
        r = run_cli(["--bbox", "116", "39", "117", "40", "--lead-time", "notanumber"])
        assert r.returncode == 2

    def test_invalid_n_frames_type_exit_2(self):
        r = run_cli(["--bbox", "116", "39", "117", "40", "--n-frames", "x"])
        assert r.returncode == 2

    def test_input_file_not_found_exit_2(self):
        r = run_cli(["--input", "nonexistent_stack.tif"])
        assert r.returncode == 2

    def test_synthetic_without_bbox_exit_2(self):
        r = run_cli(["--synthetic"])
        assert r.returncode == 2


class TestCLISynthetic:
    def test_synthetic_runs(self):
        out = "./_test_cli_nowcast"
        r = run_cli([
            "--bbox", "116", "39", "117", "40", "--synthetic",
            "--lead-time", "60", "--n-frames", "4", "--quiet",
            "--output-dir", out,
        ])
        assert r.returncode == 0, r.stderr
        assert os.path.exists(os.path.join(
            os.path.dirname(SCRIPT), "_test_cli_nowcast", "output-manifest.json"))
