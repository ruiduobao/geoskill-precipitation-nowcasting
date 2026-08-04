"""Pytest fixtures for precipitation-nowcasting tests."""
import importlib.util
import os
import sys
import pytest

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(SKILL_DIR, "geoskill-precipitation-nowcasting.py")


def load_module():
    spec = importlib.util.spec_from_file_location("precipitation_nowcasting", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules["precipitation_nowcasting"] = module
    spec.loader.exec_module(module)
    return module


mod = load_module()


@pytest.fixture
def skill_dir():
    return SKILL_DIR


@pytest.fixture
def script_path():
    return SCRIPT
