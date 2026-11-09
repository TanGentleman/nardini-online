"""The offline test environments must match the Modal images they stand in for."""

import re
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent
COMMON = (ROOT / "src/app/common.py").read_text()
WEB_IMAGE, _, WORKER_IMAGE = COMMON.partition("nardini_image =")
PIN = re.compile(r'"([A-Za-z0-9_.-]+)(?:\[[a-z,]+\])?==([^"]+)"')


def image_pins(block: str) -> dict[str, str]:
    return dict(PIN.findall(block))


def dev_pins() -> dict[str, str]:
    groups = tomllib.loads((ROOT / "pyproject.toml").read_text())["dependency-groups"]
    return dict(p.split("==") for p in groups["dev"])


def worker_pins() -> dict[str, str]:
    lines = (ROOT / "tests/requirements-worker.in").read_text().splitlines()
    return dict(line.split("==") for line in lines if "==" in line)


@pytest.mark.parametrize(
    ("image", "test_pins"),
    [(WEB_IMAGE, dev_pins), (WORKER_IMAGE, worker_pins)],
    ids=["web", "worker"],
)
def test_test_env_pins_match_image(image, test_pins):
    pins = image_pins(image)
    assert pins
    assert {name: test_pins().get(name) for name in pins} == pins


@pytest.mark.parametrize(
    ("image", "python"), [(WEB_IMAGE, "3.12"), (WORKER_IMAGE, "3.9")]
)
def test_image_python_version(image, python):
    assert f'python_version="{python}"' in image
