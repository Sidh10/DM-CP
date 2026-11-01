import importlib
from pathlib import Path

import pytest

PKG = Path(__file__).parents[1] / "genomath"
MODULES = sorted(p.stem for p in PKG.glob("*.py") if p.stem != "__init__")


@pytest.mark.parametrize("name", MODULES)
def test_module_imports(name):
    importlib.import_module(f"genomath.{name}")


@pytest.mark.parametrize("name", MODULES)
def test_core_is_ui_and_oracle_free(name):
    # genomath must stay UI-free, and library oracles belong in tests only
    src = (PKG / f"{name}.py").read_text(encoding="utf-8")
    for banned in ("streamlit", "networkx", "sympy", "Bio", "math.gcd"):
        assert f"import {banned}" not in src and f"from {banned}" not in src
