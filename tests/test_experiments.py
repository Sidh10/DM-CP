import pytest
from pathlib import Path

from experiments.exp1_mutation_sets import run as run_exp1


def test_exp1_all_expected_values_match():
    checks, _ = run_exp1()
    failed = [c.item for c in checks if c.match is False]
    assert not failed, failed


def test_exp2_all_expected_values_match():
    from experiments.exp2_counting_enumeration import run as run_exp2
    checks, _ = run_exp2()
    failed = [c.item for c in checks if c.match is False]
    assert not failed, failed


def test_exp3_all_expected_values_match():
    from experiments.exp3_algebraic_verification import run as run_exp3
    failed = [c.item for c in run_exp3() if c.match is False]
    assert not failed, failed


def test_exp4_all_expected_values_match():
    from experiments.exp4_modular_indexing import run as run_exp4
    checks, _ = run_exp4()
    failed = [c.item for c in checks if c.match is False]
    assert not failed, failed


@pytest.mark.parametrize("name", ["exp5_network_structure", "exp6_paths_and_trees", "exp7_tree_analysis",
                                 "exp8_real_data"])
def test_exp5_exp6_all_expected_values_match(name):
    import importlib
    checks, _ = importlib.import_module(f"experiments.{name}").run()
    failed = [c.item for c in checks if c.match is False]
    assert not failed, failed


def test_documentation_figures_match_the_code():
    """Every number quoted in the Stage 9 deliverables is recomputed; documentation must not drift from the code."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "check_report_figures", Path(__file__).resolve().parents[1] / "scripts" / "check_report_figures.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    try:
        module.main()
    except SystemExit as exc:
        raise AssertionError(str(exc)) from None
    assert not module.FAILED
    assert len(module.PASSED) >= 75
