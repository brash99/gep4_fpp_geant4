#!/usr/bin/env python3
"""Numerical checks for the GEp-III analyzing-power implementation."""

import importlib.util
import math
from pathlib import Path


MODULE_PATH = Path(__file__).parents[1] / "analysis" / "analyzing_power.py"
SPEC = importlib.util.spec_from_file_location("analyzing_power", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def close(actual, expected, tolerance=1.0e-10):
    if not math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance):
        raise AssertionError(f"{actual} != {expected}")


ay, pt = MODULE.analyzing_power(6.7, 5.0)
close(pt, 0.5839434764093098)
close(ay, 0.07987478631776045)

ay, pt = MODULE.analyzing_power(6.7, 8.0)
close(pt, 0.9324597764324385)
close(ay, 0.03538488162658322)

ay, pt = MODULE.analyzing_power(6.7, 10.0)
assert pt > MODULE.PT_MAX_GEV
close(ay, 0.0)

print("GEp-III analyzing-power numerical checks passed")
