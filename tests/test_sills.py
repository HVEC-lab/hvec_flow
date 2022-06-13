"""
Tests for the submodule sills.

HVEC, May 2022
"""
# Public packages
import numpy as np
import pytest as pyt
import math
# Company packages
import hvec_flow.sills as sl


tol = 5e-4  # Accepted absolute tolerance of numerical results


@pyt.mark.parametrize(
    "Hup, hsill, qcr",
    [(1, -2, 8.859),
    (1, 1, 0),
    (1, 2, 0)
    ])
def test_q_crit(Hup, hsill, qcr):
    """
    Compare result to a set of hand calculations
    """
    assert (
        math.isclose(sl.q_crit(Hup, hsill), qcr,
        abs_tol = tol, rel_tol = tol)
    )


@pyt.mark.parametrize(
    "q, dcr",
    [(1.5, 0.612)
    ])
def test_d_crit(q, dcr):
    """
    Compare result to a set of hand calculations
    """
    assert (
        math.isclose(sl.d_crit(q), dcr,
        abs_tol = tol, rel_tol = tol)
    )