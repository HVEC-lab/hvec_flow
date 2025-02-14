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


@pyt.mark.parametrize(
    "h1, h2, Bfl, hsill, Cd, expected",
    [  (0.5, -0.05, 7, -3.05, 1  , 68.984)
     , (0.5, -2.7 , 7, -3.05, 0.9, 71.842)]
)
def test_capacity(h1, h2, Bfl, hsill, Cd, expected):
    res = sl.capacity(h1, h2, Bfl, hsill, Cd)
    assert np.isclose(res, expected, atol = tol, rtol = tol)