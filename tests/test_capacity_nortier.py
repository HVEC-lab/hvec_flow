"""
Tests for sub-package capacity_nortier

HVEC-lab, 2025.
"""

import sys
import pytest as pyt
import numpy as np

sys.path.append('./hvec_flow')
import capacity_nortier as cn


@pyt.mark.parametrize(
      "hin, hout, hroof, hsill, B, mu, expected", 
      [   (-3.5, -4, 0, -3, 3, 0.9, 0)
        , (0, -0.5, 0, -3, 3, 0.9, 21.142)
        , (-0.5, 0, 0, -3, 3, 0.9, -21.142)
        , (2, 0.8333, 0, -3, 3, 0.9, 38.753)
        , (0.8333, 2, 0, -3, 3, 0.9, -38.753)
        , (0, -1.5, 0, -3, 3, 0.9, 23.919)
        , (-1.5, 0, 0, -3, 3, 0.9, -23.919)
        , (2, -0.16667, 0, -3, 3, 0.9, 50.740)
        , (-0.16667, 2, 0, -3, 3, 0.9, -50.740)
        ])
def test_Qnort(hin, hout, hroof, hsill, B, mu, expected):
    res = cn.Q(hin, hout, hroof, hsill, B, mu)
    assert np.isclose(res, expected, atol = 1e-4, rtol = 1e-4)
