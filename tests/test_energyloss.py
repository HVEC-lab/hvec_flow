from hvec_flow import energyloss as loss
import math


def test_cf():
    R = 10
    k = 0.3
    expected = 0.004467
    result = loss.cf(R, k).item()
    assert math.isclose(expected, result, rel_tol = 1e-4, abs_tol = 1e-4)


def test_Chezy():
    R = 10
    k = 0.3
    expected = 46.862
    result = loss.Ch(R, k)
    assert math.isclose(expected, result, rel_tol = 1e-4, abs_tol = 1e-4)
