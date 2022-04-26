"""
Tests for the submodule culverts.

HVEC, April 2022
"""
# Public packages
import numpy as np
import pytest as pyt

# Company packages
import hvec_flow.culverts as cv


atol = 5e-4   # Accepted absolute tolerance of numerical results


@pyt.mark.parametrize(
    "muA, delta_H, g, Q",
    [(13, 0.75, cv.g, 49.868),    # Case with positive head and positive area
    (13, -0.75, cv.g, -49.868),   # Case with negative head and positive area
    (0, 0.75, cv.g, 0),             # Case with positive head and zero area
    (25, 0, cv.g, 0),               # Case with positive area and no head
    (25, 0.1, 4, 22.36068)          # Case with fantasy value of g
    ])
def test_Q(muA, delta_H, g, Q):
    """
    Compare result to a set of hand calculations
    """
    assert np.abs(Q - cv.Q(muA, delta_H, g)) < atol


@pyt.mark.parametrize(
    "muA, delta_H, Q",
    [(13, 0.75, 49.86818),    # Case with positive head and positive area
    (13, -0.75, -49.86818),   # Case with negative head and positive area
    (0, 0.75, 0),             # Case with positive head and zero area
    (25, 0, 0)                # Case with positive area and no head
    ])
def test_Q_defaults(muA, delta_H, Q):
    """
    Compare result to a set of hand calculations
    """
    assert np.abs(Q - cv.Q(muA, delta_H)) < atol


@pyt.mark.parametrize(
    "muA, Q, g, delta_H,",
    [(13, 50, cv.g, 0.75397),    # Case with positive head and positive area
    (13, -50, cv.g, -0.75397),   # Case with negative head and positive area
    (25, 0, cv.g, 0),               # Case with positive area and no head
    (25, 30, 4, 0.18)          # Case with fantasy value of g
    ])
def test_delta_H(muA, Q, g, delta_H):
    """
    Compare result to a set of hand calculations
    """
    assert np.abs(delta_H - cv.delta_H(muA, Q, g)) < atol


@pyt.mark.parametrize(
    "muA, Q, delta_H",
    [(13, 50, 0.75397),    # Case with positive head and positive area
    (13, -50, -0.75397),   # Case with negative head and positive area
    (25, 0, 0)                # Case with positive area and no head
    ])
def test_delta_H_defaults(muA, Q, delta_H):
    """
    Compare result to a set of hand calculations
    """
    assert np.abs(delta_H - cv.delta_H(muA, Q)) < atol


@pyt.mark.parametrize(
    "Q, delta_H, g, muA",
    [(40, 0.35, cv.g, 15.264),    # Case with positive head and positive area
    (-40, -0.35, cv.g, 15.264),   # Case with negative head and positive area
    (0, 0.35, cv.g, 0),             # Case with positive head and zero area
    (40, 0.35, 3, 27.603)          # Case with fantasy value of g
    ])
def test_muA_from_Q(Q, delta_H, g, muA):
    """
    Compare result to a set of hand calculations
    """
    assert np.abs(cv.muA_from_Q(Q, delta_H, g) - muA) < atol


@pyt.mark.parametrize(
    "Q, delta_H, muA",
    [(40, 0.35, 15.264),    # Case with positive head and positive area
    (-40, -0.35, 15.264),   # Case with negative head and positive area
    (0, 0.35, 0)              # Case with positive head and zero area
    ])
def test_muA_from_Q_defaults(Q, delta_H, muA):
    """
    Compare result to a set of hand calculations
    """
    assert np.abs(cv.muA_from_Q(Q, delta_H) - muA) < atol


@pyt.mark.parametrize(
    "muA, Q",
    [([25, 35, 10], 10),
    ([1000, 10, 1000], 15),
    ([0, 500, 250], 0),
    (250, 25)
    ])
def test_muA_series(muA, Q):
    """
    Function tested by returning to the original physics:
    - When all values in input are positive, summed head difference over
        all individual elements should equal head difference calculated from
        combined effective area
    - In case of at least one zero area in the input, combined area should be 
        zero
    - Input of only one value: result should equal that value
    - No hand calculation provided
    """
    if np.min(muA) == 0:
        assert cv.muA_series(muA) == 0
    else:
        delta_H_ind = cv.delta_H(muA = muA, Q = Q) # head loss individual elements
        delta_H_direct = np.sum(delta_H_ind)  # head loss system, direct
        
        muA_system = cv.muA_series(muA)
        delta_H_from_system = cv.delta_H(muA = muA_system, Q = Q)
        assert np.abs(delta_H_direct - delta_H_from_system) < atol


@pyt.mark.parametrize(
    "muA, muA_expected",
    [([25, 35, 10], 70),
    ([1000, 10, 1000], 2010),
    ([0, 500, 250], 750),
    (250, 250)
    ])
def test_muA_par(muA, muA_expected):
    """
    For parallel systems the system effective area is
    simply the addition of the areas. Nevertheless, the method
    is included for convenience.

    No hand calculation provided.
    """
    assert np.abs(cv.muA_par(muA) - muA_expected) < atol


@pyt.mark.parametrize(
    "muA, muA_expected",
    [([25, 35, 10], cv.muA_series([25, 35, 10])), # Pure series
    ([[25, 35, 10]], cv.muA_par([25, 35, 10])),   # Pure parallel
    ([[10, 20], [50, 75, 60], 100], 28.394319197) # Arbitrary system
    ])
def test_muA_system(muA, muA_expected):
    print(cv.muA_system(muA))
    assert np.abs(cv.muA_system(muA) - muA_expected) < atol