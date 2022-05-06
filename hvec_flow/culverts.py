"""
Package hvec_flow, subpackage culverts.

Function names are kept short. Use namespace to distinguish between
culverts (this sub-package) and other structures.

The "vectorize" decorator from numpy has been used to allow for vector
input. Ensure that only one input is a vector to prevent ill-defined output.

Prepared by HVEC lab, 2022
"""

import numpy as np
from constants import g


@np.vectorize
def Q(muA, delta_H, g = g):
    """
    Discharge in a simple culvert using Bernoulli

    Parameters
    -------
    muA: array_like. Effective flow area (area including losses)
    delta_H : array_like. Difference in energy head over the culvert
    g: acceleration of gravity; default taken from sub-package constants

    Returns
    -------
    Q : float
        Discharge

    Issues
    --------
    - Provide English reference
 
    References
    --------
    Battjes (1989) - Fluid mechanics, lecture notes Delft University (in Dutch)
    """

    assert muA >= 0, 'Error in Q. Effective area is negative'
    assert g > 0, 'Error in Q. Acceleration of gravity is not positive'

    dir = np.sign(delta_H)                 # Direction of flow
    U = np.sqrt(2 * g * abs(delta_H))   # Profile averaged velocity (value)
    Q = dir * U * muA                   # Discharge with direction
    return Q


@np.vectorize
def delta_H(muA, Q, g = g):
    """
    Head difference over a simple culvert using Bernoulli

    Parameters
    -------
    muA: array_like. Effective flow area (area including losses)
    Q : array_like. Discharge
    g: acceleration of gravity; default taken from sub-package constants

    Returns
    -------
    delta_H : float
        Energy head difference

    Issues
    --------
    - Provide English reference
 
    References
    --------
    Battjes (1989) - Fluid mechanics, lecture notes Delft University (in Dutch)
    """
    assert muA >= 0, 'Error in delta_H. Effective area is negative'
    assert muA != 0, 'Error in delta_H. Head can not be calculated for zero area'
    assert g > 0, 'Error in delta_H. Acceleration of gravity is not positive'

    U = Q / muA
    cnt = U * np.abs(U)  # Flow velocity including direction
    denom = 2 * g
    delta_H = cnt / denom
    return delta_H


@np.vectorize
def muA_from_Q(Q, delta_H, g = g):
    """
    Effective culvert area from prescribed flow and energy head difference.
    
    Parameters
    -------
    Q: array_like. Discharge
    delta_H : array_like. Difference in energy head over the culvert

    Returns
    -------
    muA: float; effective flow area

    Issues
    --------
    - Provide English reference
 
    References
    --------
    Battjes (1989) - Fluid mechanics, lecture notes Delft University (in Dutch)
    """
    assert g > 0, 'Error in muA_from_Q: Acceleration of gravity is not positive'
    assert Q * delta_H >= 0, 'Error in muA_from_Q. Q and H should have the same sign'
    assert delta_H != 0, 'Error in muA_from_Q. Area can not be determined for zero head'
    
    Q_int = np.abs(Q)
    delta_H_int = np.abs(delta_H)
    U = np.sqrt(2 * g * delta_H_int)
    muA = Q_int / U    
    return muA


def muA_series(muA):
    """
    Capacity expressed in effective area; series system of culverts.

    Calculating the effective area of a series of culverts is a rapid
    and convenient way of quantifying system capacity.

    Parameters
    -------
    muA: array_like. Effective flow area (area including losses)

    Returns
    -------
    muA_series: float. Capacity of the full system expressed in 
        effective area

    Issues
    --------
    - Provide English reference
 
    References
    --------
    Based on:
    Battjes (1989) - Fluid mechanics, lecture notes Delft University (in Dutch)
    """
    muA = np.array(muA)

    assert np.min(muA) >= 0, 'Error in muA_series. All areas should be non-negative'

    if np.min(muA) > 0:
        var1 = 1/(muA ** 2)
        sum = var1.sum()
        res = 1 / np.sqrt(sum)
    else:
        res = 0
    return res


def muA_par(muA):
    """
    Effective area of parallel culverts

    Parameters
    -------
    muA: array_like. Effective flow area (area including losses)

    Returns
    -------
    muA_par: float. Capacity of the full system expressed in 
        effective area

    Issues
    --------
    - Provide English reference
 
    References
    --------
    Based on:
    Battjes (1989) - Fluid mechanics, lecture notes Delft University (in Dutch)
    """
    muA = np.array(muA)

    assert np.min(muA) >= 0, 'Error in muA_par. All areas should be non-negative'

    res = muA.sum()
    return res


def muA_system(muA):
    """
    Effective area of an arbitrary system of culverts.

    Culverts are organised in a "list of lists". The elements
    in the main list are assumed to describe parallel systems of 
    culverts. From this, a list of single-valued effective areas
    is formed by calling muA_par.

    Finally, muA_series is invoked to calculate the total effective
    area.

    Parameters
    -------
    muA: list of lists, with every sublist describing a parallel system

    Returns
    -------
    muA_system: float. Capacity of the full system expressed in 
        effective area

    Issues
    --------
    - Provide English reference
 
    References
    --------
    Based on:
    Battjes (1989) - Fluid mechanics, lecture notes Delft University (in Dutch)
    """

    # Transform all parallel systems to their singular equivalent...
    muA_eq = []
    for row in muA:
        muA_eq.append(muA_par(row))

    # ... then run muA_series 
    res = muA_series(muA_eq)
    return res    
