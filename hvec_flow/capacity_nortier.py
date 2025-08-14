"""
Simple capacity model for structures. Basis is Nortier.

HVEC-lab, 2025.
"""

import numpy as np
from constants import g


def _hcr(h1, hsill):
    """
    Internal function for critical water level
    """
    num = 2 * h1 + hsill
    res = num / 3
    return res


def Q(hin, hout, hroof, hsill, W, mu=1):
    """
    Capacity of rectangular structure based on Nortier, 1990.
    Covering all combinations of water levels and roof top levels.

    The function takes into account the direction of the head difference,
    rendering negative discharges for flow in the other direction.

    Args:
    hin: inside water level
    hout: outside water level
    hroof: level of roof
    hsill: level of sill
    W: flow width
    mu: discharge coefficient
    """
    # Detect flow direction
    dir = np.sign(hin - hout)

    # Internally, h1 is the highest water level and h2 the lowest
    # Direction is taken into account at the end.
    h1 = np.maximum(hin, hout)
    h2 = np.minimum(hin, hout)

    # Early return if water levels below sill level
    if h1 <= hsill:
        res = 0
        return res

    hcr = _hcr(h1, hsill)

    if h2 < hcr:
        deltaH = h1 - np.minimum(hcr, hroof) # Velocity head
        d = np.minimum(hcr, hroof) - hsill
    else:
        deltaH = h1 - h2
        d = np.minimum(h2, hroof) - hsill
    
    # Average flow velocity
    U = mu * np.sqrt(2 * g * deltaH)

    # Flow area
    A = d * W

    # Discharge
    res = dir * U * A

    return res
