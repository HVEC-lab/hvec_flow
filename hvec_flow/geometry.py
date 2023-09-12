"""
Package hvec_flow, subpackage geometry.

Prepared by HVEC lab, 2022
"""

from numpy import vectorize
import numpy as np

#TODO: deal with arbitary shapes using package shapely or sectionproperties


@vectorize
def d(h, z):
    """
    Local flow depth
    """
    return h - z


@vectorize
def Bsurf(h, z, Bbot, mslope_left=0, mslope_right=0):
    """
    Width water surface in trapezoidal profile

    Parameters
    ----------
    h : array_like. Water level with respect to reference
    z : array_like. Bottom level with the same reference and units as h
    Bbot: array_like. Width of channel bottom in the samen units as h and z
    mslope_left: slope of left bank (1:m), default vertical
    mslope_right: slope of right bank (1:m), default equal to left bank

    Returns
    -------
    Bsurf : float
        The area of the cross section of the flow
    """
    return Bbot + d(h, z) * (mslope_left + mslope_right)


@vectorize
def Ac(h, z, Bbot, mslope_left=0, mslope_right=0):
    """
    Cross section area of an open channel

    Parameters
    ----------
    h : array_like. Water level with respect to reference
    z : array_like. Bottom level with the same reference and units as h
    B: array_like. Width of channel in the samen units as h and z
    mslope_left: left slope, defined as 1 vertical in mslope horizontal
    mslope_right: right slope

    Returns
    -------
    Ac : float
        The area of the cross section of the flow
    """
    return (1/2)*(Bbot + Bsurf(h, z, Bbot, mslope_left, mslope_right)) * d(h, z)


@vectorize
def P(h, z, Bbot, mslope_left=0, mslope_right=0):
    """
    Perimeter of a channel. Closed conduit covered by providing a roof
    This is the internal version. The vectorized version is to be used

    Parameters
    ----------
    h : array_like. Water level with respect to reference
    z : array_like. Bottom level with the same reference and units as h
    Bbot: array_like. Width of channel in the samen units as h and z
    mslope_left: slope left bank; default vertical
    mslope_right: slope right bank; default vertical

    Returns
    -------
    P : float
        The perimeter of the flow
    """
    res = (
        Bbot +
        np.sqrt(1 + mslope_left**2) * d(h, z) +
        np.sqrt(1 + mslope_right**2) * d(h, z)
    )
    return res


@vectorize
def R(h, z, Bbot, mslope_left = 0, mslope_right = 0):
    """
    Hydraulic radius of an open channel

    Parameters
    ----------
    h : array_like. Water level with respect to reference
    z : array_like. Bottom level with the same reference and units as h
    B: array_like. Width of channel in the samen units as h and z
    
    Returns
    -------
    R : float
        The hydraulic radius of the flow

    References
    --------
    Battjes & Labeur, Unsteady flow in open channels, 2017
    """
    return Ac(h, z, Bbot, mslope_left, mslope_right)/P(h, z, Bbot, mslope_left, mslope_right)
