"""
Package hvec_flow, subpackage hydrostatics.

HVEC-lab, May 2022
"""

from hvec_flow.constants import g


def p(h, z, rhow = 1025):
    """
    Hydrostatic pressure

    Parameters
    -------
    h: array_like. Water level with respect to reference
    z: array_like. Level with respect to reference
    rhow: array_like. Water density. Default assumption is salt water

    Returns
    -------
    p : float. Hydrostatic pressure

    See Also
    --------
    NA

    Examples
    --------
    NA

    Raises
    --------
    No warning messages are specified. Informed user is assumed

    Issues
    --------
    - Definition of h different from literature. Specified with respect to 
        reference level
    - Provide English reference
 
    References
    --------
    Battjes (1989) - Fluid mechanics, lecture notes Delft University (in Dutch)

    """
    return rhow * g * (h-z)


def Fstatic(h, ztop, zbot, B = 1, rhow = 1025):
    """
    Hydrostatic horizontal force

    Parameters
    -------
    h: array_like. Water level with respect to reference
    z: array_like. Level with respect to reference
    rhow: array_like. Water density. Default assumption is salt water

    Returns
    -------
    F_static : float. Hydrostatic force

    See Also
    --------
    NA

    Examples
    --------
    NA

    Raises
    --------
    No warning messages are specified. Informed user is assumed

    Issues
    --------
    - Definition of h different from literature. Specified with respect to 
        reference level
    - Provide English reference
 
    References
    --------
    Battjes (1989) - Fluid mechanics, lecture notes Delft University (in Dutch)
    """
    ptop = p(h = h, z = ztop, rhow = rhow)
    pbot = p(h = h, z = zbot, rhow = rhow)
    d = ztop - zbot
    Fstat = ((ptop + pbot)/2) * d * B
    return Fstat
