"""
Package hvec_flow, subpackage channelflow

Prepared by HVEC lab, 2022
"""


from numpy import vectorize


@vectorize
def deq(q, cf, ib):
    """
    Equilibrium depth of flow on a slope

    Parameters
    -------
    q: array_like. Discharge per meter running length
    cf: array_like. Friction factor
    ib: array_like. Bottom slope. Negative slope gives positive discharge

    Returns
    -------
    deq: array_like. Equilibrium depth

    Issues
    --------
    - Provide English reference
 
    References
    --------
    Battjes (1989) - Fluid mechanics, lecture notes Delft University (in Dutch)
    """
    try:
        res = ((cf * q * np.abs(q))/(-ib * g))**(1/3)
    except:
        raise ValueError(
            'Calculation equilibrium depth failed. Combinations of '
            'negative flow with positive slope and vice versa are allowed')
    return res


@vectorize
def Qeq(ifr, Ac, R, cf):
    """
    Equilibrium discharge in a prismatic channel of arbitrary shape

    Parameters
    -------
    ifr: array_like. Slope of energy line (friction slope)
    Ac: array_like. Area of channel cross section
    R: array_like. Hydraulic radius
    cf: array_like. Friction factor (White-Colebrook)

    Returns
    -------
    Qeq: array_like. Equilibrium discharge

    Notes
    --------
    Works with slope in both directions by applying "sign". Convention
    is based on Battjes, so negative slope gives positive discharge.

    Issues
    --------
    - Provide English reference
 
    References
    --------
    Battjes (1989) - Fluid mechanics, lecture notes Delft University (in Dutch)
    """
    counter = np.abs(ifr) * g * Ac**2 * R
    Qeq = -np.sign(ifr) * np.sqrt(counter / cf)
    return Qeq


@vectorize
def U(Q, h, z, B, mleft = 0, mright = 0, hroof=1e12):
    """
    Profile-averaged flow velocity

    Parameters
    -------
    Q: array_like. Discharge
    h: array_like. Water level with respect to reference
    z: array_like. Bottom level with respect to reference
    B: array_like. Width of profile
    mleft: array_like. Slope left bank
    mright: array_like. Slope right bank
    hroof: level of roof with respect to reference

    Returns
    -------
    U : float. Profile-averaged velocity

    Issues
    --------
    - Provide English reference
 
    References
    --------
    Battjes (1989) - Fluid mechanics, lecture notes Delft University (in Dutch)

    """
    return Q/Ac(h, z, B, mleft, mright, hroof)


@vectorize
def u(z, U, R, k):
    """
    Local horizontal flow velocity (logarithmic velocity profile)
    
    Parameters
    -------
    z : array_like. Vertical position (reference = 0)
    U : array_like. Profile averaged velocity
    R : array_like. Hydraulic radius
    k : array_like. Roughness as defined by Nikuradse

    Returns
    -------
    u : array_like. Flow velocity at given level

    Issues
    --------
    - Provide English reference
 
    References
    --------
    Battjes (1989) - Fluid mechanics, lecture notes Delft University (in Dutch)
     """
    cfloc = cf(R, k)
    z0 = k/30
    ust = ustar(cfloc, U)
    return (ust/kappa) * np.log(z/z0)


@vectorize
def c(Ac, Bs):
    """
    Wave celerity in an arbitrary profile. Equation 3.5 in Battjes 
    and Labeur

    Parameters
    -------
    Ac : array_like. Area of the flow cross section
    Bs : array_like. Width of the water surface

    Returns
    -------
    c : array_like. Wave celerity

    Issues
    --------
    - Provide English reference
 
    References
    --------
    Battjes and Labeur - Unsteady flow in open channels,
        Cambridge University Press, 2017
     """
    return np.sqrt((g * Ac)/Bs)