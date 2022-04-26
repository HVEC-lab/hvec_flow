"""
Package hvec_flow, subpackage sills.

Prepared by HVEC lab, 2022
"""


from numpy import vectorize


@vectorize
def q_crit(H_up, h_sill):
    """
    Critical discharge per meter running length for given water level and sill level

    Parameters
    ----------
    H_up : array_like. Upstream water level with respect to reference
    h_sill : array_like. Sill level with the same reference and units as h_up

    Returns
    -------
    q_crit : float
        Critical specific discharge over a sill

    References
    --------
    Cruise, Sherif and Singh - Elementary hydraulics, 2007
    """
    return (2/3)*np.sqrt(2/3*g) * (H_up-h_sill)**(3/2)


@vectorize
def d_crit(q):
    """
    Critical depth for given specific discharge

    Parameters
    ----------
    q : array_like. Specific discharge (= discharge per unit running length)

    Returns
    -------
    d_crit : float
        Critical depth

    References
    --------
    Cruise, Sherif and Singh - Elementary hydraulics, 2007
    """
    return (3/2)*((q**2)/g)**(1/3)


@vectorize
def Q_crit(H_up, h_sill, B):
    """
    Critical total discharge for given water level and sill level

    Parameters
    ----------
    H_up : array_like. Upstream energy level with respect to reference
    h_sill : array_like. Sill level with the same reference and units as h_up
    B: array_like. Width of channel in the same units as the levels

    Returns
    -------
    Q_crit : float
        Critical total discharge

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
    NA

    References
    --------
    Cruise, Sherif and Singh - Elementary hydraulics, 2007
    """

    return q_crit(H_up, h_sill)*B


@vectorize
def h_sill_crit(h_up, q):
    """
    Critical sill level, with respect to reference, for given water level
    and discharge per meter running length

    Parameters
    -------
    h_up : array_like. Upstream water level with respect to reference
    q : array_like. Specific discharge (= discharge per unit running length)

    
    Returns
    -------
    h_sill_crit : float
        Critical total discharge

    References
    --------
    Cruise, Sherif and Singh - Elementary hydraulics, 2007
    """

    return h_up - d_crit(q)


@vectorize
def h_w_crit(h_sill, q):
    """
    Critical water level based on given sill level and specific discharge

    Parameters
    -------
    h_sill : array_like. Upstream water level with respect to reference
    q : array_like. Specific discharge (= discharge per unit running length)

    Returns
    -------
    h_w_crit : float
        Water level assuming critical discharge

    References
    --------
    Cruise, Sherif and Singh - Elementary hydraulics, 2007
    """
    return h_sill + d_crit(q)


@vectorize
def hw_from_HandQ_elem(Q, H, hbot, Bbot, mleft=0, mright=0):
    """
    Local water level for given energy level and discharge

    Parameters
    -------
    Q: array_like. Discharge
    H: array_like. Local energy level
    hbot: array_like. Bottom level with respect to reference
    Bbot: array_like. Bottom width
    mleft: array_like. Slope left bank
    mright: array_like. Slope right bank

    Returns
    -------
    hw_from_HandQ: float. Profile-averaged velocity

    References
    --------
    Cruise, Sherif and Singh - Elementary hydraulics, 2007
    """
    def f(hw, Q):  # Function to be solved
        U = Q/Ac(hw, hbot, Bbot, mleft, mright)
        Hkin = U**2/(2 * g)
        return H - hw - Hkin

    # Check discharge below critical discharge
    Qcr = Q_crit(H_up = H, h_sill = hbot, B = Bbot)
    if (Q > Qcr):
        raise ValueError("Discharge not possible at this energy level")

    # Calculate water level corresponding to discharge and energy level
    hw = root_scalar(f, args = (Q), x0 = H - hbot, bracket = [1e-6, H - hbot]).root

    return hw


@vectorize
def U_from_HandQ(Q, H, hbot, Bbot, mleft=0, mright=0):
    """
    Profile averaged velocity for given energy level and local profile

    Parameters
    -------
    H: array_like. Local energy level
    hbot: array_like. Bottom level with respect to reference
    Bbot: array_like. Bottom width
    mleft: array_like. Slope left bank
    mright: array_like. Slope right bank

    Returns
    -------
    U_from_H: float. Profile-averaged velocity

    Raises
    --------
    - Discharge higher than critical discharge for the 
        given energy level raises an error

    Issues
    --------
    - Not finished and not checked
    - Implemented for rectangular channel. Generalisation to trapezoidal
        channel is plan.

    References
    --------
    Cruise, Sherif and Singh - Elementary hydraulics, 2007
    """
    hw = hw_from_HandQ(Q, H, hbot, Bbot, mleft, mright)

    # Calculate velocity
    res = Q/Ac(hw, hbot, Bbot, mleft, mright)
    return res