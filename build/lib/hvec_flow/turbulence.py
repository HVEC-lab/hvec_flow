def keq(cf, U):
    """
    Equilibrium turbulent energy

    Parameters
    ----------
    cf: array_like. Friction factor.
    U: array_like: Profile averaged flow velocity

    Returns
    -------
    keq: equilibrium turbulence level

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
    J.A. Battjes, Vloeistofmechanica. Lecture notes 1989
    Nezu and Nakagawa, Turbulence in open-channel flows, 1993
    """
    return c0 * cf * U**2


def nu_t(y, h, U, cf):
    """
    Vertical profile of the eddy viscosity, equation 2.59 from 
    Nezu and Nakagawa

    Parameters
    ----------
    cf: array_like. Friction factor.
    d: array_like: Water depth

    Returns
    -------
    nu_t: eddy viscosity

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
    Nezu and Nakagawa, Turbulence in open-channel flows, 1993
    """
    ksi = y/h  # Relative vertical coordinate as specified by N&N
    ust = ustar(cf, U)
    return kappa * h * ust * ksi * (1 - ksi)
