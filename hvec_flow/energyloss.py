"""
Package hvec_flow, subpackage energyloss.

Prepared by HVEC lab, 2022
"""

from numpy import vectorize
import numpy as np


@vectorize
def ksi_entry(mu_contr):
    """
    Loss coefficient for entry of a culvert

    Parameters
    -------
    mu_contr : array_like. Contraction of the flow

    Returns
    -------
    ksi_entry : float
        Dimensionless loss coefficient

    Issues
    --------
    - Provide English reference

    References
    --------
    Battjes (1989) - Fluid mechanics, lecture notes Delft University (in Dutch), pg. 169
    """
    return ((1/mu_contr) - 1) ** 2



@vectorize
def cf(R, k):
    """
    Dimensionless friction factor based on the work of Colebrook.
    Based on equation 12.17 in Battjes (1989).

    Theoretical remarks:
    - We have compared the equations in Battjes (1989) and Battjes and 
    Labeur (2017), concluding that lambda in the former equals cf in
    the latter;
    - We have adopted the formula for hydraulic rough conditions only,
    envisaging an application in hydraulic engineering only

    Parameters
    -------
    R : array_like. Hydraulic radius
    k : array_like. Roughness as defined by Nikuradse

    Returns
    -------
    cf : float
        Dimensionless friction coefficient

    Issues
    --------
    - Provide English reference

    References
    --------
    Battjes (1989) - Fluid mechanics, lecture notes Delft University (in Dutch)
    Battjes and Labeur - Unsteady flow in open channels, 
        Cambridge University Press, 2017 
    """
    R = max(R, 1e-3)

    if (k>0) & (R>0):
        res = ((5.75) * np.log10((12 * R)/k))**(-2)
    else:
        res = 1e-6
        warnings.warn(
            "hvec_flow function c_f. Negative roughness or negative hydraulic radius"
            )
    return res.squeeze()


@vectorize
def ustar(cf, U):
    """
    Friction velocity
    
    Parameters
    -------
    cf : array_like. Friction factor
    U : array_like. Profile-averaged velocity

    Returns
    -------
    ustar : float. Friction velocity

    Issues
    --------
    - Provide English reference
 
    References
    --------
    Battjes (1989) - Fluid mechanics, lecture notes Delft University (in Dutch)
    """
    return np.sqrt(cf) * U


@vectorize
def Ch(R, k):
    """
    Chezy coefficient calculated from c_f.

    The Chezy coefficient is still widely used in (Dutch) practice.
    However, Chezy is not dimensionless and hence the use of the 
    friction factor c_f is preferred. 

    It is advised to use this function only to provide reference 
    information for verification purposes.

    By calculating Ch from c_f, we ensure consistency between approaches.

    Parameters
    -------
    R : array_like. Hydraulic radius
    k : array_like. Roughness as defined by Nikuradse

    Returns
    -------
    Ch : float
        Value of Chezy coefficient

    Issues
    --------
    - Provide English reference

    References
    --------
    Battjes (1989) - Fluid mechanics, lecture notes Delft University (in Dutch)
    """
    if cf(R, k)>0:
        res = np.sqrt(g/cf(R,k))
    else:
        res = float("NaN")    
    return res


@vectorize
def ksi_fr(L, R, k):
    """
    Loss coefficient for friction in a conduit.

    By definition the following holds:
    delta_H = ksi * (U**2/2*g)

    From 12.12 in Battjes (1989) then follows this equation.

    Parameters
    -------
    L: array_like. Length
    R : array_like. Hydraulic radius
    k : array_like. Roughness as defined by Nikuradse

    Returns
    -------
    ksi_fr : float
        loss coefficient

    Issues
    --------
    - Provide English reference
 
    References
    --------
    Battjes (1989) - Fluid mechanics, lecture notes Delft University (in Dutch)
    """
    return (2 * cf(R, k) * L)/R


@vectorize
def mu(ksi_tot):
    """
    Discharge coefficient calculated from summed loss coefficients

    Parameters
    -------
    ksi_tot: array_like. Total loss coefficients

    Returns
    -------
    mu : float
        Value of discharge coefficient

    Issues
    --------
     - Provide English reference
 
    References
    --------
    Nortier - Fluid mechanics, Educaboek 1989 (in Dutch)
    """
    return 1/np.sqrt(ksi_tot)


@vectorize
def ifr(cf, U, R):
    """
    Local slope of energy line due to friction
    Parameters

    -------
    cf: array_like. Friction factor
    U: array_like. Profile-averaged velocity
    R: array_like. Hydraulic radius

    Returns
    -------
    ifr : float. Slope of energy line

    Issues
    --------
     - Provide English reference
 
    References
    --------
    Battjes (1989) - Fluid mechanics, lecture notes Delft University (in Dutch)
    """
    return -(cf * U**2)/(g * R)