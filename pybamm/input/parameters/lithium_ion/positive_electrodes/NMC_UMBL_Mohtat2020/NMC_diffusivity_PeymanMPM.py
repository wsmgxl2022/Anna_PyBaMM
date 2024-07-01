from pybamm import exp, constants, Parameter


def NMC_diffusivity_PeymanMPM(sto, T):
    """
    NMC diffusivity as a function of stochiometry, in this case the
    diffusivity is taken to be a constant. The value is taken from Peyman MPM.

    References
    ----------
    .. [1] http://www.cchem.berkeley.edu/jsngrp/fortran.html

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stochiometry
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Solid diffusivity
    """
    # D_ref = 8 * 10 ** (-15)
    # E_D_s = 18550
    D_ref =  Parameter("Positive electrode diffusion coefficient [m2.s-1]")
    E_D_s = Parameter("Positive electrode diffusion coefficient activation energy [J.mol-1]")
    arrhenius = exp(E_D_s / constants.R * (1 / 298.15 - 1 / T))

    return D_ref * arrhenius
