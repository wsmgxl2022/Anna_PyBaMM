#
# Tafel classes
#
import pybamm
from .base_kinetics import BaseKinetics


class ForwardTafel(BaseKinetics):
    """
    Base submodel which implements the forward Tafel equation:

    .. math::
        j = u * j_0(c) * \\exp((ne / (2 * (1 + \\Theta T)) * \\eta_r(c))

    Parameters
    ----------
    param :
        model parameters
    domain : str
        The domain to implement the model, either: 'Negative' or 'Positive'.
    reaction : str
        The name of the reaction being implemented
    options: dict
        A dictionary of options to be passed to the model.
        See :class:`pybamm.BaseBatteryModel`
    phase : str, optional
        Phase of the particle (default is "primary")

    **Extends:** :class:`pybamm.interface.kinetics.BaseKinetics`
    """

    def __init__(self, param, domain, reaction, options, phase="primary"):
        super().__init__(param, domain, reaction, options, phase)

    def _get_kinetics(self, j0, ne, eta_r, T, u):
        alpha = self.phase_param.alpha_bv
        return (
            u * j0 * pybamm.exp((ne * alpha / (2 * (1 + self.param.Theta * T))) * eta_r)
        )

    def _get_dj_dc(self, variables):
        """See :meth:`pybamm.interface.kinetics.BaseKinetics._get_dj_dc`"""
        alpha = self.phase_param.alpha_bv
        (
            c_e,
            delta_phi,
            j0,
            ne,
            ocp,
            T,
            u,
        ) = self._get_interface_variables_for_first_order(variables)
        eta_r = delta_phi - ocp
        return (2 * u * j0.diff(c_e)) * pybamm.exp(
            (ne * alpha / (2 * (1 + self.param.Theta * T))) * eta_r
        )

    def _get_dj_ddeltaphi(self, variables):
        """See :meth:`pybamm.interface.kinetics.BaseKinetics._get_dj_ddeltaphi`"""
        alpha = self.phase_param.alpha_bv
        _, delta_phi, j0, ne, ocp, T, u = self._get_interface_variables_for_first_order(
            variables
        )
        eta_r = delta_phi - ocp
        return (2 * u * j0 * (ne / (2 * (1 + self.param.Theta * T)))) * pybamm.exp(
            (ne * alpha / (2 * (1 + self.param.Theta * T))) * eta_r
        )


# backwardtafel not used by any of the models
# class BackwardTafel(BaseKinetics):
#     """
#     Base submodel which implements the backward Tafel equation:

#     .. math::
#         j = -j_0(c) * \\exp(-\\eta_r(c))

#     Parameters
#     ----------
#     param :
#         model parameters
#     domain : str
#         The domain to implement the model, either: 'Negative' or 'Positive'.


#     **Extends:** :class:`pybamm.interface.kinetics.BaseKinetics`
#     """

#     def __init__(self, param, domain, reaction, options):
#         super().__init__(param, domain, reaction, options)

#     def _get_kinetics(self, j0, ne, eta_r, T):
#         return -j0 * pybamm.exp(-(ne / (2 * (1 + self.param.Theta * T))) * eta_r)
