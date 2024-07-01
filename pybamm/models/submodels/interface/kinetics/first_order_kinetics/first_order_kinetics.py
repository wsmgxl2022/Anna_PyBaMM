#
# First-order Butler-Volmer kinetics
#
import pybamm
from ...base_interface import BaseInterface


class FirstOrderKinetics(BaseInterface):
    """
    First-order kinetics

    Parameters
    ----------
    param :
        model parameters
    domain : str
        The domain to implement the model, either: 'Negative' or 'Positive'.
    leading_order_model : :class:`pybamm.interface.kinetics.BaseKinetics`
        The leading-order model with respect to which this is first-order
    options: dict
        A dictionary of options to be passed to the model. See
        :class:`pybamm.BaseBatteryModel`

    **Extends:** :class:`pybamm.interface.BaseInterface`
    """

    def __init__(self, param, domain, leading_order_model, options):
        super().__init__(param, domain, leading_order_model.reaction, options)
        self.leading_order_model = leading_order_model

    def get_coupled_variables(self, variables):
        Domain = self.domain
        domain = Domain.lower()
        reaction_name = self.reaction_name

        # Unpack
        c_e_0 = variables[f"Leading-order {domain} electrolyte concentration"]
        c_e = variables[self.domain + " electrolyte concentration"]
        c_e_1 = (c_e - c_e_0) / self.param.C_e

        dj_dc_0 = self.leading_order_model._get_dj_dc(variables)
        dj_ddeltaphi_0 = self.leading_order_model._get_dj_ddeltaphi(variables)

        # Update delta_phi with new phi_e and phi_s
        phi_s = variables[self.domain + " electrode potential"]
        phi_e = variables[self.domain + " electrolyte potential"]
        delta_phi = phi_s - phi_e
        variables.update(
            self._get_standard_average_surface_potential_difference_variables(
                pybamm.x_average(delta_phi)
            )
        )
        variables.update(
            self._get_standard_surface_potential_difference_variables(delta_phi)
        )

        delta_phi_0 = variables[
            f"Leading-order {domain} electrode surface potential difference"
        ]
        delta_phi_1 = (delta_phi - delta_phi_0) / self.param.C_e

        j_0 = variables[
            f"Leading-order {domain} electrode {reaction_name}"
            "interfacial current density"
        ]
        j_1 = dj_dc_0 * c_e_1 + dj_ddeltaphi_0 * delta_phi_1
        j = j_0 + self.param.C_e * j_1
        # Get exchange-current density
        j0 = self._get_exchange_current_density(variables)
        # Get open-circuit potential variables and reaction overpotential
        ocp = variables[f"{Domain} electrode {reaction_name}open circuit potential"]
        eta_r = delta_phi - ocp

        variables.update(self._get_standard_interfacial_current_variables(j))
        variables.update(self._get_standard_exchange_current_variables(j0))
        variables.update(self._get_standard_overpotential_variables(eta_r))

        # SEI film resistance not implemented in this model
        eta_sei = pybamm.Scalar(0)
        variables.update(self._get_standard_sei_film_overpotential_variables(eta_sei))

        # Add first-order averages
        j_1_bar = dj_dc_0 * pybamm.x_average(c_e_1) + dj_ddeltaphi_0 * pybamm.x_average(
            delta_phi_1
        )

        variables.update(
            {
                "First-order x-averaged "
                + self.domain.lower()
                + " electrode"
                + self.reaction_name
                + " interfacial current density": j_1_bar
            }
        )

        return variables
