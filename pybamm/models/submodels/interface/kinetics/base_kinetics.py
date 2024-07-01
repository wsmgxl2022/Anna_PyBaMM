#
# Base kinetics class
#
import pybamm
from ..base_interface import BaseInterface


class BaseKinetics(BaseInterface):
    """
    Base submodel for kinetics

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

    **Extends:** :class:`pybamm.interface.BaseInterface`
    """

    def __init__(self, param, domain, reaction, options, phase="primary"):
        super().__init__(param, domain, reaction, options=options, phase=phase)

    def get_fundamental_variables(self):
        domain = self.domain.lower()
        phase_name = self.phase_name

        if (
            self.options["total interfacial current density as a state"] == "true"
            and "main" in self.reaction
        ):
            j = pybamm.Variable(
                f"Total {domain} electrode {phase_name}"
                "interfacial current density variable",
                domain=f"{domain} electrode",
                auxiliary_domains={"secondary": "current collector"},
            )

            variables = {
                f"Total {domain} electrode {phase_name}"
                "interfacial current density variable": j,
                f"Total {domain} electrode {phase_name}"
                "interfacial current density variable": j,
                f"X-averaged total {domain} electrode {phase_name}"
                "interfacial current density variable": pybamm.x_average(j),
            }
            return variables
        else:
            return {}

    def get_coupled_variables(self, variables):
        Domain = self.domain
        domain = Domain.lower()
        reaction_name = self.reaction_name
        phase_name = self.phase_name

        if self.reaction == "lithium metal plating":  # li metal electrode (half-cell)
            delta_phi = variables[
                "Lithium metal interface surface potential difference"
            ]
        else:
            delta_phi = variables[f"{Domain} electrode surface potential difference"]
            # If delta_phi was broadcast, take only the orphan.
            if isinstance(delta_phi, pybamm.Broadcast):
                delta_phi = delta_phi.orphans[0]
        # For "particle-size distribution" models, delta_phi must then be
        # broadcast to "particle size" domain
        if (
            self.reaction == "lithium-ion main"
            and self.options["particle size"] == "distribution"
        ):
            delta_phi = pybamm.PrimaryBroadcast(delta_phi, [f"{domain} particle size"])

        # Get exchange-current density
        j0 = self._get_exchange_current_density(variables)
        # Get open-circuit potential variables and reaction overpotential
        if self.options["particle size"] == "distribution":
            ocp = variables[
                f"{Domain} electrode {reaction_name}open circuit potential distribution"
            ]
        else:
            ocp = variables[f"{Domain} electrode {reaction_name}open circuit potential"]
        # If ocp was broadcast, take only the orphan.
        if isinstance(ocp, pybamm.Broadcast):
            ocp = ocp.orphans[0]
        eta_r = delta_phi - ocp

        # Get average interfacial current density
        j_tot_av, a_j_tot_av = self._get_average_total_interfacial_current_density(
            variables
        )
        # Add SEI resistance in the negative electrode
        if self.domain == "Negative":
            if self.half_cell:
                R_sei = self.phase_param.R_sei
                L_sei = variables[f"Total {phase_name}SEI thickness"]  # on interface
                eta_sei = -j_tot_av * L_sei * R_sei
            elif self.options["SEI film resistance"] == "average":
                R_sei = self.phase_param.R_sei
                L_sei_av = variables[f"X-averaged total {phase_name}SEI thickness"]
                eta_sei = -j_tot_av * L_sei_av * R_sei
            elif self.options["SEI film resistance"] == "distributed":
                R_sei = self.phase_param.R_sei
                L_sei = variables[f"Total {phase_name}SEI thickness"]
                j_tot = variables[
                    f"Total negative electrode {phase_name}"
                    "interfacial current density variable"
                ]

                # Override print_name
                j_tot.print_name = "j_tot"

                eta_sei = -j_tot * L_sei * R_sei
            else:
                eta_sei = pybamm.Scalar(0)
            eta_r += eta_sei
            if self.options["lithium plating"] != "none":
                R_plated_Li = self.phase_param.R_plated_Li
                L_plated_Li = variables[f"X-averaged {phase_name}lithium plating thickness"]
                eta_plated_Li = -j_tot_av * L_plated_Li * R_plated_Li
            else:
                eta_plated_Li = pybamm.Scalar(0)
            eta_r += eta_plated_Li

        # Get number of electrons in reaction
        ne = self._get_number_of_electrons_in_reaction()
        # Get kinetics. Note: T and u must have the same domain as j0 and eta_r
        if self.half_cell and self.domain == "Negative":
            T = variables["X-averaged cell temperature"]
            u = variables["Lithium metal interface utilisation"]
        elif j0.domain in ["current collector", ["current collector"]]:
            T = variables["X-averaged cell temperature"]
            u = variables[f"X-averaged {domain} electrode interface utilisation"]
        elif j0.domain == [f"{domain} particle size"]:
            if j0.domains["secondary"] != [f"{domain} electrode"]:
                T = variables["X-averaged cell temperature"]
                u = variables[f"X-averaged {domain} electrode interface utilisation"]
            else:
                T = variables[f"{Domain} electrode temperature"]
                u = variables[f"{Domain} electrode interface utilisation"]

            # Broadcast T onto "particle size" domain
            T = pybamm.PrimaryBroadcast(T, [f"{domain} particle size"])
        else:
            T = variables[f"{Domain} electrode temperature"]
            u = variables[f"{Domain} electrode interface utilisation"]

        # Update j, except in the "distributed SEI resistance" model, where j will be
        # found by solving an algebraic equation.
        # (In the "distributed SEI resistance" model, we have already defined j)
        j = self._get_kinetics(j0, ne, eta_r, T, u)

        if j.domain == [f"{domain} particle size"]:
            # If j depends on particle size, get size-dependent "distribution"
            # variables first
            variables.update(
                self._get_standard_size_distribution_interfacial_current_variables(j)
            )
            variables.update(
                self._get_standard_size_distribution_exchange_current_variables(j0)
            )
            variables.update(
                self._get_standard_size_distribution_overpotential_variables(eta_r)
            )

        variables.update(self._get_standard_interfacial_current_variables(j))

        variables.update(
            self._get_standard_total_interfacial_current_variables(j_tot_av, a_j_tot_av)
        )
        variables.update(self._get_standard_exchange_current_variables(j0))
        variables.update(self._get_standard_overpotential_variables(eta_r))

        variables.update(
            self._get_standard_volumetric_current_density_variables(variables)
        )

        if self.domain == "Negative" and self.reaction in [
            "lithium-ion main",
            "lithium metal plating",
            "lead-acid main",
        ]:
            variables.update(
                self._get_standard_sei_film_overpotential_variables(eta_sei)
            )

        return variables

    def set_algebraic(self, variables):
        Domain = self.domain
        domain = Domain.lower()
        phase_name = self.phase_name

        if (
            self.options["total interfacial current density as a state"] == "true"
            and "main" in self.reaction
        ):
            j_tot_var = variables[
                f"Total {domain} electrode {phase_name}"
                "interfacial current density variable"
            ]

            # Override print_name
            j_tot_var.print_name = "j_tot"

            a_j_tot = variables[
                f"Sum of {domain} electrode {phase_name}"
                "volumetric interfacial current densities"
            ]
            a = variables[
                f"{Domain} electrode {phase_name}surface area to volume ratio"
            ]

            # Algebraic equation to set the variable j_tot_var
            # equal to the sum of currents j_tot = a_j_tot / a
            self.algebraic[j_tot_var] = j_tot_var - a_j_tot / a

    def set_initial_conditions(self, variables):
        domain = self.domain.lower()
        phase_name = self.phase_name

        if (
            self.options["total interfacial current density as a state"] == "true"
            and "main" in self.reaction
        ):
            param = self.param
            j_tot_var = variables[
                f"Total {domain} electrode {phase_name}"
                "interfacial current density variable"
            ]
            current_at_0 = (
                pybamm.FunctionParameter("Current function [A]", {"Time [s]": 0})
                / param.I_typ
                * pybamm.sign(param.I_typ)
            )
            sgn = 1 if self.domain == "Negative" else -1
            # i / (a*l), assuming a=1 initially
            j_tot_av_init = sgn * current_at_0 / self.domain_param.l

            self.initial_conditions[j_tot_var] = j_tot_av_init

    def _get_interface_variables_for_first_order(self, variables):
        # This is a bit of a hack, but we need to wrap electrolyte concentration with
        # the NotConstant class
        # to differentiate it from the electrolyte concentration inside the
        # surface potential difference when taking j.diff(c_e) later on
        Domain = self.domain
        domain = Domain.lower()

        c_e_0 = pybamm.NotConstant(
            variables["Leading-order x-averaged electrolyte concentration"]
        )
        c_e = pybamm.PrimaryBroadcast(c_e_0, self.domain_for_broadcast)
        hacked_variables = {**variables, f"{Domain} electrolyte concentration": c_e}
        delta_phi = variables[
            f"Leading-order x-averaged {domain} electrode surface potential difference"
        ]
        j0 = self._get_exchange_current_density(hacked_variables)
        ne = self._get_number_of_electrons_in_reaction()
        if self.reaction == "lead-acid main":
            ocp = self.phase_param.U(c_e_0, self.param.T_init)
        elif self.reaction == "lead-acid oxygen":
            ocp = self.phase_param.U_Ox

        T = variables["X-averaged cell temperature"]
        u = variables[f"X-averaged {domain} electrode interface utilisation"]

        return c_e_0, delta_phi, j0, ne, ocp, T, u

    def _get_j_diffusion_limited_first_order(self, variables):
        """
        First-order correction to the interfacial current density due to
        diffusion-limited effects. For a general model the correction term is zero,
        since the reaction is not diffusion-limited
        """
        return pybamm.Scalar(0)
