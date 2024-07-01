#
# Full model of electrode employing Ohm's law
#
import pybamm
from .base_ohm import BaseModel


class Full(BaseModel):
    """Full model of electrode employing Ohm's law.

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    domain : str
        Either 'Negative' or 'Positive'
    options : dict, optional
        A dictionary of options to be passed to the model.

    **Extends:** :class:`pybamm.electrode.ohm.BaseModel`
    """

    def __init__(self, param, domain, options=None):
        super().__init__(param, domain, options=options)

    def get_fundamental_variables(self):

        if self.domain == "Negative":
            phi_s = pybamm.standard_variables.phi_s_n
        elif self.domain == "Positive":
            phi_s = pybamm.standard_variables.phi_s_p

        variables = self._get_standard_potential_variables(phi_s)

        return variables

    def get_coupled_variables(self, variables):

        phi_s = variables[self.domain + " electrode potential"]
        tor = variables[self.domain + " electrode transport efficiency"]
        T = variables[self.domain + " electrode temperature"]

        sigma = self.domain_param.sigma(T)

        sigma_eff = sigma * tor
        i_s = -sigma_eff * pybamm.grad(phi_s)

        variables.update({self.domain + " electrode effective conductivity": sigma_eff})

        variables.update(self._get_standard_current_variables(i_s))

        if self.domain == "Positive":
            variables.update(self._get_standard_whole_cell_variables(variables))

        return variables

    def set_algebraic(self, variables):

        phi_s = variables[self.domain + " electrode potential"]
        i_s = variables[self.domain + " electrode current density"]

        # Variable summing all of the interfacial current densities
        sum_a_j = variables[
            f"Sum of {self.domain.lower()} electrode volumetric "
            "interfacial current densities"
        ]

        self.algebraic[phi_s] = pybamm.div(i_s) + sum_a_j

    def set_boundary_conditions(self, variables):

        phi_s = variables[self.domain + " electrode potential"]
        phi_s_cn = variables["Negative current collector potential"]
        tor = variables[self.domain + " electrode transport efficiency"]
        T = variables[self.domain + " electrode temperature"]

        if self.domain == "Negative":
            lbc = (phi_s_cn, "Dirichlet")
            rbc = (pybamm.Scalar(0), "Neumann")

        elif self.domain == "Positive":
            lbc = (pybamm.Scalar(0), "Neumann")
            sigma_eff = self.param.p.sigma(T) * tor
            i_boundary_cc = variables["Current collector current density"]
            rbc = (
                i_boundary_cc / pybamm.boundary_value(-sigma_eff, "right"),
                "Neumann",
            )

        self.boundary_conditions[phi_s] = {"left": lbc, "right": rbc}

    def set_initial_conditions(self, variables):

        phi_s = variables[self.domain + " electrode potential"]

        if self.domain == "Negative":
            phi_s_init = pybamm.Scalar(0)
        elif self.domain == "Positive":
            phi_s_init = self.param.ocv_init

        self.initial_conditions[phi_s] = phi_s_init
