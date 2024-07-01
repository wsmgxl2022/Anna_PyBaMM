#
# Class for ohmic electrodes in the surface potential formulation
#
import pybamm

from .base_ohm import BaseModel


class SurfaceForm(BaseModel):
    """A submodel for the electrode with Ohm's law in the surface potential
    formulation.

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

    def get_coupled_variables(self, variables):

        x_n = pybamm.standard_spatial_vars.x_n
        x_p = pybamm.standard_spatial_vars.x_p
        i_boundary_cc = variables["Current collector current density"]
        i_e = variables[self.domain + " electrolyte current density"]
        tor = variables[self.domain + " electrode transport efficiency"]
        phi_s_cn = variables["Negative current collector potential"]
        T = variables[self.domain + " electrode temperature"]

        conductivity = self.domain_param.sigma(T) * tor
        i_s = i_boundary_cc - i_e

        if self.domain == "Negative":
            phi_s = phi_s_cn - pybamm.IndefiniteIntegral(i_s / conductivity, x_n)

        elif self.domain == "Positive":

            phi_e_s = variables["Separator electrolyte potential"]
            delta_phi_p = variables["Positive electrode surface potential difference"]

            phi_s = -pybamm.IndefiniteIntegral(i_s / conductivity, x_p) + (
                pybamm.boundary_value(phi_e_s, "right")
                + pybamm.boundary_value(delta_phi_p, "left")
            )

        variables.update(self._get_standard_potential_variables(phi_s))
        variables.update(self._get_standard_current_variables(i_s))

        if (
            self.half_cell or "Negative electrode current density" in variables
        ) and "Positive electrode current density" in variables:
            variables.update(self._get_standard_whole_cell_variables(variables))

        return variables
