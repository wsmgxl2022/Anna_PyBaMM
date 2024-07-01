#
# Class for leading-order surface form electrolyte conductivity employing stefan-maxwell
#
import pybamm

from ..leading_order_conductivity import LeadingOrder


class BaseLeadingOrderSurfaceForm(LeadingOrder):
    """Base class for leading-order conservation of charge in the electrolyte employing
    the Stefan-Maxwell constitutive equations employing the surface potential difference
    formulation. (Leading refers to leading order in asymptotics)

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    domain : str
        The domain in which the model holds
    options : dict, optional
        A dictionary of options to be passed to the model.


    **Extends:** :class:`pybamm.electrolyte_conductivity.BaseElectrolyteConductivity`
    """

    def __init__(self, param, domain, options=None):
        super().__init__(param, domain, options)

    def get_fundamental_variables(self):

        if self.domain == "Negative":
            delta_phi_av = pybamm.standard_variables.delta_phi_n_av
        elif self.domain == "Positive":
            delta_phi_av = pybamm.standard_variables.delta_phi_p_av

        variables = self._get_standard_average_surface_potential_difference_variables(
            delta_phi_av
        )
        variables.update(
            self._get_standard_surface_potential_difference_variables(delta_phi_av)
        )
        return variables

    def get_coupled_variables(self, variables):
        # Only update coupled variables once
        if self.domain == "Negative":
            return super().get_coupled_variables(variables)
        else:
            return variables

    def set_initial_conditions(self, variables):

        delta_phi = variables[
            "X-averaged "
            + self.domain.lower()
            + " electrode surface potential difference"
        ]
        delta_phi_init = self.domain_param.prim.U_init

        self.initial_conditions = {delta_phi: delta_phi_init}

    def set_boundary_conditions(self, variables):
        if self.domain == "Negative":
            phi_e = variables["Electrolyte potential"]
            self.boundary_conditions = {
                phi_e: {
                    "left": (pybamm.Scalar(0), "Neumann"),
                    "right": (pybamm.Scalar(0), "Neumann"),
                }
            }


class LeadingOrderDifferential(BaseLeadingOrderSurfaceForm):
    """Leading-order model for conservation of charge in the electrolyte employing the
    Stefan-Maxwell constitutive equations employing the surface potential difference
    formulation and where capacitance is present.

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    options : dict, optional
        A dictionary of options to be passed to the model.

    **Extends:** :class:`BaseLeadingOrderSurfaceForm`

    """

    def __init__(self, param, domain, options=None):
        super().__init__(param, domain, options)

    def set_rhs(self, variables):
        sum_a_j = variables[
            "Sum of x-averaged "
            + self.domain.lower()
            + " electrode volumetric interfacial current densities"
        ]

        sum_a_j_av = variables[
            "X-averaged "
            + self.domain.lower()
            + " electrode total volumetric interfacial current density"
        ]
        delta_phi = variables[
            "X-averaged "
            + self.domain.lower()
            + " electrode surface potential difference"
        ]

        C_dl = self.domain_param.C_dl

        self.rhs[delta_phi] = 1 / C_dl * (sum_a_j_av - sum_a_j)


class LeadingOrderAlgebraic(BaseLeadingOrderSurfaceForm):
    """Leading-order model for conservation of charge in the electrolyte employing the
    Stefan-Maxwell constitutive equations employing the surface potential difference
    formulation.

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    options : dict, optional
        A dictionary of options to be passed to the model.


    **Extends:** :class:`BaseLeadingOrderSurfaceForm`
    """

    def __init__(self, param, domain, options=None):
        super().__init__(param, domain, options)

    def set_algebraic(self, variables):
        sum_a_j = variables[
            "Sum of x-averaged "
            + self.domain.lower()
            + " electrode volumetric interfacial current densities"
        ]

        sum_a_j_av = variables[
            "X-averaged "
            + self.domain.lower()
            + " electrode total volumetric interfacial current density"
        ]
        delta_phi = variables[
            "X-averaged "
            + self.domain.lower()
            + " electrode surface potential difference"
        ]

        self.algebraic[delta_phi] = sum_a_j_av - sum_a_j
