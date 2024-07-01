#
# Class for one-dimensional (x-direction) thermal submodel
#
import pybamm

from .base_thermal import BaseThermal


class OneDimensionalX(BaseThermal):
    """
    Class for one-dimensional (x-direction) thermal submodel.
    Note: this model assumes infinitely large electrical and thermal conductivity
    in the current collectors, so that the contribution to the Ohmic heating
    from the current collectors is zero and the boundary conditions are applied
    at the edges of the electrodes (at x=0 and x=1, in non-dimensional coordinates).
    For more information see [1]_ and [2]_.

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    options : dict, optional
        A dictionary of options to be passed to the model.

    References
    ----------
    .. [1] R Timms, SG Marquis, V Sulzer, CP Please and SJ Chapman. “Asymptotic
           Reduction of a Lithium-ion Pouch Cell Model”. SIAM Journal on Applied
           Mathematics, 81(3), 765--788, 2021
    .. [2] SG Marquis, R Timms, V Sulzer, CP Please and SJ Chapman. “A Suite of
           Reduced-Order Models of a Single-Layer Lithium-ion Pouch Cell”. Journal
           of The Electrochemical Society, 167(14):140513, 2020

    **Extends:** :class:`pybamm.thermal.BaseThermal`
    """

    def __init__(self, param, options=None):
        super().__init__(param, options=options)
        pybamm.citations.register("Timms2021")

    def get_fundamental_variables(self):
        T_n = pybamm.standard_variables.T_n
        T_s = pybamm.standard_variables.T_s
        T_p = pybamm.standard_variables.T_p
        T_cn = pybamm.BoundaryValue(T_n, "left")
        T_cp = pybamm.BoundaryValue(T_p, "right")

        T = pybamm.concatenation(T_n, T_s, T_p)
        T_x_av = self._x_average(T, T_cn, T_cp)
        T_vol_av = self._yz_average(T_x_av)

        variables = self._get_standard_fundamental_variables(
            T_cn, T_n, T_s, T_p, T_cp, T_x_av, T_vol_av
        )
        return variables

    def get_coupled_variables(self, variables):
        variables.update(self._get_standard_coupled_variables(variables))
        return variables

    def set_rhs(self, variables):
        T = variables["Cell temperature"]
        T_n = variables["Negative electrode temperature"]
        T_s = variables["Separator temperature"]
        T_p = variables["Positive electrode temperature"]

        Q = variables["Total heating"]

        # Define volumetric heat capacity
        rho_k = pybamm.concatenation(
            self.param.n.rho(T_n),
            self.param.s.rho(T_s),
            self.param.p.rho(T_p),
        )

        # Devine thermal conductivity
        lambda_k = pybamm.concatenation(
            self.param.n.lambda_(T_n),
            self.param.s.lambda_(T_s),
            self.param.p.lambda_(T_p),
        )

        # Fourier's law for heat flux
        q = -lambda_k * pybamm.grad(T)

        # N.B only y-z surface cooling is implemented for this model
        self.rhs = {
            T: (-pybamm.div(q) / self.param.delta ** 2 + self.param.B * Q)
            / (self.param.C_th * rho_k)
        }

    def set_boundary_conditions(self, variables):
        T = variables["Cell temperature"]
        T_n_left = pybamm.boundary_value(T, "left")
        T_p_right = pybamm.boundary_value(T, "right")
        T_amb = variables["Ambient temperature"]

        # N.B only y-z surface cooling is implemented for this thermal model.
        # Tab and edge cooling is not accounted for.
        self.boundary_conditions = {
            T: {
                "left": (
                    self.param.n.h_cc
                    * (T_n_left - T_amb)
                    / self.param.n.lambda_(T_n_left),
                    "Neumann",
                ),
                "right": (
                    -self.param.p.h_cc
                    * (T_p_right - T_amb)
                    / self.param.p.lambda_(T_p_right),
                    "Neumann",
                ),
            }
        }

    def set_initial_conditions(self, variables):
        T = variables["Cell temperature"]
        self.initial_conditions = {T: self.param.T_init}
