#
# Class for swelling only (no cracking)
#
import pybamm
from .base_mechanics import BaseMechanics


class SwellingOnly(BaseMechanics):
    """
    Class for swelling only (no cracking)

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    domain : str
        The domain of the model either 'Negative' or 'Positive'
    options: dict
        A dictionary of options to be passed to the model.
        See :class:`pybamm.BaseBatteryModel`
    phase : str, optional
        Phase of the particle (default is "primary")

    **Extends:** :class:`pybamm.particle_mechanics.BaseMechanics`
    """

    def __init__(self, param, domain, options, phase="primary"):
        super().__init__(param, domain, options, phase)

    def get_fundamental_variables(self):
        zero = pybamm.FullBroadcast(
            pybamm.Scalar(0), self.domain.lower() + " electrode", "current collector"
        )
        zero_av = pybamm.x_average(zero)
        variables = self._get_standard_variables(zero)
        variables.update(
            {
                self.domain + " particle cracking rate": zero,
                "X-averaged " + self.domain + " particle cracking rate": zero_av,
            }
        )
        return variables

    def get_coupled_variables(self, variables):
        variables.update(self._get_standard_surface_variables(variables))
        variables.update(self._get_mechanical_results(variables))
        return variables
