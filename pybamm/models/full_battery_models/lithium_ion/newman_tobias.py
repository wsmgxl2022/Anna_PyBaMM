#
# Newman Tobias Model
#
import pybamm
from .dfn import DFN


class NewmanTobias(DFN):
    """
    Newman-Tobias model of a lithium-ion battery based on the formulation in [1]_.
    This model assumes a uniform concentration profile in the electrolyte.
    Unlike the model posed in [1]_, this model accounts for nonlinear Butler-Volmer
    kinetics. It also tracks the average concentration in the solid phase in each
    electrode, which is equivalent to including an equation for the local state of
    charge as in [2]_. The user can pass the "particle" option to include mass
    transport in the particles.

    Parameters
    ----------
    options : dict, optional
        A dictionary of options to be passed to the model. For a detailed list of
        options see :class:`~pybamm.BatteryModelOptions`.
    name : str, optional
        The name of the model.
    build :  bool, optional
        Whether to build the model on instantiation. Default is True. Setting this
        option to False allows users to change any number of the submodels before
        building the complete model (submodels cannot be changed after the model is
        built).

    References
    ----------
    .. [1] JS Newman and CW Tobias. "Theoretical Analysis of Current Distribution
           in Porous Electrodes". Journal of The Electrochemical Society,
           109(12):A1183-A1191, 1962
    .. [2] HN Chu, SU Kim, SK Rahimian, JB Siegel and CW Monroe. "Parameterization
           of prismatic lithium–iron–phosphate cells through a streamlined
           thermal/electrochemical model". Journal of Power Sources, 453, p.227787,
           2020


    **Extends:** :class:`pybamm.lithium_ion.DFN`
    """

    def __init__(self, options=None, name="Newman-Tobias model", build=True):

        # Set default option "uniform profile" for particle submodel. Other
        # default options are those given in `pybamm.BatteryModelOptions` defined in
        # `base_battery_model.py`.
        options = options or {}
        if "particle" not in options:
            options["particle"] = "uniform profile"

        super().__init__(options, name, build)

        pybamm.citations.register("Newman1962")
        pybamm.citations.register("Chu2020")

    def set_particle_submodel(self):
        for domain in ["negative", "positive"]:
            particle = getattr(self.options, domain)["particle"]
            for phase in self.options.phases[domain]:
                if particle == "Fickian diffusion":
                    submod = pybamm.particle.FickianDiffusion(
                        self.param, domain, self.options, phase=phase, x_average=True
                    )
                elif particle in [
                    "uniform profile",
                    "quadratic profile",
                    "quartic profile",
                ]:
                    submod = pybamm.particle.XAveragedPolynomialProfile(
                        self.param, domain, self.options, phase=phase
                    )
                self.submodels[f"{domain} {phase} particle"] = submod

    def set_electrolyte_submodel(self):

        surf_form = pybamm.electrolyte_conductivity.surface_potential_form

        self.submodels[
            "electrolyte diffusion"
        ] = pybamm.electrolyte_diffusion.ConstantConcentration(self.param)

        if self.options["electrolyte conductivity"] not in ["default", "full"]:
            raise pybamm.OptionError(
                "electrolyte conductivity '{}' not suitable for Newman-Tobias".format(
                    self.options["electrolyte conductivity"]
                )
            )

        if self.options["surface form"] == "false":
            self.submodels[
                "electrolyte conductivity"
            ] = pybamm.electrolyte_conductivity.Full(self.param)
            surf_model = surf_form.Explicit
        elif self.options["surface form"] == "differential":
            surf_model = surf_form.FullDifferential
        elif self.options["surface form"] == "algebraic":
            surf_model = surf_form.FullAlgebraic

        for domain in ["Negative", "Separator", "Positive"]:
            self.submodels[
                domain.lower() + " surface potential difference"
            ] = surf_model(self.param, domain)
