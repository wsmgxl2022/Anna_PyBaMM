#
# Single Particle Model (SPM)
#
import pybamm
from .base_lithium_ion_model import BaseModel


class SPM(BaseModel):
    """
    Single Particle Model (SPM) of a lithium-ion battery, from [1]_.

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
    Examples
    --------
    >>> import pybamm
    >>> model = pybamm.lithium_ion.SPM()
    >>> model.name
    'Single Particle Model'

    References
    ----------
    .. [1] SG Marquis, V Sulzer, R Timms, CP Please and SJ Chapman. “An asymptotic
           derivation of a single particle model with electrolyte”. Journal of The
           Electrochemical Society, 166(15):A3693–A3706, 2019

    **Extends:** :class:`pybamm.lithium_ion.BaseModel`
    """

    def __init__(self, options=None, name="Single Particle Model", build=True):
        # Use 'algebraic' surface form if non-default kinetics are used
        options = options or {}
        kinetics = options.get("intercalation kinetics")
        surface_form = options.get("surface form")
        if kinetics is not None and surface_form is None:
            options["surface form"] = "algebraic"

        # For degradation models we use the "x-average", note that for side reactions
        # this is set by "x-average side reactions"
        self.x_average = True

        # Set "x-average side reactions" to "true" if the model is SPM
        x_average_side_reactions = options.get("x-average side reactions")
        if x_average_side_reactions is None and self.__class__ in [
            pybamm.lithium_ion.SPM,
            pybamm.lithium_ion.MPM,
        ]:
            options["x-average side reactions"] = "true"

        super().__init__(options, name)

        self.set_submodels(build)

        if self.__class__ != "MPM":
            pybamm.citations.register("Marquis2019")

        if (
            self.options["SEI"] not in ["none", "constant"]
            or self.options["lithium plating"] != "none"
        ):
            pybamm.citations.register("BrosaPlanella2022")

    def set_convection_submodel(self):

        self.submodels[
            "through-cell convection"
        ] = pybamm.convection.through_cell.NoConvection(self.param, self.options)
        self.submodels[
            "transverse convection"
        ] = pybamm.convection.transverse.NoConvection(self.param, self.options)

    def set_intercalation_kinetics_submodel(self):

        for domain in ["negative", "positive"]:
            if self.options["surface form"] == "false":
                self.submodels[
                    f"{domain} interface"
                ] = self.inverse_intercalation_kinetics(
                    self.param, domain, "lithium-ion main", self.options
                )
                self.submodels[
                    f"{domain} interface current"
                ] = pybamm.kinetics.CurrentForInverseButlerVolmer(
                    self.param, domain, "lithium-ion main", self.options
                )
            else:
                intercalation_kinetics = self.get_intercalation_kinetics(domain)
                phases = self.options.phases[domain]
                for phase in phases:
                    submod = intercalation_kinetics(
                        self.param, domain, "lithium-ion main", self.options, phase
                    )
                    self.submodels[f"{domain} {phase} interface"] = submod
                if len(phases) > 1:
                    self.submodels[
                        f"total {domain} interface"
                    ] = pybamm.kinetics.TotalMainKinetics(
                        self.param, domain, "lithium-ion main", self.options
                    )

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

    def set_solid_submodel(self):

        self.submodels[
            "negative electrode potential"
        ] = pybamm.electrode.ohm.LeadingOrder(
            self.param, "Negative", options=self.options
        )
        self.submodels[
            "positive electrode potential"
        ] = pybamm.electrode.ohm.LeadingOrder(
            self.param, "Positive", options=self.options
        )

    def set_electrolyte_submodel(self):

        surf_form = pybamm.electrolyte_conductivity.surface_potential_form

        if self.options["electrolyte conductivity"] not in ["default", "leading order"]:
            raise pybamm.OptionError(
                "electrolyte conductivity '{}' not suitable for SPM".format(
                    self.options["electrolyte conductivity"]
                )
            )

        if self.options["surface form"] == "false" or self.half_cell:
            self.submodels[
                "leading-order electrolyte conductivity"
            ] = pybamm.electrolyte_conductivity.LeadingOrder(
                self.param, options=self.options
            )
        if self.options["surface form"] == "false":
            surf_model = surf_form.Explicit
        elif self.options["surface form"] == "differential":
            surf_model = surf_form.LeadingOrderDifferential
        elif self.options["surface form"] == "algebraic":
            surf_model = surf_form.LeadingOrderAlgebraic

        for domain in ["Negative", "Positive"]:
            self.submodels[
                domain.lower() + " surface potential difference"
            ] = surf_model(self.param, domain)

        self.submodels[
            "electrolyte diffusion"
        ] = pybamm.electrolyte_diffusion.ConstantConcentration(self.param, self.options)
