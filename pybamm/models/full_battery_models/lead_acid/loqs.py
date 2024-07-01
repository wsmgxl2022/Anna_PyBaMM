#
# Lead-acid LOQS model
#
import pybamm
from .base_lead_acid_model import BaseModel


class LOQS(BaseModel):
    """
    Leading-Order Quasi-Static model for lead-acid, from [1]_.

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
    .. [1] V Sulzer, SJ Chapman, CP Please, DA Howey, and CW Monroe. Faster lead-acid
           battery simulations from porous-electrode theory: Part II. Asymptotic
           analysis. Journal of The Electrochemical Society 166.12 (2019), A2372–A2382.


    **Extends:** :class:`pybamm.lead_acid.BaseModel`
    """

    def __init__(self, options=None, name="LOQS model", build=True):
        super().__init__(options, name)

        self.set_external_circuit_submodel()
        self.set_open_circuit_potential_submodel()
        self.set_intercalation_kinetics_submodel()
        self.set_interface_utilisation_submodel()
        self.set_convection_submodel()
        self.set_porosity_submodel()
        self.set_active_material_submodel()
        self.set_transport_efficiency_submodels()
        self.set_electrolyte_submodel()
        self.set_electrode_submodels()
        self.set_thermal_submodel()
        self.set_side_reaction_submodels()
        self.set_current_collector_submodel()
        self.set_sei_submodel()
        self.set_lithium_plating_submodel()
        self.set_total_interface_submodel()

        if build:
            self.build_model()

        if self.options["dimensionality"] == 0:
            self.use_jacobian = False

        pybamm.citations.register("Sulzer2019asymptotic")

    def set_external_circuit_submodel(self):
        """
        Define how the external circuit defines the boundary conditions for the model,
        e.g. (not necessarily constant-) current, voltage, etc
        """
        if self.options["operating mode"] == "current":
            self.submodels[
                "leading order external circuit"
            ] = pybamm.external_circuit.LeadingOrderExplicitCurrentControl(
                self.param, self.options
            )
        elif self.options["operating mode"] == "voltage":
            self.submodels[
                "leading order external circuit"
            ] = pybamm.external_circuit.LeadingOrderVoltageFunctionControl(
                self.param, self.options
            )
        elif self.options["operating mode"] == "power":
            self.submodels[
                "leading order external circuit"
            ] = pybamm.external_circuit.LeadingOrderPowerFunctionControl(
                self.param, self.options
            )
        elif callable(self.options["operating mode"]):
            self.submodels[
                "leading order external circuit"
            ] = pybamm.external_circuit.LeadingOrderFunctionControl(
                self.param, self.options["operating mode"], self.options
            )

    def set_current_collector_submodel(self):

        if self.options["current collector"] in [
            "uniform",
            "potential pair quite conductive",
        ]:
            submodel = pybamm.current_collector.Uniform(self.param)
        elif self.options["current collector"] == "potential pair":
            if self.options["dimensionality"] == 1:
                submodel = pybamm.current_collector.PotentialPair1plus1D(self.param)
            elif self.options["dimensionality"] == 2:
                submodel = pybamm.current_collector.PotentialPair2plus1D(self.param)
        self.submodels["leading-order current collector"] = submodel

    def set_porosity_submodel(self):

        self.submodels["leading-order porosity"] = pybamm.porosity.ReactionDrivenODE(
            self.param, self.options, True
        )

    def set_transport_efficiency_submodels(self):
        self.submodels[
            "leading-order electrolyte transport efficiency"
        ] = pybamm.transport_efficiency.Bruggeman(self.param, "Electrolyte")
        self.submodels[
            "leading-order electrode transport efficiency"
        ] = pybamm.transport_efficiency.Bruggeman(self.param, "Electrode")

    def set_convection_submodel(self):

        if self.options["convection"] == "none":
            self.submodels[
                "leading-order transverse convection"
            ] = pybamm.convection.transverse.NoConvection(self.param)
            self.submodels[
                "leading-order through-cell convection"
            ] = pybamm.convection.through_cell.NoConvection(self.param)
        else:
            if self.options["convection"] == "uniform transverse":
                self.submodels[
                    "leading-order transverse convection"
                ] = pybamm.convection.transverse.Uniform(self.param)
            elif self.options["convection"] == "full transverse":
                self.submodels[
                    "leading-order transverse convection"
                ] = pybamm.convection.transverse.Full(self.param)
            self.submodels[
                "leading-order through-cell convection"
            ] = pybamm.convection.through_cell.Explicit(self.param)

    def set_intercalation_kinetics_submodel(self):

        if self.options["surface form"] == "false":
            self.submodels[
                "leading-order negative interface"
            ] = pybamm.kinetics.InverseButlerVolmer(
                self.param, "Negative", "lead-acid main", self.options
            )
            self.submodels[
                "leading-order positive interface"
            ] = pybamm.kinetics.InverseButlerVolmer(
                self.param, "Positive", "lead-acid main", self.options
            )
            self.submodels[
                "negative interface current"
            ] = pybamm.kinetics.CurrentForInverseButlerVolmer(
                self.param, "Negative", "lead-acid main"
            )
            self.submodels[
                "positive interface current"
            ] = pybamm.kinetics.CurrentForInverseButlerVolmer(
                self.param, "Positive", "lead-acid main"
            )
        else:
            self.submodels[
                "leading-order negative interface"
            ] = pybamm.kinetics.SymmetricButlerVolmer(
                self.param, "Negative", "lead-acid main", self.options, "primary"
            )

            self.submodels[
                "leading-order positive interface"
            ] = pybamm.kinetics.SymmetricButlerVolmer(
                self.param, "Positive", "lead-acid main", self.options, "primary"
            )
        # always use forward Butler-Volmer for the reaction submodel to be passed to the
        # higher order model
        self.reaction_submodels = {
            "Negative": [
                pybamm.kinetics.SymmetricButlerVolmer(
                    self.param, "Negative", "lead-acid main", self.options, "primary"
                )
            ],
            "Positive": [
                pybamm.kinetics.SymmetricButlerVolmer(
                    self.param, "Positive", "lead-acid main", self.options, "primary"
                )
            ],
        }

    def set_electrode_submodels(self):

        self.submodels[
            "leading-order negative electrode potential"
        ] = pybamm.electrode.ohm.LeadingOrder(self.param, "Negative")
        self.submodels[
            "leading-order positive electrode potential"
        ] = pybamm.electrode.ohm.LeadingOrder(self.param, "Positive")

    def set_electrolyte_submodel(self):

        surf_form = pybamm.electrolyte_conductivity.surface_potential_form

        if self.options["surface form"] == "false":
            self.submodels[
                "leading-order electrolyte conductivity"
            ] = pybamm.electrolyte_conductivity.LeadingOrder(self.param)
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
        ] = pybamm.electrolyte_diffusion.LeadingOrder(self.param)

    def set_side_reaction_submodels(self):
        if self.options["hydrolysis"] == "true":
            self.submodels[
                "leading-order oxygen diffusion"
            ] = pybamm.oxygen_diffusion.LeadingOrder(self.param)
            self.submodels[
                "leading-order positive oxygen interface"
            ] = pybamm.kinetics.ForwardTafel(
                self.param, "Positive", "lead-acid oxygen", self.options, "primary"
            )
            self.submodels[
                "leading-order negative oxygen interface"
            ] = pybamm.kinetics.DiffusionLimited(
                self.param,
                "Negative",
                "lead-acid oxygen",
                self.options,
                order="leading",
            )
        else:
            self.submodels[
                "leading-order oxygen diffusion"
            ] = pybamm.oxygen_diffusion.NoOxygen(self.param)
            self.submodels[
                "leading-order negative oxygen interface"
            ] = pybamm.kinetics.NoReaction(
                self.param, "Negative", "lead-acid oxygen", "primary"
            )
            self.submodels[
                "leading-order positive oxygen interface"
            ] = pybamm.kinetics.NoReaction(
                self.param, "Positive", "lead-acid oxygen", "primary"
            )
        self.reaction_submodels["Negative"].append(
            self.submodels["leading-order negative oxygen interface"]
        )
        self.reaction_submodels["Positive"].append(
            self.submodels["leading-order positive oxygen interface"]
        )
