#
# Lithium-ion base model class
#
import pybamm


class BaseModel(pybamm.BaseBatteryModel):
    """
    Overwrites default parameters from Base Model with default parameters for
    lithium-ion models

    **Extends:** :class:`pybamm.BaseBatteryModel`

    """

    def __init__(self, options=None, name="Unnamed lithium-ion model", build=False):
        super().__init__(options, name)
        self.param = pybamm.LithiumIonParameters(options)

        # Assess whether the submodel is a half-cell model
        self.half_cell = self.options["working electrode"] != "both"

        # Default timescale
        self._timescale = self.param.timescale

        # Set default length scales
        self._length_scales = {
            "negative electrode": self.param.L_x,
            "separator": self.param.L_x,
            "positive electrode": self.param.L_x,
            "positive particle": self.param.p.prim.R_typ,
            "positive primary particle": self.param.p.prim.R_typ,
            "positive particle size": self.param.p.prim.R_typ,
            "current collector y": self.param.L_z,
            "current collector z": self.param.L_z,
        }

        # Add negative particle domains only if not a half cell model
        if not self.half_cell:
            self.length_scales.update(
                {
                    "negative particle": self.param.n.prim.R_typ,
                    "negative primary particle": self.param.n.prim.R_typ,
                    "negative particle size": self.param.n.prim.R_typ,
                }
            )

        # Add relevant secondary length scales
        if len(self.options.phases["positive"]) >= 2:
            self._length_scales["positive secondary particle"] = self.param.p.sec.R_typ
        if not self.half_cell and len(self.options.phases["negative"]) >= 2:
            self._length_scales["negative secondary particle"] = self.param.n.sec.R_typ

        self.set_standard_output_variables()

    def set_submodels(self, build):
        self.set_external_circuit_submodel()
        self.set_porosity_submodel()
        self.set_interface_utilisation_submodel()
        self.set_crack_submodel()
        self.set_active_material_submodel()
        self.set_transport_efficiency_submodels()
        self.set_convection_submodel()
        self.set_open_circuit_potential_submodel()
        self.set_intercalation_kinetics_submodel()
        self.set_particle_submodel()
        self.set_solid_submodel()
        self.set_electrolyte_submodel()
        self.set_thermal_submodel()
        self.set_current_collector_submodel()

        self.set_sei_submodel()
        self.set_lithium_plating_submodel()
        self.set_total_interface_submodel()
        self.set_decomposition_submodel()

        if self.half_cell:
            # This also removes "negative electrode" submodels, so should be done last
            self.set_li_metal_counter_electrode_submodels()

        if build:
            self.build_model()

    @property
    def default_parameter_values(self):
        if self.half_cell:
            return pybamm.ParameterValues("Xu2019")
        else:
            return pybamm.ParameterValues("Marquis2019")

    @property
    def default_quick_plot_variables(self):
        if self.half_cell:
            return [
                "Electrolyte concentration [mol.m-3]",
                "Positive particle surface concentration [mol.m-3]",
                "Current [A]",
                "Electrolyte potential [V]",
                "Positive electrode potential [V]",
                "Terminal voltage [V]",
            ]
        else:
            return [
                "Negative particle surface concentration [mol.m-3]",
                "Electrolyte concentration [mol.m-3]",
                "Positive particle surface concentration [mol.m-3]",
                "Current [A]",
                "Negative electrode potential [V]",
                "Electrolyte potential [V]",
                "Positive electrode potential [V]",
                "Terminal voltage [V]",
            ]

    def set_standard_output_variables(self):
        super().set_standard_output_variables()

        # Particle concentration position
        var = pybamm.standard_spatial_vars
        self.variables.update(
            {"r_p": var.r_p, "r_p [m]": var.r_p * self.param.p.prim.R_typ}
        )
        if not self.half_cell:
            self.variables.update(
                {"r_n": var.r_n, "r_n [m]": var.r_n * self.param.n.prim.R_typ}
            )

    def set_degradation_variables(self):
        """Sets variables that quantify degradation (LAM, LLI, etc)"""
        param = self.param

        if self.half_cell:
            domains = ["positive"]
        else:
            domains = ["negative", "positive"]
        for domain in domains:
            self.variables[f"Total lithium in {domain} electrode [mol]"] = sum(
                self.variables[
                    f"Total lithium in {phase} phase in {domain} electrode [mol]"
                ]
                for phase in self.options.phases[domain]
            )

        # LAM
        if self.half_cell:
            n_Li_n = pybamm.Scalar(0)
            LAM_ne = pybamm.Scalar(0)
        else:
            C_n = self.variables["Negative electrode capacity [A.h]"]
            n_Li_n = self.variables["Total lithium in negative electrode [mol]"]
            LAM_ne = (1 - C_n / param.n.cap_init) * 100

        C_p = self.variables["Positive electrode capacity [A.h]"]

        LAM_pe = (1 - C_p / param.p.cap_init) * 100

        # LLI
        n_Li_e = self.variables["Total lithium in electrolyte [mol]"]
        n_Li_p = self.variables["Total lithium in positive electrode [mol]"]
        n_Li_particles = n_Li_n + n_Li_p
        n_Li = n_Li_particles + n_Li_e

        # LLI is usually defined based only on the percentage lithium lost from
        # particles
        LLI = (1 - n_Li_particles / param.n_Li_particles_init) * 100
        LLI_tot = (1 - n_Li / param.n_Li_init) * 100

        self.variables.update(
            {
                "LAM_ne [%]": LAM_ne,
                "LAM_pe [%]": LAM_pe,
                "LLI [%]": LLI,
                "Loss of active material in negative electrode [%]": LAM_ne,
                "Loss of active material in positive electrode [%]": LAM_pe,
                "Loss of lithium inventory [%]": LLI,
                "Loss of lithium inventory, including electrolyte [%]": LLI_tot,
                # Total lithium
                "Total lithium [mol]": n_Li,
                "Total lithium in particles [mol]": n_Li_particles,
                # Lithium lost
                "Total lithium lost [mol]": param.n_Li_init - n_Li,
                "Total lithium lost from particles [mol]": param.n_Li_particles_init
                - n_Li_particles,
                "Total lithium lost from electrolyte [mol]": param.n_Li_e_init - n_Li_e,
            }
        )

        # Lithium lost to side reactions
        # Different way of measuring LLI but should give same value
        LLI_sei = self.variables["Loss of lithium to SEI [mol]"]
        if self.half_cell:
            LLI_sei_cracks = pybamm.Scalar(0)
            LLI_pl = pybamm.Scalar(0)
        else:
            LLI_sei_cracks = self.variables["Loss of lithium to SEI on cracks [mol]"]
            LLI_pl = self.variables["Loss of lithium to lithium plating [mol]"]

        LLI_reactions = LLI_sei + LLI_sei_cracks + LLI_pl
        self.variables.update(
            {
                "Total lithium lost to side reactions [mol]": LLI_reactions,
                "Total capacity lost to side reactions [A.h]": LLI_reactions
                * param.F
                / 3600,
            }
        )

    def set_summary_variables(self):
        """
        Sets the default summary variables.
        """
        summary_variables = [
            "Time [s]",
            "Time [h]",
            "Throughput capacity [A.h]",
            "Throughput energy [W.h]",
            "Positive electrode capacity [A.h]",
            # LAM, LLI
            "Loss of active material in positive electrode [%]",
            "X-averaged positive electrode active material volume fraction",
            "Loss of lithium inventory [%]",
            "Loss of lithium inventory, including electrolyte [%]",
            # Total lithium
            "Total lithium [mol]",
            "Total lithium in electrolyte [mol]",
            "Total lithium in positive electrode [mol]",
            "Total lithium in particles [mol]",
            # Lithium lost
            "Total lithium lost [mol]",
            "Total lithium lost from particles [mol]",
            "Total lithium lost from electrolyte [mol]",
            "Loss of lithium to SEI [mol]",
            "Loss of capacity to SEI [A.h]",
            "Total lithium lost to side reactions [mol]",
            "Total capacity lost to side reactions [A.h]",
            # Resistance
            "Local ECM resistance [Ohm]",
        ]

        if not self.half_cell:
            summary_variables += [
                "Negative electrode capacity [A.h]",
                "Loss of active material in negative electrode [%]",
                "X-averaged negative electrode active material volume fraction",
                "Total lithium in negative electrode [mol]",
                "Loss of lithium to lithium plating [mol]",
                "Loss of capacity to lithium plating [A.h]",
                "X-averaged SEI thickness [m]",
                "Loss of lithium to SEI on cracks [mol]",
                "Loss of capacity to SEI on cracks [A.h]",
                "X-averaged lithium plating concentration [mol.m-3]",
                "X-averaged lithium plating thickness [m]",
            ]

        self.summary_variables = summary_variables

    def set_open_circuit_potential_submodel(self):
        for domain in ["negative", "positive"]:
            domain_options = getattr(self.options, domain)
            for phase in self.options.phases[domain]:
                ocp_option = getattr(domain_options, phase)["open circuit potential"]
                ocp_submodels = pybamm.open_circuit_potential
                if ocp_option == "single":
                    ocp_model = ocp_submodels.SingleOpenCircuitPotential
                elif ocp_option == "current sigmoid":
                    ocp_model = ocp_submodels.CurrentSigmoidOpenCircuitPotential
                self.submodels[f"{domain} {phase} open circuit potential"] = ocp_model(
                    self.param, domain, "lithium-ion main", self.options, phase
                )

    def set_sei_submodel(self):
        if self.half_cell:
            reaction_loc = "interface"
        elif self.options["x-average side reactions"] == "true":
            reaction_loc = "x-average"
        else:
            reaction_loc = "full electrode"

        phases = self.options.phases["negative"]
        for phase in phases:
            if self.options["SEI"] == "none":
                submodel = pybamm.sei.NoSEI(self.param, self.options, phase)
            elif self.options["SEI"] == "constant":
                submodel = pybamm.sei.ConstantSEI(self.param, self.options, phase)
            else:
                submodel = pybamm.sei.SEIGrowth(
                    self.param, reaction_loc, self.options, phase, cracks=False
                )
            self.submodels[f"{phase} sei"] = submodel
            # Do not set "sei on cracks" submodel for half-cells
            # For full cells, "sei on cracks" submodel must be set, even if it is zero
            if reaction_loc != "interface":
                if (
                    self.options["SEI"] in ["none", "constant"]
                    or self.options["SEI on cracks"] == "false"
                ):
                    submodel = pybamm.sei.NoSEI(
                        self.param, self.options, phase, cracks=True
                    )
                else:
                    submodel = pybamm.sei.SEIGrowth(
                        self.param, reaction_loc, self.options, phase, cracks=True
                    )
                self.submodels[f"{phase} sei on cracks"] = submodel

        if len(phases) > 1:
            self.submodels["total sei"] = pybamm.sei.TotalSEI(self.param, self.options)
            self.submodels["total sei on cracks"] = pybamm.sei.TotalSEI(
                self.param, self.options, cracks=True
            )

    def set_lithium_plating_submodel(self):
        if self.options["lithium plating"] == "none":
            self.submodels["lithium plating"] = pybamm.lithium_plating.NoPlating(
                self.param, self.options
            )
        else:
            x_average = self.options["x-average side reactions"] == "true"
            self.submodels["lithium plating"] = pybamm.lithium_plating.Plating(
                self.param, x_average, self.options
            )

    def set_total_interface_submodel(self):
        self.submodels["total interface"] = pybamm.interface.TotalInterfacialCurrent(
            self.param, "lithium-ion", self.options
        )

    def set_crack_submodel(self):
        for domain in ["Negative", "Positive"]:
            crack = getattr(self.options, domain.lower())["particle mechanics"]
            if crack == "none":
                self.submodels[
                    domain.lower() + " particle mechanics"
                ] = pybamm.particle_mechanics.NoMechanics(
                    self.param, domain, options=self.options, phase="primary"
                )
            elif crack == "swelling only":
                self.submodels[
                    domain.lower() + " particle mechanics"
                ] = pybamm.particle_mechanics.SwellingOnly(
                    self.param, domain, options=self.options, phase="primary"
                )
            elif crack == "swelling and cracking":
                self.submodels[
                    domain.lower() + " particle mechanics"
                ] = pybamm.particle_mechanics.CrackPropagation(
                    self.param,
                    domain,
                    self.x_average,
                    options=self.options,
                    phase="primary",
                )

    def set_active_material_submodel(self):
        for domain in ["negative", "positive"]:
            lam = getattr(self.options, domain)["loss of active material"]
            phases = self.options.phases[domain]
            for phase in phases:
                if lam == "none":
                    submod = pybamm.active_material.Constant(
                        self.param, domain, self.options, phase
                    )
                else:
                    submod = pybamm.active_material.LossActiveMaterial(
                        self.param, domain, self.options, self.x_average
                    )
                self.submodels[f"{domain} {phase} active material"] = submod

            # Submodel for the total active material, summing up each phase
            if len(phases) > 1:
                self.submodels[
                    f"{domain} total active material"
                ] = pybamm.active_material.Total(self.param, domain, self.options)

    def set_porosity_submodel(self):
        if (
            self.options["pressure driven porosity"] == "true"
        ):
            self.submodels["porosity"] = pybamm.porosity.PressureDriven(
                self.param, self.options, self.x_average
            )
        elif (
            self.options["SEI porosity change"] == "false"
            and self.options["lithium plating porosity change"] == "false"
        ):
            self.submodels["porosity"] = pybamm.porosity.Constant(
                self.param, self.options
            )
        elif (
            self.options["SEI porosity change"] == "true"
            or self.options["lithium plating porosity change"] == "true"
        ):
            x_average = self.options["x-average side reactions"] == "true"
            self.submodels["porosity"] = pybamm.porosity.ReactionDriven(
                self.param, self.options, x_average
            )
        

    def set_decomposition_submodel(self):
        if "true" in self.options["decomposition"]:
            self.submodels[
                "anode decomposition"
            ] = pybamm.decomposition.AnodeDecomposition(self.param)
            self.submodels[
                "cathode decomposition"
            ] = pybamm.decomposition.CathodeDecomposition(self.param)
            self.submodels["SEI decomposition"] = pybamm.decomposition.SeiDecomposition(
                self.param
            )
        else:
            self.submodels[
                "anode decomposition"
            ] = pybamm.decomposition.NoAnodeDecomposition(self.param)
            self.submodels[
                "cathode decomposition"
            ] = pybamm.decomposition.NoCathodeDecomposition(self.param)
            self.submodels[
                "SEI decomposition"
            ] = pybamm.decomposition.NoSeiDecomposition(self.param)

    def set_li_metal_counter_electrode_submodels(self):
        self.submodels[
            "counter electrode open circuit potential"
        ] = pybamm.open_circuit_potential.SingleOpenCircuitPotential(
            self.param, "Negative", "lithium metal plating", self.options, "primary"
        )

        if (
            self.options["SEI"] in ["none", "constant"]
            and self.options["intercalation kinetics"] == "symmetric Butler-Volmer"
            and self.options["surface form"] == "false"
        ):
            # only symmetric Butler-Volmer can be inverted
            self.submodels[
                "counter electrode potential"
            ] = pybamm.electrode.ohm.LithiumMetalExplicit(self.param, self.options)
            self.submodels[
                "counter electrode interface"
            ] = pybamm.kinetics.InverseButlerVolmer(
                self.param, "Negative", "lithium metal plating", self.options
            )  # assuming symmetric reaction for now so we can take the inverse
            self.submodels[
                "counter electrode interface current"
            ] = pybamm.kinetics.CurrentForInverseButlerVolmerLithiumMetal(
                self.param, "Negative", "lithium metal plating", self.options
            )
        else:
            self.submodels[
                "counter electrode potential"
            ] = pybamm.electrode.ohm.LithiumMetalSurfaceForm(self.param, self.options)
            neg_intercalation_kinetics = self.get_intercalation_kinetics("Negative")
            self.submodels["counter electrode interface"] = neg_intercalation_kinetics(
                self.param, "Negative", "lithium metal plating", self.options, "primary"
            )

        # For half-cell models, remove negative electrode submodels
        # that are not needed before building
        # We do this whether the working electrode is 'positive' or 'negative' since
        # the half-cell models are always defined assuming the positive electrode is
        # the working electrode

        # This should be done before `self.build_model`, which is the expensive part

        # Models added specifically for the counter electrode have been labelled with
        # "counter electrode" so as not to be caught by this check
        self.submodels = {
            k: v
            for k, v in self.submodels.items()
            if not (k.startswith("negative") or k == "lithium plating")
        }
