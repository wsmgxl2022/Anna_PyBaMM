#
# Test for the standard lead acid parameters
#
import pybamm
from tests import get_discretisation_for_testing

import unittest


class TestStandardParametersLeadAcid(unittest.TestCase):
    def test_scipy_constants(self):
        param = pybamm.LeadAcidParameters()
        self.assertAlmostEqual(param.R.evaluate(), 8.314, places=3)
        self.assertAlmostEqual(param.F.evaluate(), 96485, places=0)

    def test_print_parameters(self):
        parameters = pybamm.LeadAcidParameters()
        parameter_values = pybamm.lead_acid.BaseModel().default_parameter_values
        output_file = "lead_acid_parameters.txt"
        parameter_values.print_parameters(parameters, output_file)
        # test print_parameters with dict and without C-rate
        del parameter_values["Nominal cell capacity [A.h]"]
        parameters = {"C_e": parameters.C_e, "sigma_n": parameters.n.sigma}
        parameter_values.print_parameters(parameters)

    def test_parameters_defaults_lead_acid(self):
        # Load parameters to be tested
        parameters = pybamm.LeadAcidParameters()
        parameter_values = pybamm.lead_acid.BaseModel().default_parameter_values
        param_eval = parameter_values.print_parameters(parameters)
        param_eval = {k: v[0] for k, v in param_eval.items()}

        # Diffusional C-rate should be smaller than C-rate
        self.assertLess(param_eval["C_e"], param_eval["C_rate"])

        # Dimensionless electrode conductivities should be large
        self.assertGreater(parameter_values.evaluate(parameters.n.sigma(0)), 10)
        self.assertGreater(parameter_values.evaluate(parameters.p.sigma(0)), 10)

        # Dimensionless oxygen exchange current density should be small
        self.assertGreater(
            1e-10, parameter_values.evaluate(parameters.p.prim.j0_Ox(1, 0))
        )

        # Rescaled dimensionless electrode conductivities should still be large
        self.assertGreater(parameter_values.evaluate(parameters.n.sigma_prime(0)), 10)
        self.assertGreater(parameter_values.evaluate(parameters.p.sigma_prime(0)), 10)
        # Dimensionless double-layer capacity should be small
        self.assertLess(param_eval["n.C_dl"], 1e-3)
        self.assertLess(param_eval["p.C_dl"], 1e-3)
        # Volume change positive in negative electrode and negative in positive
        # electrode
        self.assertLess(param_eval["n.DeltaVsurf"], 0)
        self.assertGreater(param_eval["p.DeltaVsurf"], 0)

    def test_concatenated_parameters(self):
        # create
        param = pybamm.LeadAcidParameters()
        s_param = param.s_plus_S
        self.assertIsInstance(s_param, pybamm.Concatenation)
        self.assertEqual(
            s_param.domain, ["negative electrode", "separator", "positive electrode"]
        )

        # process parameters and discretise
        parameter_values = pybamm.ParameterValues("Sulzer2019")
        disc = get_discretisation_for_testing()
        processed_s = disc.process_symbol(parameter_values.process_symbol(s_param))

        # test output
        combined_submeshes = disc.mesh.combine_submeshes(
            "negative electrode", "separator", "positive electrode"
        )
        self.assertEqual(processed_s.shape, (combined_submeshes.npts, 1))

    def test_current_functions(self):
        # create current functions
        param = pybamm.LeadAcidParameters()
        dimensional_current_density = param.dimensional_current_density_with_time
        dimensionless_current_density = param.current_with_time

        # process
        parameter_values = pybamm.ParameterValues(
            {
                "Electrode height [m]": 0.1,
                "Electrode width [m]": 0.1,
                "Negative electrode thickness [m]": 1,
                "Separator thickness [m]": 1,
                "Positive electrode thickness [m]": 1,
                "Typical electrolyte concentration [mol.m-3]": 1,
                "Number of electrodes connected in parallel to make a cell": 8,
                "Typical current [A]": 2,
                "Current function [A]": 2,
            }
        )
        dimensional_current_density_eval = parameter_values.process_symbol(
            dimensional_current_density
        )
        dimensionless_current_density_eval = parameter_values.process_symbol(
            dimensionless_current_density
        )
        self.assertAlmostEqual(
            dimensional_current_density_eval.evaluate(t=3), 2 / (8 * 0.1 * 0.1)
        )
        self.assertEqual(dimensionless_current_density_eval.evaluate(t=3), 1)

    def test_thermal_parameters(self):
        values = pybamm.lead_acid.BaseModel().default_parameter_values
        param = pybamm.LeadAcidParameters()
        T = 1  # dummy temperature as the values are constant

        # Density
        self.assertAlmostEqual(values.evaluate(param.n.rho_cc(T)), 0.8810, places=2)
        self.assertAlmostEqual(values.evaluate(param.n.rho(T)), 0.8810, places=2)
        self.assertAlmostEqual(values.evaluate(param.s.rho(T)), 0.7053, places=2)
        self.assertAlmostEqual(values.evaluate(param.p.rho(T)), 1.4393, places=2)
        self.assertAlmostEqual(values.evaluate(param.p.rho_cc(T)), 1.4393, places=2)
        self.assertAlmostEqual(values.evaluate(param.rho(T)), 1, places=2)

        # Thermal conductivity
        self.assertAlmostEqual(values.evaluate(param.n.lambda_cc(T)), 1.6963, places=2)
        self.assertAlmostEqual(values.evaluate(param.n.lambda_(T)), 1.6963, places=2)
        self.assertAlmostEqual(values.evaluate(param.s.lambda_(T)), 0.0019, places=2)
        self.assertAlmostEqual(values.evaluate(param.p.lambda_(T)), 1.6963, places=2)
        self.assertAlmostEqual(values.evaluate(param.p.lambda_cc(T)), 1.6963, places=2)

    def test_functions_lead_acid(self):
        # Load parameters to be tested
        param = pybamm.LeadAcidParameters()
        parameters = {
            "D_e_1": param.D_e(pybamm.Scalar(1), pybamm.Scalar(0)),
            "kappa_e_0": param.kappa_e(pybamm.Scalar(0), pybamm.Scalar(0)),
            "chi_1": param.chi(pybamm.Scalar(1), pybamm.Scalar(0)),
            "chi_0.5": param.chi(pybamm.Scalar(0.5), pybamm.Scalar(0)),
            "U_n_1": param.n.prim.U(pybamm.Scalar(1), pybamm.Scalar(0)),
            "U_n_0.5": param.n.prim.U(pybamm.Scalar(0.5), pybamm.Scalar(0)),
            "U_p_1": param.p.prim.U(pybamm.Scalar(1), pybamm.Scalar(0)),
            "U_p_0.5": param.p.prim.U(pybamm.Scalar(0.5), pybamm.Scalar(0)),
        }
        # Process
        parameter_values = pybamm.ParameterValues("Sulzer2019")
        param_eval = parameter_values.print_parameters(parameters)
        param_eval = {k: v[0] for k, v in param_eval.items()}

        # Known values for dimensionless functions
        self.assertEqual(param_eval["D_e_1"], 1)
        self.assertEqual(param_eval["kappa_e_0"], 0)
        # Known monotonicity for dimensionless functions
        self.assertGreater(param_eval["chi_1"], param_eval["chi_0.5"])
        self.assertLess(param_eval["U_n_1"], param_eval["U_n_0.5"])
        self.assertGreater(param_eval["U_p_1"], param_eval["U_p_0.5"])

    def test_update_initial_state_of_charge(self):
        # Load parameters to be tested
        parameters = pybamm.LeadAcidParameters()
        parameter_values = pybamm.lead_acid.BaseModel().default_parameter_values
        param_eval = parameter_values.print_parameters(parameters)
        param_eval = {k: v[0] for k, v in param_eval.items()}

        # Update initial state of charge
        parameter_values.update({"Initial State of Charge": 0.2})
        param_eval_update = parameter_values.print_parameters(parameters)
        param_eval_update = {k: v[0] for k, v in param_eval_update.items()}

        # Test that relevant parameters have changed as expected
        self.assertLess(param_eval_update["q_init"], param_eval["q_init"])
        self.assertLess(param_eval_update["c_e_init"], param_eval["c_e_init"])
        self.assertLess(
            param_eval_update["n.epsilon_init"], param_eval["n.epsilon_init"]
        )
        self.assertEqual(
            param_eval_update["s.epsilon_init"], param_eval["s.epsilon_init"]
        )
        self.assertLess(
            param_eval_update["p.epsilon_init"], param_eval["p.epsilon_init"]
        )
        self.assertGreater(
            param_eval_update["n.curlyU_init"], param_eval["n.curlyU_init"]
        )
        self.assertGreater(
            param_eval_update["p.curlyU_init"], param_eval["p.curlyU_init"]
        )


if __name__ == "__main__":
    print("Add -v for more debug output")
    import sys

    if "-v" in sys.argv:
        debug = True
    pybamm.settings.debug_mode = True
    unittest.main()
