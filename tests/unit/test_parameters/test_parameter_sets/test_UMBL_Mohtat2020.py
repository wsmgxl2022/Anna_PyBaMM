#
# Tests for LG M50 parameter set loads
#
import pybamm
import unittest


class TestMohtat(unittest.TestCase):
    def test_load_params(self):
        negative_electrode = pybamm.ParameterValues({}).read_parameters_csv(
            pybamm.get_parameters_filepath(
                "input/parameters/lithium_ion/"
                + "negative_electrodes/graphite_UMBL_Mohtat2020/parameters.csv"
            )
        )
        self.assertEqual(negative_electrode["Negative electrode porosity"], "0.3")

        positive_electrode = pybamm.ParameterValues({}).read_parameters_csv(
            pybamm.get_parameters_filepath(
                "input/parameters/lithium_ion/"
                + "positive_electrodes/NMC_UMBL_Mohtat2020/parameters.csv"
            )
        )
        self.assertEqual(positive_electrode["Positive electrode porosity"], "0.3")

        electrolyte = pybamm.ParameterValues({}).read_parameters_csv(
            pybamm.get_parameters_filepath(
                "input/parameters/lithium_ion/electrolytes/LiPF6_Mohtat2020/"
                + "parameters.csv"
            )
        )
        self.assertEqual(electrolyte["Cation transference number"], "0.38")

        cell = pybamm.ParameterValues({}).read_parameters_csv(
            pybamm.get_parameters_filepath(
                "input/parameters/lithium_ion/cells/UMBL_Mohtat2020/parameters.csv"
            )
        )
        self.assertAlmostEqual(
            cell["Negative current collector thickness [m]"], 2.5e-05
        )

    def test_standard_lithium_parameters(self):

        parameter_values = pybamm.ParameterValues("Mohtat2020")

        model = pybamm.lithium_ion.DFN()
        sim = pybamm.Simulation(model, parameter_values=parameter_values)
        sim.set_parameters()
        sim.build()


if __name__ == "__main__":
    print("Add -v for more debug output")
    import sys

    if "-v" in sys.argv:
        debug = True
    pybamm.settings.debug_mode = True
    unittest.main()
