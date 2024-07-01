#
# Tests for the lithium-ion electrode-specific SOH model
#
import pybamm
import unittest


class TestElectrodeSOH(unittest.TestCase):
    def test_known_solution(self):

        param = pybamm.LithiumIonParameters()
        parameter_values = pybamm.ParameterValues("Mohtat2020")

        esoh_solver = pybamm.lithium_ion.ElectrodeSOHSolver(parameter_values, param)

        Vmin = 3
        Vmax = 4.2
        Cn = parameter_values.evaluate(param.n.cap_init)
        Cp = parameter_values.evaluate(param.p.cap_init)
        n_Li = parameter_values.evaluate(param.n_Li_particles_init)

        inputs = {"V_max": Vmax, "V_min": Vmin, "n_Li": n_Li, "C_n": Cn, "C_p": Cp}

        # Solve the model and check outputs
        sol = esoh_solver.solve(inputs)

        self.assertAlmostEqual(sol["Up(y_100) - Un(x_100)"].data[0], Vmax, places=5)
        self.assertAlmostEqual(sol["Up(y_0) - Un(x_0)"].data[0], Vmin, places=5)
        self.assertAlmostEqual(sol["n_Li_100"].data[0], n_Li, places=5)
        self.assertAlmostEqual(sol["n_Li_0"].data[0], n_Li, places=5)

        # Solve with split esoh and check outputs
        ics = esoh_solver._set_up_solve(inputs)
        sol_split = esoh_solver._solve_split(inputs, ics)
        for key in sol.all_models[0].variables:
            self.assertAlmostEqual(sol[key].data[0], sol_split[key].data[0], places=5)

    def test_error(self):

        param = pybamm.LithiumIonParameters()
        parameter_values = pybamm.ParameterValues("Mohtat2020")

        esoh_solver = pybamm.lithium_ion.ElectrodeSOHSolver(parameter_values, param)

        Vmin = 3
        Vmax = 4.2
        Cn = parameter_values.evaluate(param.n.cap_init)
        Cp = parameter_values.evaluate(param.p.cap_init)
        n_Li = parameter_values.evaluate(param.n_Li_particles_init) * 10

        inputs = {"V_max": Vmax, "V_min": Vmin, "n_Li": n_Li, "C_n": Cn, "C_p": Cp}

        # Solve the model and check outputs
        with self.assertRaisesRegex(ValueError, "should be between 0 and 1"):
            esoh_solver.solve(inputs)


class TestElectrodeSOHHalfCell(unittest.TestCase):
    def test_known_solution(self):
        model = pybamm.lithium_ion.ElectrodeSOHHalfCell("positive")

        param = pybamm.LithiumIonParameters({"working electrode": "positive"})
        parameter_values = pybamm.ParameterValues("Xu2019")
        sim = pybamm.Simulation(model, parameter_values=parameter_values)

        V_min = 3
        V_max = 4.2
        C_w = parameter_values.evaluate(param.p.cap_init)

        # Solve the model and check outputs
        sol = sim.solve([0], inputs={"V_min": V_min, "V_max": V_max, "C_w": C_w})
        self.assertAlmostEqual(sol["Uw(x_100)"].data[0], V_max, places=5)
        self.assertAlmostEqual(sol["Uw(x_0)"].data[0], V_min, places=5)


class TestSetInitialSOC(unittest.TestCase):
    def test_known_solutions(self):

        param = pybamm.LithiumIonParameters()
        parameter_values = pybamm.ParameterValues("Mohtat2020")

        V_min = parameter_values.evaluate(param.voltage_low_cut_dimensional)
        V_max = parameter_values.evaluate(param.voltage_high_cut_dimensional)
        C_n = parameter_values.evaluate(param.n.cap_init)
        C_p = parameter_values.evaluate(param.p.cap_init)
        n_Li = parameter_values.evaluate(param.n_Li_particles_init)

        esoh_solver = pybamm.lithium_ion.ElectrodeSOHSolver(parameter_values)

        inputs = {"V_min": V_min, "V_max": V_max, "C_n": C_n, "C_p": C_p, "n_Li": n_Li}

        # Solve the model and check outputs
        esoh_sol = esoh_solver.solve(inputs)

        x, y = pybamm.lithium_ion.get_initial_stoichiometries(1, parameter_values)
        self.assertAlmostEqual(x, esoh_sol["x_100"].data[0])
        self.assertAlmostEqual(y, esoh_sol["y_100"].data[0])
        x, y = pybamm.lithium_ion.get_initial_stoichiometries(0, parameter_values)
        self.assertAlmostEqual(x, esoh_sol["x_0"].data[0])
        self.assertAlmostEqual(y, esoh_sol["y_0"].data[0])

    def test_error(self):
        with self.assertRaisesRegex(
            ValueError, "Initial SOC should be between 0 and 1"
        ):
            pybamm.lithium_ion.get_initial_stoichiometries(2, None)

        parameter_values = pybamm.ParameterValues("Chen2020")
        param = pybamm.LithiumIonParameters()
        esoh_solver = pybamm.lithium_ion.ElectrodeSOHSolver(parameter_values, param)

        C_n = parameter_values.evaluate(param.n.cap_init)
        C_p = parameter_values.evaluate(param.p.cap_init)
        n_Li = parameter_values.evaluate(param.n_Li_particles_init)

        inputs = {"V_min": 0, "V_max": 6, "C_n": C_n, "C_p": C_p, "n_Li": n_Li}
        with self.assertRaisesRegex(ValueError, "lower bound of the voltage"):
            esoh_solver._check_esoh_feasible(inputs)
        inputs = {"V_min": 3, "V_max": 6, "C_n": C_n, "C_p": C_p, "n_Li": n_Li}
        with self.assertRaisesRegex(ValueError, "upper bound of the voltage"):
            esoh_solver._check_esoh_feasible(inputs)


if __name__ == "__main__":
    print("Add -v for more debug output")
    import sys

    if "-v" in sys.argv:
        debug = True
    pybamm.settings.debug_mode = True
    unittest.main()
