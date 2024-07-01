#
# Tests for the Casadi Solver class
#
import pybamm
import unittest
import numpy as np
from tests import get_mesh_for_testing, get_discretisation_for_testing
from scipy.sparse import eye


class TestCasadiSolver(unittest.TestCase):
    def test_bad_mode(self):
        with self.assertRaisesRegex(ValueError, "invalid mode"):
            pybamm.CasadiSolver(mode="bad mode")

    def test_model_solver(self):
        # Create model
        model = pybamm.BaseModel()
        var = pybamm.Variable("var")
        model.rhs = {var: 0.1 * var}
        model.initial_conditions = {var: 1}
        # No need to set parameters; can use base discretisation (no spatial operators)

        # create discretisation
        disc = pybamm.Discretisation()
        model_disc = disc.process_model(model, inplace=False)
        # Solve
        solver = pybamm.CasadiSolver(mode="fast", rtol=1e-8, atol=1e-8)
        t_eval = np.linspace(0, 1, 100)
        solution = solver.solve(model_disc, t_eval)
        np.testing.assert_array_equal(solution.t, t_eval)
        np.testing.assert_array_almost_equal(
            solution.y.full()[0], np.exp(0.1 * solution.t), decimal=5
        )

        # Safe mode (enforce events that won't be triggered)
        model.events = [pybamm.Event("an event", var + 1)]
        disc.process_model(model)
        solver = pybamm.CasadiSolver(rtol=1e-8, atol=1e-8)
        solution = solver.solve(model, t_eval)
        np.testing.assert_array_equal(solution.t, t_eval)
        np.testing.assert_array_almost_equal(
            solution.y.full()[0], np.exp(0.1 * solution.t), decimal=5
        )

        # Safe mode, without grid (enforce events that won't be triggered)
        solver = pybamm.CasadiSolver(mode="safe without grid", rtol=1e-8, atol=1e-8)
        solution = solver.solve(model, t_eval)
        np.testing.assert_array_equal(solution.t, t_eval)
        np.testing.assert_array_almost_equal(
            solution.y.full()[0], np.exp(0.1 * solution.t), decimal=5
        )

        # Fast with events
        # with an ODE model this behaves exactly the same as "fast"
        solver = pybamm.CasadiSolver(mode="fast with events", rtol=1e-8, atol=1e-8)
        solution = solver.solve(model, t_eval)
        np.testing.assert_array_equal(solution.t, t_eval)
        np.testing.assert_array_almost_equal(
            solution.y.full()[0], np.exp(0.1 * solution.t), decimal=5
        )

    def test_model_solver_python(self):
        # Create model
        pybamm.set_logging_level("ERROR")
        model = pybamm.BaseModel()
        model.convert_to_format = "python"
        var = pybamm.Variable("var")
        model.rhs = {var: 0.1 * var}
        model.initial_conditions = {var: 1}
        # No need to set parameters; can use base discretisation (no spatial operators)

        # create discretisation
        disc = pybamm.Discretisation()
        disc.process_model(model)
        # Solve
        solver = pybamm.CasadiSolver(mode="fast", rtol=1e-8, atol=1e-8)
        t_eval = np.linspace(0, 1, 100)
        solution = solver.solve(model, t_eval)
        np.testing.assert_array_equal(solution.t, t_eval)
        np.testing.assert_array_almost_equal(
            solution.y.full()[0], np.exp(0.1 * solution.t), decimal=5
        )
        pybamm.set_logging_level("WARNING")

    def test_model_solver_failure(self):
        # Create model
        model = pybamm.BaseModel()
        var = pybamm.Variable("var")
        model.rhs = {var: -pybamm.sqrt(var)}
        model.initial_conditions = {var: 1}
        # add events so that safe mode is used (won't be triggered)
        model.events = [pybamm.Event("10", var - 10)]
        # No need to set parameters; can use base discretisation (no spatial operators)

        # create discretisation
        disc = pybamm.Discretisation()
        model_disc = disc.process_model(model, inplace=False)

        solver = pybamm.CasadiSolver(extra_options_call={"regularity_check": False})
        # Solve with failure at t=2
        t_eval = np.linspace(0, 20, 100)
        with self.assertRaises(pybamm.SolverError):
            solver.solve(model_disc, t_eval)
        # Solve with failure at t=0
        model.initial_conditions = {var: 0}
        model_disc = disc.process_model(model, inplace=False)
        t_eval = np.linspace(0, 20, 100)
        with self.assertRaises(pybamm.SolverError):
            solver.solve(model_disc, t_eval)

    def test_model_solver_events(self):
        # Create model
        model = pybamm.BaseModel()
        whole_cell = ["negative electrode", "separator", "positive electrode"]
        var1 = pybamm.Variable("var1", domain=whole_cell)
        var2 = pybamm.Variable("var2", domain=whole_cell)
        model.rhs = {var1: 0.1 * var1}
        model.algebraic = {var2: 2 * var1 - var2}
        model.initial_conditions = {var1: 1, var2: 2}
        model.events = [
            pybamm.Event("var1 = 1.5", pybamm.min(1.5 - var1)),
            pybamm.Event("var2 = 2.5", pybamm.min(2.5 - var2)),
        ]
        disc = get_discretisation_for_testing()
        disc.process_model(model)

        # Solve using "safe" mode
        solver = pybamm.CasadiSolver(mode="safe", rtol=1e-8, atol=1e-8)
        t_eval = np.linspace(0, 5, 100)
        solution = solver.solve(model, t_eval)
        np.testing.assert_array_less(solution.y.full()[0, :-1], 1.5)
        np.testing.assert_array_less(solution.y.full()[-1, :-1], 2.5)
        np.testing.assert_equal(solution.t_event[0], solution.t[-1])
        np.testing.assert_array_equal(solution.y_event[:, 0], solution.y.full()[:, -1])
        np.testing.assert_array_almost_equal(
            solution.y.full()[0], np.exp(0.1 * solution.t), decimal=5
        )
        np.testing.assert_array_almost_equal(
            solution.y.full()[-1], 2 * np.exp(0.1 * solution.t), decimal=5
        )

        # Solve using "safe" mode with debug off
        pybamm.settings.debug_mode = False
        solver = pybamm.CasadiSolver(mode="safe", rtol=1e-8, atol=1e-8, dt_max=1)
        t_eval = np.linspace(0, 5, 100)
        solution = solver.solve(model, t_eval)
        np.testing.assert_array_less(solution.y.full()[0], 1.5)
        np.testing.assert_array_less(solution.y.full()[-1], 2.5 + 1e-10)
        # test the last entry is exactly 2.5
        np.testing.assert_array_almost_equal(solution.y[-1, -1], 2.5, decimal=2)
        np.testing.assert_array_almost_equal(
            solution.y.full()[0], np.exp(0.1 * solution.t), decimal=5
        )
        np.testing.assert_array_almost_equal(
            solution.y.full()[-1], 2 * np.exp(0.1 * solution.t), decimal=5
        )
        pybamm.settings.debug_mode = True

        # Try dt_max=0 to enforce using all timesteps
        solver = pybamm.CasadiSolver(dt_max=0, rtol=1e-8, atol=1e-8)
        t_eval = np.linspace(0, 5, 100)
        solution = solver.solve(model, t_eval)
        np.testing.assert_array_less(solution.y.full()[0], 1.5)
        np.testing.assert_array_less(solution.y.full()[-1], 2.5 + 1e-10)
        np.testing.assert_array_almost_equal(
            solution.y.full()[0], np.exp(0.1 * solution.t), decimal=5
        )
        np.testing.assert_array_almost_equal(
            solution.y.full()[-1], 2 * np.exp(0.1 * solution.t), decimal=5
        )

        # Solve using "fast with events" mode
        model = pybamm.BaseModel()
        var1 = pybamm.Variable("var1")
        var2 = pybamm.Variable("var2")
        model.rhs = {var1: 0.1 * var1}
        model.algebraic = {var2: 2 * var1 - var2}
        model.initial_conditions = {var1: 1, var2: 2}
        model.events = [
            pybamm.Event("var1 = 1.5", 1.5 - var1),
            pybamm.Event("var2 = 2.5", 2.5 - var2),
            pybamm.Event("var1 = 1.5 switch", var1 - 2, pybamm.EventType.SWITCH),
            pybamm.Event("var2 = 2.5 switch", var2 - 3, pybamm.EventType.SWITCH),
        ]

        solver = pybamm.CasadiSolver(mode="fast with events", rtol=1e-8, atol=1e-8)
        t_eval = np.linspace(0, 5, 100)
        solution = solver.solve(model, t_eval)
        np.testing.assert_array_less(solution.y.full()[0, :-1], 1.5)
        np.testing.assert_array_less(solution.y.full()[-1, :-1], 2.5)
        np.testing.assert_equal(solution.t_event[0], solution.t[-1])
        np.testing.assert_array_almost_equal(
            solution.y_event[:, 0].flatten(), [1.25, 2.5], decimal=5
        )
        np.testing.assert_array_almost_equal(
            solution.y.full()[0], np.exp(0.1 * solution.t), decimal=5
        )
        np.testing.assert_array_almost_equal(
            solution.y.full()[-1], 2 * np.exp(0.1 * solution.t), decimal=5
        )

        # Test when an event returns nan
        model = pybamm.BaseModel()
        var = pybamm.Variable("var")
        model.rhs = {var: 0.1 * var}
        model.initial_conditions = {var: 1}
        model.events = [
            pybamm.Event("event", 1.02 - var),
            pybamm.Event("sqrt event", pybamm.sqrt(1.0199 - var)),
        ]
        disc = pybamm.Discretisation()
        disc.process_model(model)
        solver = pybamm.CasadiSolver(rtol=1e-8, atol=1e-8)
        solution = solver.solve(model, t_eval)
        np.testing.assert_array_less(solution.y.full()[0], 1.02 + 1e-10)
        np.testing.assert_array_almost_equal(solution.y[0, -1], 1.02, decimal=2)

    def test_model_step(self):
        # Create model
        model = pybamm.BaseModel()
        domain = ["negative electrode", "separator", "positive electrode"]
        var = pybamm.Variable("var", domain=domain)
        model.rhs = {var: 0.1 * var}
        model.initial_conditions = {var: 1}
        # No need to set parameters; can use base discretisation (no spatial operators)

        # create discretisation
        mesh = get_mesh_for_testing()
        spatial_methods = {"macroscale": pybamm.FiniteVolume()}
        disc = pybamm.Discretisation(mesh, spatial_methods)
        disc.process_model(model)

        solver = pybamm.CasadiSolver(rtol=1e-8, atol=1e-8)

        # Step once
        dt = 1
        step_sol = solver.step(None, model, dt)
        np.testing.assert_array_equal(step_sol.t, [0, dt])
        np.testing.assert_array_almost_equal(
            step_sol.y.full()[0], np.exp(0.1 * step_sol.t)
        )

        # Step again (return 5 points)
        step_sol_2 = solver.step(step_sol, model, dt, npts=5)
        np.testing.assert_array_equal(
            step_sol_2.t, np.concatenate([np.array([0]), np.linspace(dt, 2 * dt, 5)])
        )
        np.testing.assert_array_almost_equal(
            step_sol_2.y.full()[0], np.exp(0.1 * step_sol_2.t)
        )

        # Check steps give same solution as solve
        t_eval = step_sol.t
        solution = solver.solve(model, t_eval)
        np.testing.assert_array_almost_equal(solution.y.full()[0], step_sol.y.full()[0])

    def test_model_step_with_input(self):
        # Create model
        model = pybamm.BaseModel()
        var = pybamm.Variable("var")
        a = pybamm.InputParameter("a")
        model.rhs = {var: a * var}
        model.initial_conditions = {var: 1}
        model.variables = {"a": a}
        # No need to set parameters; can use base discretisation (no spatial operators)

        # create discretisation
        disc = pybamm.Discretisation()
        disc.process_model(model)

        solver = pybamm.CasadiSolver(rtol=1e-8, atol=1e-8)

        # Step with an input
        dt = 0.1
        step_sol = solver.step(None, model, dt, npts=5, inputs={"a": 0.1})
        np.testing.assert_array_equal(step_sol.t, np.linspace(0, dt, 5))
        np.testing.assert_allclose(step_sol.y.full()[0], np.exp(0.1 * step_sol.t))

        # Step again with different inputs
        step_sol_2 = solver.step(step_sol, model, dt, npts=5, inputs={"a": -1})
        np.testing.assert_array_equal(step_sol_2.t, np.linspace(0, 2 * dt, 9))
        np.testing.assert_array_equal(
            step_sol_2["a"].entries, np.array([0.1, 0.1, 0.1, 0.1, 0.1, -1, -1, -1, -1])
        )
        np.testing.assert_allclose(
            step_sol_2.y.full()[0],
            np.concatenate(
                [
                    np.exp(0.1 * step_sol_2.t[:5]),
                    np.exp(0.1 * step_sol_2.t[4]) * np.exp(-(step_sol_2.t[5:] - dt)),
                ]
            ),
        )

    def test_model_step_events(self):
        # Create model
        model = pybamm.BaseModel()
        var1 = pybamm.Variable("var1")
        var2 = pybamm.Variable("var2")
        model.rhs = {var1: 0.1 * var1}
        model.algebraic = {var2: 2 * var1 - var2}
        model.initial_conditions = {var1: 1, var2: 2}
        model.events = [
            pybamm.Event("var1 = 1.5", pybamm.min(1.5 - var1)),
            pybamm.Event("var2 = 2.5", pybamm.min(2.5 - var2)),
        ]
        disc = pybamm.Discretisation()
        disc.process_model(model)

        # Solve
        step_solver = pybamm.CasadiSolver(rtol=1e-8, atol=1e-8)
        dt = 0.05
        time = 0
        end_time = 5
        step_solution = None
        while time < end_time:
            step_solution = step_solver.step(step_solution, model, dt=dt, npts=10)
            time += dt
        np.testing.assert_array_less(step_solution.y.full()[0, :-1], 1.5)
        np.testing.assert_array_less(step_solution.y.full()[-1, :-1], 2.5)
        np.testing.assert_equal(step_solution.t_event[0], step_solution.t[-1])
        np.testing.assert_array_equal(
            step_solution.y_event[:, 0], step_solution.y.full()[:, -1]
        )
        np.testing.assert_array_almost_equal(
            step_solution.y.full()[0], np.exp(0.1 * step_solution.t), decimal=5
        )
        np.testing.assert_array_almost_equal(
            step_solution.y.full()[-1], 2 * np.exp(0.1 * step_solution.t), decimal=4
        )

    def test_model_solver_with_inputs(self):
        # Create model
        model = pybamm.BaseModel()
        domain = ["negative electrode", "separator", "positive electrode"]
        var = pybamm.Variable("var", domain=domain)
        model.rhs = {var: -pybamm.InputParameter("rate") * var}
        model.initial_conditions = {var: 1}
        model.events = [pybamm.Event("var=0.5", pybamm.min(var - 0.5))]
        # No need to set parameters; can use base discretisation (no spatial
        # operators)

        # create discretisation
        mesh = get_mesh_for_testing()
        spatial_methods = {"macroscale": pybamm.FiniteVolume()}
        disc = pybamm.Discretisation(mesh, spatial_methods)
        disc.process_model(model)
        # Solve
        solver = pybamm.CasadiSolver(rtol=1e-8, atol=1e-8)
        t_eval = np.linspace(0, 10, 100)
        solution = solver.solve(model, t_eval, inputs={"rate": 0.1})
        self.assertLess(len(solution.t), len(t_eval))
        np.testing.assert_allclose(
            solution.y.full()[0], np.exp(-0.1 * solution.t), rtol=1e-04
        )

        # Without grid
        solver = pybamm.CasadiSolver(mode="safe without grid", rtol=1e-8, atol=1e-8)
        t_eval = np.linspace(0, 10, 100)
        solution = solver.solve(model, t_eval, inputs={"rate": 0.1})
        self.assertLess(len(solution.t), len(t_eval))
        np.testing.assert_allclose(
            solution.y.full()[0], np.exp(-0.1 * solution.t), rtol=1e-04
        )
        solution = solver.solve(model, t_eval, inputs={"rate": 1.1})
        self.assertLess(len(solution.t), len(t_eval))
        np.testing.assert_allclose(
            solution.y.full()[0], np.exp(-1.1 * solution.t), rtol=1e-04
        )

    def test_model_solver_dae_inputs_in_initial_conditions(self):
        # Create model
        model = pybamm.BaseModel()
        var1 = pybamm.Variable("var1")
        var2 = pybamm.Variable("var2")
        model.rhs = {var1: pybamm.InputParameter("rate") * var1}
        model.algebraic = {var2: var1 - var2}
        model.initial_conditions = {
            var1: pybamm.InputParameter("ic 1"),
            var2: pybamm.InputParameter("ic 2"),
        }

        # Solve
        solver = pybamm.CasadiSolver(rtol=1e-8, atol=1e-8)
        t_eval = np.linspace(0, 5, 100)
        solution = solver.solve(
            model, t_eval, inputs={"rate": -1, "ic 1": 0.1, "ic 2": 2}
        )
        np.testing.assert_array_almost_equal(
            solution.y.full()[0], 0.1 * np.exp(-solution.t), decimal=5
        )
        np.testing.assert_array_almost_equal(
            solution.y.full()[-1], 0.1 * np.exp(-solution.t), decimal=5
        )

        # Solve again with different initial conditions
        solution = solver.solve(
            model, t_eval, inputs={"rate": -0.1, "ic 1": 1, "ic 2": 3}
        )
        np.testing.assert_array_almost_equal(
            solution.y.full()[0], 1 * np.exp(-0.1 * solution.t), decimal=5
        )
        np.testing.assert_array_almost_equal(
            solution.y.full()[-1], 1 * np.exp(-0.1 * solution.t), decimal=5
        )

    def test_model_solver_with_external(self):
        # Create model
        model = pybamm.BaseModel()
        domain = ["negative electrode", "separator", "positive electrode"]
        var1 = pybamm.Variable("var1", domain=domain)
        var2 = pybamm.Variable("var2", domain=domain)
        model.rhs = {var1: -var2}
        model.initial_conditions = {var1: 1}
        model.external_variables = [var2]
        model.variables = {"var1": var1, "var2": var2}
        # No need to set parameters; can use base discretisation (no spatial
        # operators)

        # create discretisation
        mesh = get_mesh_for_testing()
        spatial_methods = {"macroscale": pybamm.FiniteVolume()}
        disc = pybamm.Discretisation(mesh, spatial_methods)
        disc.process_model(model)
        # Solve
        solver = pybamm.CasadiSolver(rtol=1e-8, atol=1e-8)
        t_eval = np.linspace(0, 10, 100)
        solution = solver.solve(model, t_eval, external_variables={"var2": 0.5})
        np.testing.assert_allclose(
            solution.y.full()[0], 1 - 0.5 * solution.t, rtol=1e-06
        )

    def test_model_solver_with_non_identity_mass(self):
        model = pybamm.BaseModel()
        var1 = pybamm.Variable("var1", domain="negative electrode")
        var2 = pybamm.Variable("var2", domain="negative electrode")
        model.rhs = {var1: var1}
        model.algebraic = {var2: 2 * var1 - var2}
        model.initial_conditions = {var1: 1, var2: 2}
        disc = get_discretisation_for_testing()
        disc.process_model(model)

        # FV discretisation has identity mass. Manually set the mass matrix to
        # be a diag of 10s here for testing. Note that the algebraic part is all
        # zeros
        mass_matrix = 10 * model.mass_matrix.entries
        model.mass_matrix = pybamm.Matrix(mass_matrix)

        # Note that mass_matrix_inv is just the inverse of the ode block of the
        # mass matrix
        mass_matrix_inv = 0.1 * eye(int(mass_matrix.shape[0] / 2))
        model.mass_matrix_inv = pybamm.Matrix(mass_matrix_inv)

        # Solve
        solver = pybamm.CasadiSolver(rtol=1e-8, atol=1e-8)
        t_eval = np.linspace(0, 1, 100)
        solution = solver.solve(model, t_eval)
        np.testing.assert_array_equal(solution.t, t_eval)
        np.testing.assert_allclose(solution.y.full()[0], np.exp(0.1 * solution.t))
        np.testing.assert_allclose(solution.y.full()[-1], 2 * np.exp(0.1 * solution.t))

    def test_dae_solver_algebraic_model(self):
        model = pybamm.BaseModel()
        var = pybamm.Variable("var")
        model.algebraic = {var: var + 1}
        model.initial_conditions = {var: 0}

        disc = pybamm.Discretisation()
        disc.process_model(model)

        solver = pybamm.CasadiSolver()
        t_eval = np.linspace(0, 1)
        with self.assertRaisesRegex(
            pybamm.SolverError, "Cannot use CasadiSolver to solve algebraic model"
        ):
            solver.solve(model, t_eval)

    def test_interpolant_extrapolate(self):
        model = pybamm.lithium_ion.DFN()
        param = pybamm.ParameterValues("NCA_Kim2011")
        experiment = pybamm.Experiment(
            ["Charge at 1C until 4.6 V"], period="10 seconds"
        )

        param["Upper voltage cut-off [V]"] = 4.8

        sim = pybamm.Simulation(
            model,
            parameter_values=param,
            experiment=experiment,
            solver=pybamm.CasadiSolver(
                mode="safe",
                dt_max=0.001,
                extrap_tol=1e-3,
                extra_options_setup={"max_num_steps": 500},
            ),
        )
        with self.assertRaisesRegex(pybamm.SolverError, "interpolation bounds"):
            sim.solve()

        ci = param["Initial concentration in positive electrode [mol.m-3]"]
        param["Initial concentration in positive electrode [mol.m-3]"] = 0.8 * ci

        sim = pybamm.Simulation(
            model,
            parameter_values=param,
            experiment=experiment,
            solver=pybamm.CasadiSolver(mode="safe", dt_max=0.05),
        )
        with self.assertRaisesRegex(pybamm.SolverError, "interpolation bounds"):
            sim.solve()

    def test_casadi_safe_no_termination(self):
        model = pybamm.BaseModel()
        v = pybamm.Variable("v")
        model.rhs = {v: -1}
        model.initial_conditions = {v: 1}
        model.events.append(
            pybamm.Event(
                "Triggered event",
                v - 0.5,
                pybamm.EventType.INTERPOLANT_EXTRAPOLATION,
            )
        )
        model.events.append(
            pybamm.Event(
                "Ignored event",
                v + 10,
                pybamm.EventType.INTERPOLANT_EXTRAPOLATION,
            )
        )
        solver = pybamm.CasadiSolver(mode="safe")
        solver.set_up(model)

        with self.assertRaisesRegex(pybamm.SolverError, "interpolation bounds"):
            solver.solve(model, t_eval=[0, 1])


class TestCasadiSolverODEsWithForwardSensitivityEquations(unittest.TestCase):
    def test_solve_sensitivity_scalar_var_scalar_input(self):
        # Create model
        model = pybamm.BaseModel()
        var = pybamm.Variable("var")
        p = pybamm.InputParameter("p")
        model.rhs = {var: p * var}
        model.initial_conditions = {var: 1}
        model.variables = {"var squared": var ** 2}

        # Solve
        # Make sure that passing in extra options works
        solver = pybamm.CasadiSolver(mode="fast", rtol=1e-10, atol=1e-10)
        t_eval = np.linspace(0, 1, 80)
        solution = solver.solve(
            model, t_eval, inputs={"p": 0.1}, calculate_sensitivities=True
        )
        np.testing.assert_array_equal(solution.t, t_eval)
        np.testing.assert_allclose(solution.y[0], np.exp(0.1 * solution.t))
        np.testing.assert_allclose(
            solution.sensitivities["p"],
            (solution.t * np.exp(0.1 * solution.t))[:, np.newaxis],
        )
        np.testing.assert_allclose(
            solution["var squared"].data, np.exp(0.1 * solution.t) ** 2
        )
        np.testing.assert_allclose(
            solution["var squared"].sensitivities["p"],
            (2 * np.exp(0.1 * solution.t) * solution.t * np.exp(0.1 * solution.t))[
                :, np.newaxis
            ],
        )

        # More complicated model
        # Create model
        model = pybamm.BaseModel()
        var = pybamm.Variable("var")
        p = pybamm.InputParameter("p")
        q = pybamm.InputParameter("q")
        r = pybamm.InputParameter("r")
        s = pybamm.InputParameter("s")
        model.rhs = {var: p * q}
        model.initial_conditions = {var: r}
        model.variables = {"var times s": var * s}

        # Solve
        # Make sure that passing in extra options works
        solver = pybamm.CasadiSolver(rtol=1e-10, atol=1e-10)
        t_eval = np.linspace(0, 1, 80)
        solution = solver.solve(
            model,
            t_eval,
            inputs={"r": -1, "s": 0.5, "q": 2, "p": 0.1},
            calculate_sensitivities=True,
        )

        np.testing.assert_allclose(solution.y[0], -1 + 0.2 * solution.t)
        np.testing.assert_allclose(
            solution.sensitivities["p"],
            (2 * solution.t)[:, np.newaxis],
        )
        np.testing.assert_allclose(
            solution.sensitivities["q"],
            (0.1 * solution.t)[:, np.newaxis],
        )
        np.testing.assert_allclose(solution.sensitivities["r"], 1)
        np.testing.assert_allclose(solution.sensitivities["s"], 0)
        np.testing.assert_allclose(
            solution.sensitivities["all"],
            np.hstack(
                [
                    solution.sensitivities["p"],
                    solution.sensitivities["q"],
                    solution.sensitivities["r"],
                    solution.sensitivities["s"],
                ]
            ),
        )
        np.testing.assert_allclose(
            solution["var times s"].data, 0.5 * (-1 + 0.2 * solution.t)
        )
        np.testing.assert_allclose(
            solution["var times s"].sensitivities["p"],
            0.5 * (2 * solution.t)[:, np.newaxis],
        )
        np.testing.assert_allclose(
            solution["var times s"].sensitivities["q"],
            0.5 * (0.1 * solution.t)[:, np.newaxis],
        )
        np.testing.assert_allclose(solution["var times s"].sensitivities["r"], 0.5)
        np.testing.assert_allclose(
            solution["var times s"].sensitivities["s"],
            (-1 + 0.2 * solution.t)[:, np.newaxis],
        )
        np.testing.assert_allclose(
            solution["var times s"].sensitivities["all"],
            np.hstack(
                [
                    solution["var times s"].sensitivities["p"],
                    solution["var times s"].sensitivities["q"],
                    solution["var times s"].sensitivities["r"],
                    solution["var times s"].sensitivities["s"],
                ]
            ),
        )

    def test_solve_sensitivity_vector_var_scalar_input(self):
        var = pybamm.Variable("var", "negative electrode")
        model = pybamm.BaseModel()
        # Set length scales to avoid warning
        model.length_scales = {"negative electrode": 1}
        param = pybamm.InputParameter("param")
        model.rhs = {var: -param * var}
        model.initial_conditions = {var: 2}
        model.variables = {"var": var}

        # create discretisation
        disc = get_discretisation_for_testing()
        disc.process_model(model)
        n = disc.mesh["negative electrode"].npts

        # Solve - scalar input
        solver = pybamm.CasadiSolver()
        t_eval = np.linspace(0, 1)
        solution = solver.solve(
            model, t_eval, inputs={"param": 7}, calculate_sensitivities=["param"]
        )
        np.testing.assert_array_almost_equal(
            solution["var"].data,
            np.tile(2 * np.exp(-7 * t_eval), (n, 1)),
            decimal=4,
        )
        np.testing.assert_array_almost_equal(
            solution["var"].sensitivities["param"],
            np.repeat(-2 * t_eval * np.exp(-7 * t_eval), n)[:, np.newaxis],
            decimal=4,
        )

        # More complicated model
        # Create model
        model = pybamm.BaseModel()
        # Set length scales to avoid warning
        model.length_scales = {"negative electrode": 1}
        var = pybamm.Variable("var", "negative electrode")
        p = pybamm.InputParameter("p")
        q = pybamm.InputParameter("q")
        r = pybamm.InputParameter("r")
        s = pybamm.InputParameter("s")
        model.rhs = {var: p * q}
        model.initial_conditions = {var: r}
        model.variables = {"var times s": var * s}

        # Discretise
        disc.process_model(model)

        # Solve
        # Make sure that passing in extra options works
        solver = pybamm.CasadiSolver(
            rtol=1e-10,
            atol=1e-10,
        )
        t_eval = np.linspace(0, 1, 80)
        solution = solver.solve(
            model,
            t_eval,
            inputs={"p": 0.1, "q": 2, "r": -1, "s": 0.5},
            calculate_sensitivities=True,
        )
        np.testing.assert_allclose(solution.y, np.tile(-1 + 0.2 * solution.t, (n, 1)))
        np.testing.assert_allclose(
            solution.sensitivities["p"],
            np.repeat(2 * solution.t, n)[:, np.newaxis],
        )
        np.testing.assert_allclose(
            solution.sensitivities["q"],
            np.repeat(0.1 * solution.t, n)[:, np.newaxis],
        )
        np.testing.assert_allclose(solution.sensitivities["r"], 1)
        np.testing.assert_allclose(solution.sensitivities["s"], 0)
        np.testing.assert_allclose(
            solution.sensitivities["all"],
            np.hstack(
                [
                    solution.sensitivities["p"],
                    solution.sensitivities["q"],
                    solution.sensitivities["r"],
                    solution.sensitivities["s"],
                ]
            ),
        )
        np.testing.assert_allclose(
            solution["var times s"].data, np.tile(0.5 * (-1 + 0.2 * solution.t), (n, 1))
        )
        np.testing.assert_allclose(
            solution["var times s"].sensitivities["p"],
            np.repeat(0.5 * (2 * solution.t), n)[:, np.newaxis],
        )
        np.testing.assert_allclose(
            solution["var times s"].sensitivities["q"],
            np.repeat(0.5 * (0.1 * solution.t), n)[:, np.newaxis],
        )
        np.testing.assert_allclose(solution["var times s"].sensitivities["r"], 0.5)
        np.testing.assert_allclose(
            solution["var times s"].sensitivities["s"],
            np.repeat(-1 + 0.2 * solution.t, n)[:, np.newaxis],
        )
        np.testing.assert_allclose(
            solution["var times s"].sensitivities["all"],
            np.hstack(
                [
                    solution["var times s"].sensitivities["p"],
                    solution["var times s"].sensitivities["q"],
                    solution["var times s"].sensitivities["r"],
                    solution["var times s"].sensitivities["s"],
                ]
            ),
        )

    def test_solve_sensitivity_scalar_var_vector_input(self):
        var = pybamm.Variable("var", "negative electrode")
        model = pybamm.BaseModel()
        # Set length scales to avoid warning
        model.length_scales = {"negative electrode": 1}

        param = pybamm.InputParameter("param", "negative electrode")
        model.rhs = {var: -param * var}
        model.initial_conditions = {var: 2}
        model.variables = {
            "var": var,
            "integral of var": pybamm.Integral(var, pybamm.standard_spatial_vars.x_n),
        }

        # create discretisation
        mesh = get_mesh_for_testing(xpts=5)
        spatial_methods = {"macroscale": pybamm.FiniteVolume()}
        disc = pybamm.Discretisation(mesh, spatial_methods)
        disc.process_model(model)
        n = disc.mesh["negative electrode"].npts

        # Solve - constant input
        solver = pybamm.CasadiSolver(mode="fast", rtol=1e-10, atol=1e-10)
        t_eval = np.linspace(0, 1)
        solution = solver.solve(
            model,
            t_eval,
            inputs={"param": 7 * np.ones(n)},
            calculate_sensitivities=True,
        )
        l_n = mesh["negative electrode"].edges[-1]
        np.testing.assert_array_almost_equal(
            solution["var"].data,
            np.tile(2 * np.exp(-7 * t_eval), (n, 1)),
            decimal=4,
        )

        np.testing.assert_array_almost_equal(
            solution["var"].sensitivities["param"],
            np.vstack([np.eye(n) * -2 * t * np.exp(-7 * t) for t in t_eval]),
        )
        np.testing.assert_array_almost_equal(
            solution["integral of var"].data,
            2 * np.exp(-7 * t_eval) * l_n,
            decimal=4,
        )
        np.testing.assert_array_almost_equal(
            solution["integral of var"].sensitivities["param"],
            np.tile(-2 * t_eval * np.exp(-7 * t_eval) * l_n / n, (n, 1)).T,
        )

        # Solve - linspace input
        p_eval = np.linspace(1, 2, n)
        solution = solver.solve(
            model, t_eval, inputs={"param": p_eval}, calculate_sensitivities=True
        )
        l_n = mesh["negative electrode"].edges[-1]
        np.testing.assert_array_almost_equal(
            solution["var"].data, 2 * np.exp(-p_eval[:, np.newaxis] * t_eval), decimal=4
        )
        np.testing.assert_array_almost_equal(
            solution["var"].sensitivities["param"],
            np.vstack([np.diag(-2 * t * np.exp(-p_eval * t)) for t in t_eval]),
        )

        np.testing.assert_array_almost_equal(
            solution["integral of var"].data,
            np.sum(
                2
                * np.exp(-p_eval[:, np.newaxis] * t_eval)
                * mesh["negative electrode"].d_edges[:, np.newaxis],
                axis=0,
            ),
        )
        np.testing.assert_array_almost_equal(
            solution["integral of var"].sensitivities["param"],
            np.vstack([-2 * t * np.exp(-p_eval * t) * l_n / n for t in t_eval]),
        )

    def test_solve_sensitivity_then_no_sensitivity(self):
        # Create model
        model = pybamm.BaseModel()
        var = pybamm.Variable("var")
        p = pybamm.InputParameter("p")
        model.rhs = {var: p * var}
        model.initial_conditions = {var: 1}
        model.variables = {"var squared": var ** 2}

        # Solve
        # Make sure that passing in extra options works
        solver = pybamm.CasadiSolver(mode="fast", rtol=1e-10, atol=1e-10)
        t_eval = np.linspace(0, 1, 80)
        solution = solver.solve(
            model, t_eval, inputs={"p": 0.1}, calculate_sensitivities=True
        )

        # check sensitivities
        np.testing.assert_allclose(
            solution.sensitivities["p"],
            (solution.t * np.exp(0.1 * solution.t))[:, np.newaxis],
        )

        solution = solver.solve(model, t_eval, inputs={"p": 0.1})

        np.testing.assert_array_equal(solution.t, t_eval)
        np.testing.assert_allclose(solution.y, np.exp(0.1 * solution.t).reshape(1, -1))
        np.testing.assert_allclose(
            solution["var squared"].data, np.exp(0.1 * solution.t) ** 2
        )

    def test_solve_sensitivity_subset(self):
        # Create model
        model = pybamm.BaseModel()
        var = pybamm.Variable("var")
        p = pybamm.InputParameter("p")
        q = pybamm.InputParameter("q")
        r = pybamm.InputParameter("r")
        model.rhs = {var: p * q}
        model.initial_conditions = {var: r}

        # only calculate the sensitivities of a subset of parameters
        solver = pybamm.CasadiSolver(rtol=1e-10, atol=1e-10)
        t_eval = np.linspace(0, 1, 80)
        solution = solver.solve(
            model,
            t_eval,
            inputs={"q": 2, "r": -1, "p": 0.1},
            calculate_sensitivities=["q", "p"],
        )
        np.testing.assert_allclose(solution.y[0], -1 + 0.2 * solution.t)
        np.testing.assert_allclose(
            solution.sensitivities["p"],
            (2 * solution.t)[:, np.newaxis],
        )
        np.testing.assert_allclose(
            solution.sensitivities["q"],
            (0.1 * solution.t)[:, np.newaxis],
        )
        self.assertTrue("r" not in solution.sensitivities)
        np.testing.assert_allclose(
            solution.sensitivities["all"],
            np.hstack(
                [
                    solution.sensitivities["p"],
                    solution.sensitivities["q"],
                ]
            ),
        )

        solution = solver.solve(
            model,
            t_eval,
            inputs={"q": 2, "r": -1, "p": 0.1},
            calculate_sensitivities=["r"],
        )
        np.testing.assert_allclose(solution.y[0], -1 + 0.2 * solution.t)
        self.assertTrue("p" not in solution.sensitivities)
        self.assertTrue("q" not in solution.sensitivities)
        np.testing.assert_allclose(solution.sensitivities["r"], 1)
        np.testing.assert_allclose(
            solution.sensitivities["all"],
            np.hstack(
                [
                    solution.sensitivities["r"],
                ]
            ),
        )


class TestCasadiSolverDAEsWithForwardSensitivityEquations(unittest.TestCase):
    def test_solve_sensitivity_scalar_var_scalar_input(self):
        # Create model
        model = pybamm.BaseModel()
        var1 = pybamm.Variable("var1")
        p = pybamm.InputParameter("p")
        var1 = pybamm.Variable("var1")
        var2 = pybamm.Variable("var2")
        model.rhs = {var1: p * var1}
        model.algebraic = {var2: 2 * var1 - var2}
        model.initial_conditions = {var1: 1, var2: 2}
        model.variables = {"var2 squared": var2 ** 2}

        # Solve
        # Make sure that passing in extra options works
        solver = pybamm.CasadiSolver(mode="fast", rtol=1e-10, atol=1e-10)
        t_eval = np.linspace(0, 1, 80)
        solution = solver.solve(
            model, t_eval, inputs={"p": 0.1}, calculate_sensitivities=True
        )
        np.testing.assert_array_equal(solution.t, t_eval)
        np.testing.assert_allclose(solution.y[0], np.exp(0.1 * solution.t))
        np.testing.assert_allclose(
            solution.sensitivities["p"],
            np.stack(
                (
                    solution.t * np.exp(0.1 * solution.t),
                    2 * solution.t * np.exp(0.1 * solution.t),
                )
            )
            .transpose()
            .reshape(-1, 1),
            atol=1e-7,
        )
        np.testing.assert_allclose(
            solution["var2 squared"].data, 4 * np.exp(2 * 0.1 * solution.t)
        )
        np.testing.assert_allclose(
            solution["var2 squared"].sensitivities["p"],
            (8 * solution.t * np.exp(2 * 0.1 * solution.t))[:, np.newaxis],
            atol=1e-7,
        )

    def test_solve_sensitivity_algebraic(self):
        # Create model
        model = pybamm.BaseModel()
        var = pybamm.Variable("var")
        p = pybamm.InputParameter("p")
        model.algebraic = {var: var - p * pybamm.t}
        model.initial_conditions = {var: 0}
        model.variables = {"var squared": var ** 2}

        # Solve
        # Make sure that passing in extra options works
        solver = pybamm.CasadiAlgebraicSolver(tol=1e-10)
        t_eval = np.linspace(0, 1, 80)
        solution = solver.solve(
            model, t_eval, inputs={"p": 0.1}, calculate_sensitivities=True
        )
        np.testing.assert_array_equal(solution.t, t_eval)
        np.testing.assert_allclose(solution.y[0], 0.1 * solution.t)
        np.testing.assert_allclose(
            solution.sensitivities["p"], solution.t.reshape(-1, 1), atol=1e-7
        )
        np.testing.assert_allclose(
            solution["var squared"].data, (0.1 * solution.t) ** 2
        )
        np.testing.assert_allclose(
            solution["var squared"].sensitivities["p"],
            (2 * 0.1 * solution.t ** 2).reshape(-1, 1),
            atol=1e-7,
        )

    def test_solve_sensitivity_subset(self):
        # Create model
        model = pybamm.BaseModel()
        var = pybamm.Variable("var")
        var2 = pybamm.Variable("var2")
        p = pybamm.InputParameter("p")
        q = pybamm.InputParameter("q")
        r = pybamm.InputParameter("r")
        model.rhs = {var: p * q}
        model.algebraic = {var2: 2 * var - var2}
        model.initial_conditions = {var: r, var2: 2 * r}

        # only calculate the sensitivities of a subset of parameters
        solver = pybamm.CasadiSolver(rtol=1e-10, atol=1e-10)
        t_eval = np.linspace(0, 1, 80)
        solution = solver.solve(
            model,
            t_eval,
            inputs={"p": 0.1, "q": 2, "r": -1, "s": 0.5},
            calculate_sensitivities=["p", "q"],
        )
        np.testing.assert_allclose(solution.y[0], -1 + 0.2 * solution.t)
        np.testing.assert_allclose(solution.y[-1], 2 * (-1 + 0.2 * solution.t))
        np.testing.assert_allclose(
            solution.sensitivities["p"][::2],
            (2 * solution.t)[:, np.newaxis],
        )
        np.testing.assert_allclose(
            solution.sensitivities["q"][::2],
            (0.1 * solution.t)[:, np.newaxis],
        )
        self.assertTrue("r" not in solution.sensitivities)
        self.assertTrue("s" not in solution.sensitivities)
        np.testing.assert_allclose(
            solution.sensitivities["all"],
            np.hstack(
                [
                    solution.sensitivities["p"],
                    solution.sensitivities["q"],
                ]
            ),
        )


if __name__ == "__main__":
    print("Add -v for more debug output")
    import sys

    if "-v" in sys.argv:
        debug = True
    pybamm.settings.debug_mode = True
    unittest.main()
