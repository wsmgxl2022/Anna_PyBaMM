#
# Tests for the Solution class
#
import pybamm
import unittest
import numpy as np


class TestSolution(unittest.TestCase):
    def test_append(self):
        model = pybamm.lithium_ion.SPMe()
        # create geometry
        geometry = model.default_geometry

        # load parameter values and process model and geometry
        param = model.default_parameter_values
        param.process_model(model)
        param.process_geometry(geometry)

        # set mesh
        mesh = pybamm.Mesh(geometry, model.default_submesh_types, model.default_var_pts)

        # discretise model
        disc = pybamm.Discretisation(mesh, model.default_spatial_methods)
        disc.process_model(model)

        # solve model
        t_eval = np.linspace(0, 3600, 100)
        solver = model.default_solver
        solution = solver.solve(model, t_eval)

        # step model
        old_t = 0
        step_solver = model.default_solver
        step_solution = None
        # dt should be dimensional
        solution_times_dimensional = solution.t * model.timescale_eval
        for t in solution_times_dimensional[1:]:
            dt = t - old_t
            step_solution = step_solver.step(step_solution, model, dt=dt, npts=10)
            if t == solution_times_dimensional[1]:
                # Create voltage variable
                step_solution.update("Terminal voltage")
            old_t = t

        # Check both give the same answer
        np.testing.assert_array_almost_equal(
            solution["Terminal voltage"](solution.t[:-1] * model.timescale_eval),
            step_solution["Terminal voltage"](solution.t[:-1] * model.timescale_eval),
            decimal=4,
        )

    def test_append_external_variables(self):
        model = pybamm.lithium_ion.SPM(
            {
                "thermal": "lumped",
                "external submodels": ["thermal", "negative primary particle"],
            }
        )
        # create geometry
        geometry = model.default_geometry

        # load parameter values and process model and geometry
        param = model.default_parameter_values
        param.process_model(model)
        param.process_geometry(geometry)

        # set mesh
        mesh = pybamm.Mesh(geometry, model.default_submesh_types, model.default_var_pts)

        # discretise model
        disc = pybamm.Discretisation(mesh, model.default_spatial_methods)
        disc.process_model(model)

        # solve model
        solver = model.default_solver

        Nr = model.default_var_pts["r_n"]

        T_av = 0
        c_s_n_av = np.ones((Nr, 1)) * 0.6
        external_variables = {
            "Volume-averaged cell temperature": T_av,
            "X-averaged negative particle concentration": c_s_n_av,
        }

        # Step
        dt = 0.1
        sol_step = None
        for _ in range(5):
            sol_step = solver.step(
                sol_step, model, dt, external_variables=external_variables
            )
        np.testing.assert_array_equal(
            sol_step.all_inputs[0]["Volume-averaged cell temperature"], 0
        )
        np.testing.assert_array_equal(
            sol_step.all_inputs[0]["X-averaged negative particle concentration"], 0.6
        )

        # Solve
        t_eval = np.linspace(0, 3600)
        sol = solver.solve(model, t_eval, external_variables=external_variables)
        np.testing.assert_array_equal(
            sol.all_inputs[0]["Volume-averaged cell temperature"], 0
        )
        np.testing.assert_array_equal(
            sol.all_inputs[0]["X-averaged negative particle concentration"], 0.6
        )


if __name__ == "__main__":
    print("Add -v for more debug output")
    import sys

    if "-v" in sys.argv:
        debug = True
    pybamm.settings.debug_mode = True
    unittest.main()
