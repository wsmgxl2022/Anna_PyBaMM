#
# Tests for the Dummy Solver class
#
import pybamm
import numpy as np
import unittest
import sys


class TestDummySolver(unittest.TestCase):
    def test_dummy_solver(self):
        model = pybamm.BaseModel()
        v = pybamm.Scalar(1)
        model.variables = {"v": v}

        disc = pybamm.Discretisation()
        disc.process_model(model)

        solver = pybamm.DummySolver()
        t_eval = np.linspace(0, 1)
        sol = solver.solve(model, t_eval)
        np.testing.assert_array_equal(sol.t, t_eval)
        np.testing.assert_array_equal(sol.y, np.zeros((1, t_eval.size)))
        np.testing.assert_array_equal(sol["v"].data, np.ones(t_eval.size))

    def test_dummy_solver_step(self):
        model = pybamm.BaseModel()
        v = pybamm.Scalar(1)
        model.variables = {"v": v}

        disc = pybamm.Discretisation()
        disc.process_model(model)

        solver = pybamm.DummySolver()
        t_eval = np.linspace(0, 1)

        sol = None
        for dt in np.diff(t_eval):
            sol = solver.step(sol, model, dt)

        np.testing.assert_array_equal(sol.t, t_eval)
        np.testing.assert_array_equal(sol.y, np.zeros((1, t_eval.size)))
        np.testing.assert_array_equal(sol["v"].data, np.ones(t_eval.size))


if __name__ == "__main__":
    print("Add -v for more debug output")

    if "-v" in sys.argv:
        debug = True
    pybamm.settings.debug_mode = True
    unittest.main()
