#
# Tests particle size distribution parameters are loaded into a parameter set
# and give expected values
#
import pybamm

import unittest
import numpy as np


class TestSizeDistributionParameters(unittest.TestCase):
    def test_parameter_values(self):
        values = pybamm.lithium_ion.BaseModel().default_parameter_values
        param = pybamm.LithiumIonParameters()

        # add distribution parameter values
        values = pybamm.get_size_distribution_parameters(values)

        # check dimensionless parameters

        # min and max radii
        np.testing.assert_almost_equal(values.evaluate(param.n.prim.R_min), 0.0, 3)
        np.testing.assert_almost_equal(values.evaluate(param.p.prim.R_min), 0.0, 3)
        np.testing.assert_almost_equal(values.evaluate(param.n.prim.R_max), 2.5, 3)
        np.testing.assert_almost_equal(values.evaluate(param.p.prim.R_max), 2.5, 3)

        # standard deviations
        np.testing.assert_almost_equal(values.evaluate(param.n.prim.sd_a), 0.3, 3)
        np.testing.assert_almost_equal(values.evaluate(param.p.prim.sd_a), 0.3, 3)

        # check function parameters (size distributions) evaluate
        R_test = pybamm.Scalar(1.0)
        values.evaluate(param.n.prim.f_a_dist(R_test))
        values.evaluate(param.p.prim.f_a_dist(R_test))


if __name__ == "__main__":
    print("Add -v for more debug output")
    import sys

    if "-v" in sys.argv:
        debug = True
    pybamm.settings.debug_mode = True
    unittest.main()
