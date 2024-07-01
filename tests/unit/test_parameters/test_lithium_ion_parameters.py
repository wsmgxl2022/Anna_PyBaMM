#
# Tests lithium ion parameters load and give expected values
#
import pybamm

import unittest
import numpy as np


class TestDimensionlessParameterValues(unittest.TestCase):
    def test_print_parameters(self):
        parameters = pybamm.LithiumIonParameters()
        parameter_values = pybamm.lithium_ion.BaseModel().default_parameter_values
        output_file = "lithium_ion_parameters.txt"
        parameter_values.print_parameters(parameters, output_file)
        # test print_parameters with dict and without C-rate
        del parameter_values["Nominal cell capacity [A.h]"]
        parameters = {"C_e": parameters.C_e, "sigma_n": parameters.n.sigma}
        parameter_values.print_parameters(parameters)

    def test_lithium_ion(self):
        """This test checks that all the dimensionless parameters are being calculated
        correctly for the specific set of parameters for LCO from dualfoil. The values
        are those converted from those in Scott's transfer which previous versions of
        the DFN work with. A 1C rate corresponds to a 24A/m^2 current density"""
        values = pybamm.lithium_ion.BaseModel().default_parameter_values
        param = pybamm.LithiumIonParameters()

        c_rate = param.i_typ / 24  # roughly for the numbers I used before

        # particle geometry
        # Note: in general these can be functions, but are constant for this
        # set, so we just arbitrarily evaluate at 0

        # a_n_typ
        np.testing.assert_almost_equal(
            values.evaluate(param.n.prim.a_typ), 0.18 * 10 ** (6), 2
        )
        # R_n dimensional
        np.testing.assert_almost_equal(
            values.evaluate(param.n.prim.R_typ), 1 * 10 ** (-5), 2
        )

        # a_R_n = a_n_typ * R_n_typ
        np.testing.assert_almost_equal(values.evaluate(param.n.prim.a_R), 1.8, 2)

        # a_p_typ
        np.testing.assert_almost_equal(
            values.evaluate(param.p.prim.a_typ), 0.15 * 10 ** (6), 2
        )

        # R_p dimensional
        np.testing.assert_almost_equal(
            values.evaluate(param.p.prim.R_typ), 1 * 10 ** (-5), 2
        )

        # a_p = a_p_typ * R_p_typ
        np.testing.assert_almost_equal(values.evaluate(param.p.prim.a_R), 1.5, 2)

        # j0_m
        np.testing.assert_almost_equal(
            values.evaluate(
                param.n.prim.j0_dimensional(
                    param.c_e_typ, param.n.prim.c_max / 2, param.T_ref
                )
            ),
            values.evaluate(
                2 * 10 ** (-5) * param.c_e_typ ** 0.5 * param.n.prim.c_max / 2
            ),
            8,
        )

        # j0_p
        np.testing.assert_almost_equal(
            values.evaluate(
                param.p.prim.j0_dimensional(
                    param.c_e_typ, param.p.prim.c_max / 2, param.T_ref
                )
            ),
            values.evaluate(
                6 * 10 ** (-7) * param.c_e_typ ** 0.5 * param.p.prim.c_max / 2
            ),
            8,
        )

        # particle dynamics
        # neg diffusion coefficient
        np.testing.assert_almost_equal(
            values.evaluate(
                pybamm.xyz_average(
                    pybamm.r_average(
                        param.n.prim.D_dimensional(param.n.prim.c_init, param.T_ref)
                    )
                )
            ),
            3.9 * 10 ** (-14),
            2,
        )

        # neg diffusion timescale
        np.testing.assert_almost_equal(
            values.evaluate(param.n.prim.tau_diffusion), 2.5641 * 10 ** (3), 2
        )

        # tau_n / tau_d (1/gamma_n in Scott's transfer)
        np.testing.assert_almost_equal(
            values.evaluate(param.n.prim.C_diff / c_rate), 0.11346, 3
        )

        # pos diffusion coefficient
        np.testing.assert_almost_equal(
            values.evaluate(
                pybamm.xyz_average(
                    pybamm.r_average(
                        param.p.prim.D_dimensional(param.p.prim.c_init, param.T_ref)
                    )
                )
            ),
            1 * 10 ** (-13),
            2,
        )

        # pos diffusion timescale
        np.testing.assert_almost_equal(
            values.evaluate(param.p.prim.tau_diffusion), 1 * 10 ** (3), 2
        )

        # tau_p / tau_d (1/gamma_p in Scott's transfer)
        np.testing.assert_almost_equal(
            values.evaluate(param.p.prim.C_diff / c_rate), 0.044249, 3
        )

        # electrolyte dynamics
        # typical diffusion coefficient (we should change the typ value in paper to
        # match this one. We take this parameter excluding the exp(-0.65) in the
        # paper at the moment
        np.testing.assert_almost_equal(
            values.evaluate(param.D_e_dimensional(param.c_e_typ, param.T_ref)),
            5.34 * 10 ** (-10) * np.exp(-0.65),
            10,
        )

        # electrolyte diffusion timescale (accounting for np.exp(-0.65) in
        # diffusion_typ). Change value in paper to this.
        np.testing.assert_almost_equal(
            values.evaluate(param.tau_diffusion_e), 181.599, 3
        )

        # C_e
        np.testing.assert_almost_equal(values.evaluate(param.C_e / c_rate), 0.008, 3)

        # electrolyte conductivity
        np.testing.assert_almost_equal(
            values.evaluate(param.kappa_e_dimensional(param.c_e_typ, param.T_ref)),
            1.1045,
            3,
        )

        # potential scale
        # F R / T (should be equal to old 1 / Lambda)
        old_Lambda = 38
        np.testing.assert_almost_equal(
            values.evaluate(param.potential_scale), 1 / old_Lambda, 3
        )

        # electrode conductivities
        # neg dimensional
        np.testing.assert_almost_equal(
            values.evaluate(param.n.sigma_dimensional(param.T_ref)), 100, 3
        )

        # neg dimensionless (old sigma_n / old_Lambda ) (this is different to values
        # in paper so check again, it is close enough though for now)
        np.testing.assert_almost_equal(
            values.evaluate(param.n.sigma(param.T_ref) * c_rate), 475.7, 1
        )

        # neg dimensionless rescaled
        np.testing.assert_almost_equal(
            values.evaluate(param.n.sigma_prime(param.T_ref) * c_rate), 0.7814, 1
        )

        # pos dimensional
        np.testing.assert_almost_equal(
            values.evaluate(param.p.sigma_dimensional(param.T_ref)), 10, 3
        )

        # pos dimensionless (old sigma_n / old_Lambda ) (this is different to values in
        # paper so check again, it is close enough for now though)
        np.testing.assert_almost_equal(
            values.evaluate(param.p.sigma(param.T_ref) * c_rate), 47.57, 1
        )

        # pos dimensionless rescaled
        np.testing.assert_almost_equal(
            values.evaluate(param.p.sigma_prime(param.T_ref) * c_rate), 0.07814, 1
        )

    def test_thermal_parameters(self):
        values = pybamm.lithium_ion.BaseModel().default_parameter_values
        param = pybamm.LithiumIonParameters()
        c_rate = param.i_typ / 24
        T = 1  # dummy temperature as the values are constant

        # Density
        np.testing.assert_almost_equal(values.evaluate(param.n.rho_cc(T)), 1.9019, 2)
        np.testing.assert_almost_equal(values.evaluate(param.n.rho(T)), 0.6403, 2)
        np.testing.assert_almost_equal(values.evaluate(param.s.rho(T)), 0.1535, 2)
        np.testing.assert_almost_equal(values.evaluate(param.p.rho(T)), 1.2605, 2)
        np.testing.assert_almost_equal(values.evaluate(param.p.rho_cc(T)), 1.3403, 2)

        # Thermal conductivity
        np.testing.assert_almost_equal(values.evaluate(param.n.lambda_cc(T)), 6.7513, 2)
        np.testing.assert_almost_equal(values.evaluate(param.n.lambda_(T)), 0.0296, 2)
        np.testing.assert_almost_equal(values.evaluate(param.s.lambda_(T)), 0.0027, 2)
        np.testing.assert_almost_equal(values.evaluate(param.p.lambda_(T)), 0.0354, 2)
        np.testing.assert_almost_equal(values.evaluate(param.p.lambda_cc(T)), 3.9901, 2)

        # other thermal parameters

        # note: in paper this is 0.0534 * c_rate which conflicts with this
        # if we do C_th * c_rate we get 0.0534 so probably error in paper
        # np.testing.assert_almost_equal(
        #     values.evaluate(param.C_th / c_rate), 0.0253, 2
        # )

        np.testing.assert_almost_equal(values.evaluate(param.Theta / c_rate), 0.008, 2)

        # np.testing.assert_almost_equal(
        #     values.evaluate(param.B / c_rate), 36.216, 2
        # )

        np.testing.assert_equal(values.evaluate(param.T_init), 0)

        # test timescale
        # np.testing.assert_almost_equal(
        #     values.evaluate(param.tau_th_yz), 1.4762 * 10 ** (3), 2
        # )

        # thermal = pybamm.thermal_parameters
        # np.testing.assert_almost_equal(
        # values.evaluate(thermal.rho_eff_dim), 1.8116 * 10 ** (6), 2
        # )
        # np.testing.assert_almost_equal(
        #     values.evaluate(thermal.lambda_eff_dim), 59.3964, 2
        # )

    def test_parameter_functions(self):
        values = pybamm.lithium_ion.BaseModel().default_parameter_values
        param = pybamm.LithiumIonParameters()

        T_test = pybamm.Scalar(0)

        c_e_test = pybamm.Scalar(1)
        values.evaluate(param.D_e(c_e_test, T_test))
        values.evaluate(param.kappa_e(c_e_test, T_test))

    def test_timescale(self):
        param = pybamm.LithiumIonParameters({"timescale": 2.5})
        self.assertEqual(param.timescale.evaluate(), 2.5)


if __name__ == "__main__":
    print("Add -v for more debug output")
    import sys

    if "-v" in sys.argv:
        debug = True
    pybamm.settings.debug_mode = True
    unittest.main()
