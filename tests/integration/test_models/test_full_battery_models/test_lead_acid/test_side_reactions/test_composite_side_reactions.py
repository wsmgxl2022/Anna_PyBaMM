#
# Tests for the lead-acid Full model
#
import pybamm
import tests

import unittest


class TestLeadAcidCompositeSideReactions(unittest.TestCase):
    def test_basic_processing_differential(self):
        options = {"hydrolysis": "true", "surface form": "differential"}
        model = pybamm.lead_acid.Composite(options)
        modeltest = tests.StandardModelTest(model)
        modeltest.test_all(skip_output_tests=True)

    def test_basic_processing_algebraic(self):
        options = {"hydrolysis": "true", "surface form": "algebraic"}
        model = pybamm.lead_acid.Composite(options)
        modeltest = tests.StandardModelTest(model)
        modeltest.test_all(skip_output_tests=True)

    def test_basic_processing_charge(self):
        options = {"hydrolysis": "true", "surface form": "differential"}
        model = pybamm.lead_acid.Composite(options)
        parameter_values = model.default_parameter_values
        parameter_values.update(
            {"Current function [A]": -1, "Initial State of Charge": 0.5}
        )
        modeltest = tests.StandardModelTest(model, parameter_values=parameter_values)
        modeltest.test_all(skip_output_tests=True)

    def test_basic_processing_zero_current(self):
        options = {"hydrolysis": "true", "surface form": "differential"}
        model = pybamm.lead_acid.Composite(options)
        parameter_values = model.default_parameter_values
        parameter_values.update({"Current function [A]": 0})
        modeltest = tests.StandardModelTest(model, parameter_values=parameter_values)
        modeltest.test_all(skip_output_tests=True)

    def test_basic_processing_extended_differential(self):
        options = {"hydrolysis": "true", "surface form": "differential"}
        model = pybamm.lead_acid.CompositeExtended(options)
        modeltest = tests.StandardModelTest(model)
        modeltest.test_all(skip_output_tests=True)


if __name__ == "__main__":
    print("Add -v for more debug output")
    import sys

    if "-v" in sys.argv:
        debug = True
    unittest.main()
