#
# Tests for the half-cell lithium-ion SPMe model
#
import pybamm
import unittest
from tests import BaseIntegrationTestLithiumIonHalfCell


class TestSPMeHalfCell(BaseIntegrationTestLithiumIonHalfCell, unittest.TestCase):
    def setUp(self):
        self.model = pybamm.lithium_ion.SPMe


if __name__ == "__main__":
    print("Add -v for more debug output")
    import sys

    if "-v" in sys.argv:
        debug = True
    unittest.main()
