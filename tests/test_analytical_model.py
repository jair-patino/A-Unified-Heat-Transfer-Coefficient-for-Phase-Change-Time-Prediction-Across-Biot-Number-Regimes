import unittest
from src.analytical_model import global_U, phase_change_time

class TestAnalyticalModel(unittest.TestCase):
    def test_global_U(self):
        # Test for small Bi
        heff = 10
        k = 0.6
        Lc = 0.01
        Phi = 1.0
        U = global_U(heff, k, Lc, Phi)
        Bi = heff * Lc / k
        expected = heff / (1 + Phi * Bi)
        self.assertAlmostEqual(U, expected)

    # Other tests...

if __name__ == '__main__':
    unittest.main()
