# This contains tests for the functions


import unittest
import numpy as np
import sys
import os

# Add the parent directory of the examples folder to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from boolnetanalyzer.utils import generate_indegrees
from boolnetanalyzer.bn_tools import is_degenerated,num_attractors
from boolnetanalyzer import randomize as rn
class TestIsDegenerated(unittest.TestCase):
    def test_is_degenerated(self):
        self.assertFalse(is_degenerated(np.array([0, 1, 1, 0])))
        self.assertFalse(is_degenerated(np.array([0, 1, 1, 1])))
        self.assertTrue(is_degenerated(np.array([1, 1, 1, 1])))

class TestIndegreeDistribution(unittest.TestCase):
    def setUp(self):
        # Set a seed for reproducibility if needed for each test
        np.random.seed(0)

    def test_indegree_distribution_poisson(self):
        result = generate_indegrees(10, 'poisson', 2)
        expected = np.array([3, 2, 5, 1, 1, 1, 7, 1, 3, 3])
        np.testing.assert_array_equal(result, expected)
    
    def test_indegree_distribution_uniform(self):
        result = generate_indegrees(10, 'uniform', 3)
        expected = np.array([1, 2, 1, 2, 2, 3, 1, 3, 1, 1])
        np.testing.assert_array_equal(result, expected)
    
    def test_indegree_distribution_constant(self):
        result = generate_indegrees(10, 'constant', 2)
        expected = np.array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2])
        np.testing.assert_array_equal(result, expected)
    
    def test_indegree_distribution_vector(self):
        result = generate_indegrees(10, indegree_vector=[1, 2, 1, 2, 1, 1, 3, 3, 3, 1])
        expected = np.array([1, 2, 1, 2, 1, 1, 3, 3, 3, 1])
        np.testing.assert_array_equal(result, expected)


class TestNumAttractors(unittest.TestCase):
    def test_num_attractors(self):
        # Define the expected results
        expected_attractors = [[11, 22], [27, 30, 15, 21]]
        expected_basin_sizes = [7, 20]
        expected_attr_dict = {0: 0, 11: 0, 22: 0, 1: 1, 19: 1, 27: 1, 30: 1, 15: 1, 21: 1, 2: 1, 3: 1, 4: 1, 5: 1, 17: 1, 6: 1, 7: 1, 8: 1, 14: 1, 9: 0, 18: 0, 10: 1, 12: 1, 13: 1, 16: 0, 20: 0, 23: 1, 24: 1}
        expected_update_dict = {0: 11, 11: 22, 22: 11, 1: 19, 19: 27, 27: 30, 30: 15, 15: 21, 21: 27, 2: 3, 3: 19, 4: 1, 5: 17, 17: 19, 6: 1, 7: 17, 8: 14, 14: 5, 9: 18, 18: 11, 10: 6, 12: 5, 13: 17, 16: 11, 20: 11, 23: 19, 24: 14}

        # Run the function with the test input
        attractors, basin_sizes, attr_dict, update_dict = num_attractors(
            rn.boolnet(5, seed=9), initial_sample_points=[i for i in range(25)]
        )

        # Assert the results match the expected values
        self.assertEqual(attractors, expected_attractors)
        self.assertEqual(basin_sizes, expected_basin_sizes)
        self.assertEqual(attr_dict, expected_attr_dict)
        self.assertEqual(update_dict, expected_update_dict)
if __name__ == '__main__':
    unittest.main()

