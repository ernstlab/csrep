import unittest
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../scripts/')))
import helper

class TestHelperMethods(unittest.TestCase):
    def test_partition_file_list(self):
        data = [1,2,3,4]
        num_cores = 4
        exp_result = [[1], [2], [3], [4]]
        obs_result = helper.partition_file_list(data, num_cores)
        self.assertCountEqual(exp_result, obs_result)
        data = [1,2,3]
        exp_result = [[1], [2], [3], []]
        obs_result = helper.partition_file_list(data, num_cores)
        self.assertCountEqual(exp_result, obs_result)
        data = [1,2]
        exp_result = [[1], [2], [], []]
        obs_result = helper.partition_file_list(data, num_cores)
        self.assertCountEqual(exp_result, obs_result)
        data = [1,2,3,4,5,6]
        exp_result = [[1,2], [3,4], [5,6], []]
        obs_result = helper.partition_file_list(data, num_cores)
        self.assertCountEqual(exp_result, obs_result)
        return

if __name__ == "__main__":
    unittest.main()