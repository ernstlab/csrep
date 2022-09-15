import unittest
import os
import sys
import shutil
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../scripts/')))
import average_pred_results as avg
import helper


def create_fake_all_ct_segment_folder(all_ct_segment_folder, total_num_files):
    helper.make_dir(all_ct_segment_folder)
    for i in range(total_num_files):
        fn = os.path.join(all_ct_segment_folder, 'chr22_{}_combined_segment.bed.gz'.format(i))
        f = open(fn, 'w+')  # create this file if it does not exist
        f.close()
    return

def create_fake_avg_result_folder(avg_folder, num_existing_files):
    helper.make_dir(avg_folder)
    for i in range(num_existing_files):
        fn = os.path.join(avg_folder, 'chr22_{}_avg_pred.txt.gz'.format(i))
        f = open(fn, 'w+')
        f.close()
    return

class TestAvgMethods(unittest.TestCase):
    total_num_files = 6
    def test_get_genomic_positions_list(self):
        testdata_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), '../testdata/'))
        all_ct_segment_folder = os.path.join(testdata_folder, 'all_ct_segments')
        avg_folder = os.path.join(testdata_folder, 'avg_folder')
        create_fake_all_ct_segment_folder(all_ct_segment_folder, self.total_num_files)
        create_fake_avg_result_folder(avg_folder, 1)
        exp_result = list(map(lambda x: 'chr22_{}'.format(x), range(1, self.total_num_files)))
        replace_existing_file = 0
        obs_result = avg.get_genomic_positions_list(all_ct_segment_folder, avg_folder, replace_existing_file)
        self.assertCountEqual(exp_result, obs_result)
        replace_existing_file = 1
        exp_result = list(map(lambda x: 'chr22_{}'.format(x), range(self.total_num_files)))
        obs_result = avg.get_genomic_positions_list(all_ct_segment_folder, avg_folder, replace_existing_file)
        self.assertCountEqual(exp_result, obs_result)
        # after all the tests, we will remove all the fake folders we created to test the functions
        shutil.rmtree(all_ct_segment_folder)
        shutil.rmtree(avg_folder)
        return

if __name__ == "__main__":
    unittest.main()