#!/usr/bin/env python

'''
This script will take in CSREP or base_count's summary chromatin state map (state assignment probabilities) for two groups of samples, and subtract group2's summary state maps from group1's (group1 - group2). The resulting output is a matrix showing the differential chromatin scores, rows: genomic bins, columns: states


command-line arguments:
python calculate_diff_rep_state_map.py
group1_folder: where the representative state maps for group1 are stored. The files in this fodler correspond to different regions in the genome
group2_folder: where the representative state maps for group2 are stored. The files in this fodler correspond to different regions in the genome
output_folder: where the difference between the representative state maps of two groups are stored. The files in this folder correspon to different regions in the genome
num_chromHMM_state
redo_existing_file: 0 or 1: 1 (yes, rewrite all the existing files in the output_folder, or 0 (no, only write files that have not been produced)
Output is a matrix, rows: genomic bins, columns: states
'''

import pandas as pd 
import numpy as np 
import sys
import os
import glob
import helper
import multiprocessing as mp
NUM_CORES = 4


def read_rep_df(fn, num_chromHMM_state):
	right_colnames = list(map(lambda x: 'state_' + str(x+1), range(num_chromHMM_state)))
	helper.check_file_exist(fn)
	try:
		df = pd.read_csv(fn, header = 0, index_col = 0, sep = '\t')
		assert list(df.columns) == right_colnames, 'index_col = 0 is wrong'
	except:
		df = pd.read_csv(fn, header = 0, index_col = None, sep = '\t')
		assert list(df.columns) == right_colnames, 'index_col = 0 is wrong'
	return df

def get_diff_rep_state_one_genomic_region(region_fn, group1_folder, group2_folder, output_folder, num_chromHMM_state):
	# this function will take data of representative chromatin maps in one regions in the genome (usually 10M bp), for both group1 and group2
	# this function should return a dataframes of rows: genomic region, columns: states, values: differential chromatin state assigment probabilities between the two regions
	g1_fn = os.path.join(group1_folder, region_fn)
	g2_fn = os.path.join(group2_folder, region_fn)
	helper.check_file_exist(g1_fn)
	helper.check_file_exist(g2_fn)
	g1_df = read_rep_df(g1_fn, num_chromHMM_state)
	g2_df = read_rep_df(g2_fn, num_chromHMM_state)
	diff_df = g1_df - g2_df # the difference in chromatin state assignment probabilities between the two groups
	output_fn = os.path.join(output_folder, region_fn)
	diff_df.to_csv(output_fn, header = True, index = True, sep = '\t', compression = 'gzip')
	return 

def get_diff_rep_state_one_process(region_fn_list, group1_folder, group2_folder, output_folder, num_chromHMM_state):
	# this function will pass through each of the file in region_fn_list and calculate the diff chromatin state map for each region, then write into output_folder
	for region_fn in region_fn_list:
		get_diff_rep_state_one_genomic_region(region_fn, group1_folder, group2_folder, output_folder, num_chromHMM_state)
	return 

def get_genomic_positions_list(group1_folder, output_folder, redo_existing_file):
	gen_pos_segment_fn_list = glob.glob(group1_folder + '/chr*')
	gen_pos_list = [(x.split('/')[-1]).split('_avg_pred.txt.gz')[0] for x in gen_pos_segment_fn_list] # get the list of all genomic positions available: [chr9_11, chr9_12, etc.]
	if redo_existing_file == 1:
		return gen_pos_list
	else: # repace_existing_file is 0
		existing_fn_list = glob.glob(output_folder + '/*_avg_pred.txt.gz')
		existing_gen_pos_list = list(map(lambda x: (x.split('/')[-1]).split('_avg_pred.txt.gz')[0], existing_fn_list))
		gen_pos_list = np.setdiff1d(gen_pos_list, existing_gen_pos_list)
		print("Number of positions that will be calculated in average_pred_results: " + str(len(gen_pos_list)))
		print(gen_pos_list)
		return gen_pos_list

def get_diff_rep_state_whole_genome(group1_folder, group2_folder, output_folder, num_chromHMM_state, redo_existing_file):
	# this function will:
	# 1. get list of regions(files) inside folder diff_folder, partition into NUM_CORES processes
	g1_fn_list = os.listdir(group1_folder)
	g2_fn_list = os.listdir(group2_folder)
	print(len(g1_fn_list))
	print(len(g2_fn_list))
	assert len(g1_fn_list) == len(g2_fn_list), 'Number of files from group 1 and group 2 are not equal'
	gen_pos_list = get_genomic_positions_list(group1_folder, output_folder, redo_existing_file) # the list of regions that we still need to calculate the differential scores for
	input_fn_list = list(map(lambda x: x + '_avg_pred.txt.gz', gen_pos_list))
	print('Number of files to calculate differential scores for: {}'.format(len(input_fn_list)))
	print(gen_pos_list)
	partition_fn_list = helper.partition_file_list(input_fn_list, NUM_CORES) # partition the list of fn in group1 (list of genomeic regions) into subgroups, so that we can add them into processes later
	# 2. For each group of region_fn we will call a process to calculate the difference between two groups' chromatin state assignment
	processes = [mp.Process(target = get_diff_rep_state_one_process, args = (partition_fn_list[i], group1_folder, group2_folder, output_folder, num_chromHMM_state)) for i in range(NUM_CORES)]
	for p in processes:
		p.start()
	for i, p in enumerate(processes):
		p.join()
		print ("Process " + str(i) + " is finished!")
	return

def main():
	if len(sys.argv) != 6:
		usage()
	group1_folder = sys.argv[1] # where the representative state maps for group1 are stored. The files in this fodler correspond to different regions in the genome
	helper.check_dir_exist(group1_folder)
	group2_folder = sys.argv[2] # where the representative state maps for group2 are stored. The files in this fodler correspond to different regions in the genome
	helper.check_dir_exist(group2_folder)
	output_folder = sys.argv[3]
	helper.make_dir(output_folder)
	num_chromHMM_state = helper.get_command_line_integer(sys.argv[4])
	redo_existing_file = helper.get_command_line_integer(sys.argv[5])
	assert redo_existing_file in [0,1], 'redo_existing_file should be 1 (yes, rewrite all the existing files in the output_folder, or 0 (no, only write files that have not been produced)'
	print ("Done getting command line argument")
	# call the function that will do the necessary work for this code
	get_diff_rep_state_whole_genome(group1_folder, group2_folder, output_folder, num_chromHMM_state, redo_existing_file)
	print ('Done!')

def usage():
	print ("python calculate_diff_rep_state_map.py")
	print ("group1_folder: where the representative state maps for group1 are stored. The files in this fodler correspond to different regions in the genome")
	print ("group2_folder: where the representative state maps for group2 are stored. The files in this fodler correspond to different regions in the genome")
	print ("output_folder: where the difference between the representative state maps of two groups are stored. The files in this folder correspon to different regions in the genome")
	print ('num_chromHMM_state')
	print ('redo_existing_file: 0 or 1: 1 (yes, rewrite all the existing files in the output_folder, or 0 (no, only write files that have not been produced)')
	print ("This code will calculate the difference of chromatin state assignment probabilities between two groups. The two groups should already had their representative chromatin state maps being calcualted. Output is a matrix, rows: genomic bins, columns: states")
	exit(1)

if __name__ == '__main__':
	main()
