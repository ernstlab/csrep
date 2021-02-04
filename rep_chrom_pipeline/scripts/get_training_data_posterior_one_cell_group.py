import pandas as pd
import numpy as np
import os
import sys
import glob
import helper
import multiprocessing as mp

def get_training_data_posterior_one_cell_group_one_genomic_window(all_ct_chrom_pos_folder, chromosome, train_cell_types, num_chromHMM_state, this_region_sample_pos_df):
	this_region_sample_pos_df['bin_index_from_chrom'] = this_region_sample_pos_df['start_bp'] / helper.NUM_BP_PER_BIN
	this_region_sample_pos_df['bin_index_from_chrom'] = (this_region_sample_pos_df['bin_index_from_chrom']).astype(int)
	input_fn_list = list(map(lambda x: os.path.join(all_ct_chrom_pos_folder, x, x + '_18_core_K27ac_' + chromosome + '_posterior.txt.gz'), train_cell_types))# get the file names of files that contain the posterior data for the listeed cell types and the listed chromosome
	result_df = pd.DataFrame(columns = ['chrom', 'start_bp', 'end_bp'])
	result_df[result_df.columns] = this_region_sample_pos_df[['chrom', 'start_bp', 'end_bp']]
	for ct_index, ct in enumerate( train_cell_types):
		this_ct_input_fn = input_fn_list[ct_index]
		this_ct_input_df = pd.read_csv(this_ct_input_fn, sep = '\t', skiprows = 1, header = 0, index_col = None)
		this_ct_input_df = this_ct_input_df.iloc[this_region_sample_pos_df['bin_index_from_chrom'], :] # only get data from the positions that we have decided that will contribute to the training data
		this_ct_input_df.reset_index(drop = True, inplace = True)
		this_ct_input_df.columns = list(map(lambda x: ct + '_S' + str(x[1:]), list(this_ct_input_df.columns)))
		result_df[this_ct_input_df.columns] = this_ct_input_df
	return result_df

def one_job_run_get_training_data_posterior_one_cell_group(all_ct_chrom_pos_folder, genome_pos_list, output_fn, train_cell_types, num_chromHMM_state, training_sample_pos_df):
	'''
	segment_fn_list and output_fn_list: the orders of regions in these two lists are similar (look at function get_training_data_posterior_one_cell_group).
	Each element corresponds to a region on the genome
	'''
	result_df = pd.DataFrame()
	for (window_index, chromosome) in enumerate(genome_pos_list): 
		this_region_sample_pos_df = training_sample_pos_df.get_group(chromosome)
		this_region_sample_pos_df.reset_index(inplace = True, drop = True)
		this_region_out_df = get_training_data_posterior_one_cell_group_one_genomic_window(all_ct_chrom_pos_folder, chromosome, train_cell_types, num_chromHMM_state, this_region_sample_pos_df)
		result_df = result_df.append(this_region_out_df)
	result_df.to_csv(output_fn, header = True, index = False, sep = '\t')
	return

def partition_file_list(file_list, num_cores):
    results = [] # list of lists of file names
    num_files_per_core = int(len(file_list) / num_cores)
    for core_i in range(num_cores):
        if core_i < (num_cores - 1):
            this_core_files = file_list[core_i * num_files_per_core : (core_i + 1) * num_files_per_core]
        elif core_i == (num_cores - 1):
            this_core_files = file_list[core_i * num_files_per_core :]
        results.append(this_core_files)
    return results

def join_output_from_processes(train_data_outDir, num_cores):
	draft_fn_list = list(map(lambda x: os.path.join(train_data_outDir, 'draft_' + str(x+1) + '.txt'), range(num_cores)))
	result_df = pd.DataFrame()
	for core_i, draft_fn in enumerate(draft_fn_list):
		draft_df = pd.read_csv(draft_fn, sep = '\t', header = 0, index_col = None)
		result_df = result_df.append(draft_df)
		os.remove(draft_fn)
	output_fn = os.path.join(train_data_outDir, 'training_data_chromHMM_posterior.txt.gz')
	result_df.to_csv(output_fn, header = True, index = False, sep = '\t')
	return result_df

def get_training_data_posterior_one_cell_group (all_ct_chrom_pos_folder, train_data_outDir, train_cell_types, num_chromHMM_state, training_sample_pos_df):
	# 1. Get list of segmentation files corresponding to different windows on the genome.
	genome_pos_list = list(map(lambda x: 'chr' + str(x), helper.CHROMOSOME_LIST)) # chr1 --> chrX, chrY is excluded
	# 2. partition the list of file names into groups, for later putting into jobs for multiple processes
	num_cores = 4
	partition_genome_pos_list = partition_file_list(genome_pos_list, num_cores)#genome_pos_list[:3], num_cores) # the input segment fn are actually by chromosome, because this file is all about chromHMM posterior data
	draft_fn_list = list(map(lambda x: os.path.join(train_data_outDir, 'draft_' + str(x+1) + '.txt'), range(num_cores)))
	processes = [mp.Process(target = one_job_run_get_training_data_posterior_one_cell_group, args = (all_ct_chrom_pos_folder,  partition_genome_pos_list[i], draft_fn_list[i], train_cell_types, num_chromHMM_state, training_sample_pos_df)) for i in range(num_cores)] 
	for p in processes:
		p.start()
	for i, p in enumerate(processes):
		p.join()
		print("Process " + str(i) + " is finished!")
	result_df = join_output_from_processes(train_data_outDir, num_cores)
	print("Done combining the data for all processes!")
	return

def get_training_sample_pos_df(training_sample_pos_fn):
	training_sample_pos_df = pd.read_csv(training_sample_pos_fn, header = None, index_col = None, sep = '\t')
	training_sample_pos_df.columns = ['chrom', 'start_bp', 'end_bp']
	training_sample_pos_df = training_sample_pos_df.groupby('chrom')
	return (training_sample_pos_df)

def main():
	num_mandatory_args = 6
	if len(sys.argv) != num_mandatory_args:
		usage()
	print("training chromHMM_pos model")
	all_ct_chrom_pos_folder = sys.argv[1]
	helper.check_dir_exist(all_ct_chrom_pos_folder)
	train_data_outDir = sys.argv[2]
	helper.make_dir(train_data_outDir)
	try:
		num_chromHMM_state = int(sys.argv[3])
		assert num_chromHMM_state > 0, "num_chromHMM_state needs to be positive"
	except:
		print("num_chromHMM_state or num_train_ct is not valid")
		usage()
	# remove the response ct from all_ct list
	all_ct_fn = sys.argv[4]
	all_ct = [line.strip() for line in open(all_ct_fn, "r").readlines()]
	train_cell_types = all_ct
	training_sample_pos_fn = sys.argv[5]
	helper.check_file_exist(training_sample_pos_fn)
	training_sample_pos_df = get_training_sample_pos_df(training_sample_pos_fn)
	print("Done getting command line arguments train_chromHMM_pos_model.py")	
	# 2. no need to train model
	# 3. process training data and predict segmentation at each position 
	get_training_data_posterior_one_cell_group (all_ct_chrom_pos_folder, train_data_outDir, train_cell_types, num_chromHMM_state, training_sample_pos_df)
	print("Done predicting whole genome")

def usage():
	print("python get_training_data_posterior_one_cell_group.py")
	print("all_ct_chrom_pos_folder: where chromHMM posterior data of all cell types are stored, for the entire genome. We will take a portion of this data folder out into the output")
	print("train_data_outDir: where the data of training data for this cell group is stored")
	print("num_chromHMM_state: Number of chromHMM states that are shared across different cell types")
	print("all_ct_fn: filename of all ct being used in the model to produce the representative data of all these samples/ct")
	print("training_sample_pos_df: where the position of the training data are stored")
	exit(1)
main()
