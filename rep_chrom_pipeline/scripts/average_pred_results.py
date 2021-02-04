'''
This function will only train the cell types for one validation cell type. We need to call a bash script to create jobs for do multiple training jobs for multiple validation cell types
'''

import os, os.path
import sys
import helper
import pandas as pd
import glob
import numpy as np

def get_genomic_positions_list(all_ct_segment_folder, outDir, replace_existing_file):
	gen_pos_segment_fn_list = glob.glob(all_ct_segment_folder + '/chr*')
	gen_pos_list = [(x.split('/')[-1]).split('_combined_segment.bed.gz')[0] for x in gen_pos_segment_fn_list] # get the list of all genomic positions available: [chr9_11, chr9_12, etc.]
	if replace_existing_file == 1:
		return gen_pos_list
	else: # repace_existing_file is 0
		existing_fn_list = glob.glob(outDir + '*_avg_pred.txt.gz')
		existing_gen_pos_list = list(map(lambda x: (x.split('/')[-1]).split('_avg_pred.txt.gz')[0], existing_fn_list))
		gen_pos_list = np.setdiff1d(gen_pos_list, existing_gen_pos_list)
		print("These positions will be calculated in average_pred_results")
		print(gen_pos_list)
		return gen_pos_list

def put_one_result_file_to_df(fn):
	return pd.read_csv(fn, header = 0, sep = '\t')

def average_multiple_result_files(result_fn_list, output_fn):
	for fn in result_fn_list:
		if not os.path.isfile(fn):
			print('File: ' + fn + ' DOES NOT EXIST. The average of this region was not calculated')
			return
	result_df_list = list(map(put_one_result_file_to_df, result_fn_list[:2])) # read all the files data and put them into  a data frame
	avg_df = pd.concat(result_df_list).groupby(level = 0).mean() # get the average across all the df. So what we get is a df : rows: genomic positions, columns: states, each entry is the average of the respective cells in all the input dfa
	# now we normalize the resulting average probabilities, such that the row sum is always 0 (i,e, the probabilities of state assignments sum up to 1 over all states in a position)
	row_sum = avg_df.sum(axis = 1) # row sum, so that we can divide each entry in a row by the row sum corresponding to that row
	avg_df = avg_df.div(row_sum, axis = 0)
	avg_df.to_csv(output_fn, compression = 'gzip', header = True, index = False, sep = '\t') # save and compression to file
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

def averaging_prediction_one_process(gen_pos_list, outDir, pred_dir_list):
	for gene_window in gen_pos_list: # loop through each genomic window and then get the avrage result across different predictions for all positions in this window
		this_window_output_fn = os.path.join(outDir, gene_window + "_avg_pred.txt.gz")
		this_window_pred_fn_list = [os.path.join(x, gene_window + "_pred_out.txt.gz") for x in pred_dir_list]
		# calculate the average prediction results across different prediction cell types for this window, and save the results
		average_multiple_result_files(this_window_pred_fn_list, this_window_output_fn)
		print("Done averaging region: " + str(gene_window))
	return 

def averaging_predictions_to_validate_one_ct(outDir, all_ct_pred_folder, ct_list, gen_pos_list):
	# validate_ct_dir = the directory where the data that are associated with the ct used for validation  cell type
	pred_dir_list = glob.glob(all_ct_pred_folder + "/pred_*")
	pred_ct_list = [(x.split('/')[-1]).split('_')[-1] for x in pred_dir_list] # path/to/pred_E034 --> E034

	# comparing sets (union, interescetion...)
	assert set(pred_ct_list) == set(ct_list), 'pred_ct_list is not the same as ct_list'
	# get the folder where the results of averaging across different prediction cell types will be stored
	num_cores = 4
	partition_genome_pos_list = partition_file_list(gen_pos_list, num_cores)
	print(("Outputting average predictions here: " + outDir))
	processes = [mp.Process(target = averaging_prediction_one_process, args = (partition_genome_pos_list[i], outDir, pred_dir_list)) for i in range(num_cores)]
	for p in processes:
		p.start()
	for i, p in enumerate(processes):
		p.join()
		print("Process " + str(i) + " is finished!")
	return 



def main():
	if len(sys.argv) != 6:
		usage()
	outDir = sys.argv[1]
	helper.make_dir(outDir)
	all_ct_list_fn = sys.argv[2]
	all_ct_segment_folder = sys.argv[3]
	all_ct_pred_folder = sys.argv[4]
	replace_existing_file = helper.get_command_line_integer(sys.argv[5])
	assert replace_existing_file in range(2), 'get_command_line_integer can only be 0 or 1'
	print ("Done getting command line arguments in  average_pred_results.py")
	# get the list of all genomic positions used to segment the genome for our model training (we exclude chromosome Y in all analysis). This will only look at chrom X, deciding whether we want to rewrite some of the existing file
	gen_pos_list = get_genomic_positions_list(all_ct_segment_folder, outDir, replace_existing_file)
	# get all cell types
	ct_list =  list(pd.read_csv(all_ct_list_fn, sep = '\n', header = None)[0]) # -->  a list with each entry being the cell type in this cell group
	# call all cell types
	print ("Averaging pred results")
	averaging_predictions_to_validate_one_ct(outDir, all_ct_pred_folder, ct_list, gen_pos_list)


def usage():
	print( "python average_pred_results.py")
	print( "outDir: Where all the output data of averaging will be stored")
	print( "all_ct_list_fn: each line of this file is the name of a cell type that is associated with the cell group that we are running the pipeline on. An example is /u/home/h/havu73/project-ernst/diff_pete/roadmap/blood.list")
	print( "all_ct_segment_folder: folder with genomic positions for all ct")
	print( "all_ct_pred_folder: folder containing outdir and pred_E### folders")
	print("replace_existing_file: 0 or 1. 0--> we only calculate the average of regions where we have not calculated, 1 --> we replace all the files for all the genomic positions, no matter whether we have calculated them or not")
	exit(1)

main()
