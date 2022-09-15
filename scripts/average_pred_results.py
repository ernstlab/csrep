#!/usr/bin/env python

'''
This script will take in the CSREP predictions of chrom-state assignment probabilites for individual samples, average the predictions across samples, and print out the outputs as files. 
This script will do the averaging of predictions in parallel, hence speeding up the computing time. 

command-line arguments:
python average_pred_results.py
outDir: Where all the output data of averaging will be stored
all_ct_list_fn: each line of this file is the name of a cell type that is associated with the cell group that we are running the pipeline on. An example is /u/home/h/havu73/project-ernst/diff_pete/roadmap/blood.list
all_ct_segment_folder: folder with genomic positions for all ct
all_ct_pred_folder: folder containing outdir and pred_E### folders
replace_existing_file: 0 or 1. 0--> we only calculate the average of regions where we have not calculated, 1 --> we replace all the files for all the genomic positions, no matter whether we have calculated them or not
num_chromHMM_state
'''

import os, os.path
import sys
import helper
import pandas as pd
import glob
import numpy as np
import multiprocessing as mp
def get_genomic_positions_list(all_ct_segment_folder, outDir, replace_existing_file):
	"""
	Function to get a list of genomic regions (ex: chr1_0 as the first 10MB in chrom 1) for which we will calculate the CSREP summary chromatin state maps for.

	Args:
	    all_ct_segment_folder: folder with chrom state  data of all input samples, each file in this folder represents input data for one genomic region (<= 10Mb)/
	    outDir: output_folder of the summary chromatin state maap. Each file in the folder corresponds to 1 genomic region (<=10Mb). 
	    replace_existing_file: 0/1. If there are existing output files in outDir, whether or not (0 / 1) we want to replace them by recalculating the summary chromatin state map for the corresponding genomic regions

	Returns:
		A list of strings: genomic regions (example: chro1_0) for which we will calculate the summary chromatin state map for

	Raises:
	    KeyError: Raises an exception: No exceptions. 
	"""	
	gen_pos_segment_fn_list = glob.glob(all_ct_segment_folder + '/chr*')
	gen_pos_list = [(x.split('/')[-1]).split('_combined_segment.bed.gz')[0] for x in gen_pos_segment_fn_list] # get the list of all genomic positions available: [chr9_11, chr9_12, etc.]
	if replace_existing_file == 1:
		return gen_pos_list
	else: # repace_existing_file is 0
		existing_fn_list = glob.glob(outDir + '/*_avg_pred.txt.gz')
		existing_gen_pos_list = list(map(lambda x: (x.split('/')[-1]).split('_avg_pred.txt.gz')[0], existing_fn_list))
		gen_pos_list = np.setdiff1d(gen_pos_list, existing_gen_pos_list)
		print("Number of positions that will be calculated in average_pred_results: " + str(len(gen_pos_list)))
		print(gen_pos_list)
		return gen_pos_list


def read_rep_df(fn, num_chromHMM_state):
	"""
	Read a file that shows predictions of chromatin state assignment probabilities into a dataframe. 
	Args:
	    fn: filename to open and read. 
	    num_chromHMM_state: number of chromatin states in the model

	Returns:
		df: a dataframe, read in from fn

	Raises:
	    KeyError: Raises an exception: No exceptions. 
	"""	
	right_colnames = list(map(lambda x: 'state_' + str(x+1), range(num_chromHMM_state)))
	helper.check_file_exist(fn)
	try:
		df = pd.read_csv(fn, header = 0, index_col = 0, sep = '\t')
		assert list(df.columns) == right_colnames, 'index_col = 0 is wrong'
	except:
		df = pd.read_csv(fn, header = 0, index_col = None, sep = '\t')
		assert list(df.columns) == right_colnames, 'index_col = 0 is wrong'
	return df

def average_multiple_result_files(result_fn_list, output_fn, num_chromHMM_state):
	"""
	Given multiple files showing the predictions of chromatin state assignment probabilties in multiple samples, this function will calculate the average of all predictions
	Args:
	    result_fn_list: list of filename showing predictions of chrom state assignment probabilities in multiple samples. 
	    output_fn: fn showing the average chrom-state assignment probabilities
	    num_chromHMM_state: number of chromatin states in the model

	Returns:
		None, the function just do average operations and save the results into output_fn

	Raises:
	    KeyError: Raises an exception: No exceptions. 
	"""		
	for fn in result_fn_list:
		if not os.path.isfile(fn):
			print('File: ' + fn + ' DOES NOT EXIST. The average of this region was not calculated')
			return
	result_df_list = list(map(lambda x: read_rep_df(x, num_chromHMM_state), result_fn_list)) # read all the files data and put them into  a data frame
	avg_df = pd.concat(result_df_list).groupby(level = 0).mean() # get the average across all the df. So what we get is a df : rows: genomic positions, columns: states, each entry is the average of the respective cells in all the input dfa
	# avg_df.to_csv(output_fn, compression = 'gzip', header = True, index = True, sep = '\t') # save and compression to file
	# now we normalize the resulting average probabilities, such that the row sum is always 0 (i,e, the probabilities of state assignments sum up to 1 over all states in a position)
	row_sum = avg_df.sum(axis = 1) # row sum, so that we can divide each entry in a row by the row sum corresponding to that row
	row_norm_avg_df = avg_df.div(row_sum, axis = 0) # Avg_df should already be row-normarlized, meaning the sum of values in each row is 1.0, but we still keep this code here
	row_norm_avg_df.to_csv(output_fn, compression = 'gzip', header = True, index = False, sep = '\t') # save and compression to file
	return


def averaging_prediction_one_process(gen_pos_list, outDir, pred_dir_list, num_chromHMM_state):
	"""
	This function will call on function average_multiple_result_files multiple times to average the predictions of chrom-state-assignment probs across samples for multiple regions of the genome
	Args:
	    gen_pos_list: list of genomic positions (ex: chr1_0 as the first 10Mb region of chrom1) to do generate average results for. 
	    outDir: outptu folder containing output files for multiple genomic regions
	    pred_dir_list: list of folder paths that contain the predictions results for multiple samples
	    num_chromHMM_state: number of chromatin state in the model

	Returns:
		None, the function only calls on average_multiple_result_files to print out prediction output into files 

	Raises:
	    KeyError: Raises an exception: No exceptions. 
	"""		
	for gene_window in gen_pos_list: # loop through each genomic window and then get the avrage result across different predictions for all positions in this window
		this_window_output_fn = os.path.join(outDir, gene_window + "_avg_pred.txt.gz")
		this_window_pred_fn_list = [os.path.join(x, gene_window + "_pred_out.txt.gz") for x in pred_dir_list]
		# calculate the average prediction results across different prediction cell types for this window, and save the results
		average_multiple_result_files(this_window_pred_fn_list, this_window_output_fn, num_chromHMM_state)
		print("Done averaging region: " + str(gene_window))
	return 

def averaging_predictions_all_processes(outDir, all_ct_pred_folder, ct_list, gen_pos_list, num_chromHMM_state):
	"""
	This function will call on function averaging_prediction_one_process to do averaging of predicted chrom-state-assignment probabilities across samples for multiple processes in parallel --> speed up the averaging task by doing it in parallel for multiple genomic regions at a time
	Args:
	    outDir: outptu folder containing output files for multiple genomic regions
	    all_ct_pred_folder: folder where there are subfolders pred_<ct_index> --> each subfolder contains the prediction results for one sample
	    ct_list: list of ct_index (or sample ID, for example for Roadmap data: E003). 
	    gen_pos_list: list of genomic positions (ex: chr1_0 as the first 10Mb region of chrom1) to do generate average results for. 
	    num_chromHMM_state: number of chromatin state in the model

	Returns:
		None, the function only calls on averaging_prediction_one_process to print out prediction output into files 
			
	Raises:
	    KeyError: Raises an exception: No exceptions. 
	"""			
	# validate_ct_dir = the directory where the data that are associated with the ct used for validation  cell type
	pred_dir_list = glob.glob(all_ct_pred_folder + "/pred_*")
	pred_ct_list = [(x.split('/')[-1]).split('_')[-1] for x in pred_dir_list] # path/to/pred_E034 --> E034

	# comparing sets (union, interescetion...)
	assert set(pred_ct_list) == set(ct_list), 'pred_ct_list is not the same as ct_list'
	# get the folder where the results of averaging across different prediction cell types will be stored
	num_cores = 4
	partition_genome_pos_list = helper.partition_file_list(gen_pos_list, num_cores)
	print(("Outputting average predictions here: " + outDir))
	processes = [mp.Process(target = averaging_prediction_one_process, args = (partition_genome_pos_list[i], outDir, pred_dir_list, num_chromHMM_state)) for i in range(num_cores)]
	for p in processes:
		p.start()
	for i, p in enumerate(processes):
		p.join()
		print("Process " + str(i) + " is finished!")
	return 



def main():
	if len(sys.argv) != 7:
		usage()
	outDir = sys.argv[1]
	helper.make_dir(outDir)
	all_ct_list_fn = sys.argv[2]
	all_ct_segment_folder = sys.argv[3]
	all_ct_pred_folder = sys.argv[4]
	replace_existing_file = helper.get_command_line_integer(sys.argv[5])
	num_chromHMM_state = helper.get_command_line_integer(sys.argv[6])
	assert replace_existing_file in range(2), 'get_command_line_integer can only be 0 or 1'
	print ("Done getting command line arguments in  average_pred_results.py")
	# get the list of all genomic positions used to segment the genome for our model training (we exclude chromosome Y in all analysis). This will only look at chrom X, deciding whether we want to rewrite some of the existing file
	gen_pos_list = get_genomic_positions_list(all_ct_segment_folder, outDir, replace_existing_file)
	# get all cell types
	ct_list = list(pd.read_csv(all_ct_list_fn, sep = '\n', header = None)[0]) # -->  a list with each entry being the cell type in this cell group
	# call all cell types
	print ("Averaging pred results")
	averaging_predictions_all_processes(outDir, all_ct_pred_folder, ct_list, gen_pos_list, num_chromHMM_state)


def usage():
	print( "python average_pred_results.py")
	print( "outDir: Where all the output data of averaging will be stored")
	print( "all_ct_list_fn: each line of this file is the name of a cell type that is associated with the cell group that we are running the pipeline on. An example is /u/home/h/havu73/project-ernst/diff_pete/roadmap/blood.list")
	print( "all_ct_segment_folder: folder with genomic positions for all ct")
	print( "all_ct_pred_folder: folder containing outdir and pred_E### folders")
	print("replace_existing_file: 0 or 1. 0--> we only calculate the average of regions where we have not calculated, 1 --> we replace all the files for all the genomic positions, no matter whether we have calculated them or not")
	print('num_chromHMM_state')
	exit(1)

if __name__ == '__main__':
	main()
