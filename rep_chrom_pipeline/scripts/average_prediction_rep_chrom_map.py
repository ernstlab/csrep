import os
import sys
sys.path.append("/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/")
import helper
from subprocess import call # to call other python script like from command line arguments
import pandas as pd
import glob
import numpy as np

def get_genomic_positions_list(region_list_fn):
	gen_pos_list =  list(pd.read_csv(region_list_fn, sep = '\n', header = None)[0])
	return gen_pos_list 

def average_multiple_result_files(result_fn_list, output_fn):
	result_df_list = map(put_one_result_file_to_df, result_fn_list) # read all the files data and put them into  a data frame
	avg_df = pd.concat(result_df_list).groupby(level = 0).mean() # get the average across all the df. So what we get is a df : rows: genomic positions, columns: states, each entry is the average of the respective cells in all the input dfa
	avg_df.to_csv(output_fn, compression = 'gzip', header = True, index = True, sep = '\t') # save and compression to file
	return

def averaging_predictions_to_validate_one_ct(all_ct_prediction_folder, ct_list, outDir, num_chromHMM_state, gen_pos_list):
	# num_pred_ct: number of ct whose predictions we use to average out and get the predictions for the validate cell type
	# all_ct_prediction_folder = the directory where the data that are associated with the ct used for validation  cell type
	pred_dir_list = glob.glob(all_ct_prediction_folder + "/pred_*")
	pred_ct_list = list(map(lambda x: (x.split('/')[-1]).split('_')[-1], pred_dir_list)) # path/to/pred_E034 --> E034
	diff_set = (set(pred_ct_list)).difference(set(ct_list)) # difference betweeen the ct that we saw in the list of /pred_<ct> folders and the list of ct that we saw in the ct_list inputted
	assert len(diff_set) == 0, 'Number of pred_dir_list is not the same as number of specificed ct used to predict the model'
	# get the folder where the results of averaging across different prediction cell types will be stored
	validate_outDir = os.path.join(all_ct_prediction_folder, 'average_predictions')
	helper.make_dir(validate_outDir)
	for gene_window in gen_pos_list: # loop through each genomic window and then get the avrage result across different predictions for all positions in this window
		this_window_output_fn = os.path.join(validate_outDir, gene_window + "_avg_pred.txt.gz")
		this_window_pred_fn_list = map(lambda x: os.path.join(x, gene_window + "_pred_out.txt.gz"), pred_dir_list)
		# calculate the average prediction results across different prediction cell types for this window, and save the results
		average_multiple_result_files(this_window_pred_fn_list, this_window_output_fn)
	return 


def main():
	if len(sys.argv) != 6:
		usage()
	outDir = sys.argv[2]
	helper.make_dir(outDir)
	all_ct_prediction_folder = sys.argv[3]
	helper.check_dir_exist(all_ct_prediction_folder)
	num_chromHMM_state = helper.get_command_line_integer(sys.argv[4]) 
	all_ct_list_fn = sys.argv[5]
	helper.check_file_exist(all_ct_list_fn)
	region_list_fn = sys.argv[6]
	helper.check_file_exist(region_list_fn)
	print ("Done getting command line arguments")
	# get the list of all genomic positions used to segment the genome for our model training (we exclude chromosome Y in all analysis)
	gen_pos_list = get_genomic_positions_list(region_list_fn)
	# get all cell types
	ct_list = list(pd.read_csv(all_ct_list_fn, sep = '\n', header = None)[0]) # -->  a list with each entry being the cell type in this cell group
	# call all cell types
	averaging_predictions_to_validate_one_ct(all_ct_prediction_folder, ct_list, outDir, num_chromHMM_state, gen_pos_list)


def usage():
	print ("python run_cross_validation_pipeline.py")
	print ("outDir: Where all the output data (averaged predictions that can be used to represent predictions from all cell types) will be stored in the most structured form")
	print ("all_ct_prediction_folder: folder that contains lots of subfolders that contains the prediction data for all cell types are stored")
	print ("num_chromHMM_state: Num chromHMM states that we will predict the genomic positions upon")
	print ("all_ct_list_fn: each line of this file is the name of a cell type that is associated with the cell group that we are running the pipeline on. An example is /u/home/h/havu73/project-ernst/diff_pete/roadmap/blood.list")
	print ("region_list_fn: list of regions in the genome. This list corresponds to the files inside each prediction folder. Format of each line in this file is 'chr10_0' ")
	exit(1)

main()