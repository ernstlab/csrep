import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
import os
import sys
import time
import glob
import helper
import multiprocessing as mp
def get_X_colnames (train_cell_types, num_chromHMM_state):
	results = []
	for ct in train_cell_types:
		results += [ct + "_S" + str(x + 1) for x in range(num_chromHMM_state)]
	return results

def get_XY_segmentation_and_Cpos_data (train_cell_types, response_ct, num_chromHMM_state, train_posterior_data_fn):# given the segmentation data of multiple cell types, we want to extract information from the cell types that we will use as predictor (X) and response (Y). Predictor will be binarized with ct_state combinations as columns. Response will be just state labels 1 --> num_chromHMM_state. Rows of these two data frames correspond to different positions on the genome. Genome positions in all cell types' data are ordered exactly as in /u/home/h/havu73/project-ernst/diff_pete/roadmap/sample_genome_regions.gz
	input_Cpos_df = pd.read_csv(train_posterior_data_fn, header = 0, index_col = None, sep = '\t')
	X_colnames = get_X_colnames (train_cell_types, num_chromHMM_state)
	assert set(X_colnames).issubset(input_Cpos_df.columns), 'The colnames for train_cell_types are not as present as expected' # check that we have the columns that we need#
	Xtrain_cPos_df = input_Cpos_df[X_colnames]
	Y_colnames = list(map(lambda x: response_ct + '_S' + str(x+1), range(num_chromHMM_state)))
	assert set(Y_colnames).issubset(input_Cpos_df.columns), 'The colnames for response_ct are not as present as expected'
	Y_df = input_Cpos_df[Y_colnames] # now it has <response_ct>_S<state> --> we have to process it and figure out the state with the highest posterior probability that gets assigned at these places
	max_state_Y = Y_df.idxmax(axis = 1) # apply to each row --> the state that is maximum-posterior-probability for each genomic positions
	Y_return = max_state_Y.apply(lambda x: int((x.split('_')[1])[1:])) # convert E042_S15 to 15 as an integer
	return Xtrain_cPos_df, Y_return # X is binarized, Y is just state label 1 --> num_chromHMM_state

def train_multinomial_logistic_regression(X_df, Y_df, num_chromHMM_state):
	# give the Xtrain_cPos_df and Y_df obtained from get_XY_segmentation_data --> train a logistic regression object
	regression_machine = LogisticRegression(random_state = 0, solver = 'lbfgs', multi_class = 'multinomial', max_iter = 10000).fit(X_df, Y_df)
	return regression_machine 

def get_predictorX_chromHMM_pos_data(all_ct_chrom_pos_folder, chromosome, train_cell_types, num_chromHMM_state):
	# given the train_posterior_data_fn where all ct's posterior data are store, we want to out the files corresponding to the train_cell_types, and put them into a list so that we can average them out, which gives us the prediction for the chromatin state assignments for states across the chromsome. 
	# Return: a list of df, each corresponds to a train_cell_type. inside each df: rows: states, columns: regions on the chromosome that we are reading data in. 
	cPos_df_list = []
	one_train_ct_colnames = list(map(lambda x: 'state_' + str(x+1), range(num_chromHMM_state)))
	result_df = pd.DataFrame()
	for ct_index, ct in enumerate(train_cell_types):
		this_ct_cPos_fn = os.path.join(all_ct_chrom_pos_folder, ct, ct + '_18_core_K27ac_' + chromosome + '_posterior.txt.gz') # example: /u/home/h/havu73/project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded/posterior/E003/E003_18_core_K27ac_chr9_posterior.txt.gz
		this_ct_cPos_df = pd.read_csv(this_ct_cPos_fn, header = 0, skiprows = 1, index_col = None, sep = '\t')
		this_ct_cPos_df.columns = list(map(lambda x: ct + '_S' + str(x+1), range(num_chromHMM_state)))
		result_df[this_ct_cPos_df.columns] = this_ct_cPos_df
	return result_df

def save_predictor_df_into_multiple_files(response_df, chromosome, predict_outDir, num_chromHMM_state):
	total_segment = response_df.shape[0] # number of rows is the number of semgents here.
	num_output_fn = int(np.ceil(total_segment / helper.NUM_BIN_PER_WINDOW))
	print(num_output_fn)
	for window_index in range(num_output_fn - 1): #except the last window, which we will save after this
		output_fn = os.path.join(predict_outDir, chromosome + '_' + str(window_index) + '_pred_out.txt.gz')
		this_window_start_row_index = int(window_index * helper.NUM_BIN_PER_WINDOW)
		this_windw_end_row_index = int((window_index + 1) * helper.NUM_BIN_PER_WINDOW) 
		this_window_response_df = response_df.iloc[range(this_window_start_row_index, this_windw_end_row_index), :]
		this_window_response_df = this_window_response_df.reset_index(drop = True) #  now the index is always from 0 to the end of the window
		this_window_response_df.to_csv(output_fn, header = True, index = False, sep = '\t', compression = 'gzip')
		print("Done saving file: " + output_fn)
	last_output_fn = os.path.join(predict_outDir, chromosome + '_' + str(num_output_fn - 1) + '_pred_out.txt.gz')
	last_window_start_row_index = int((num_output_fn - 1) * helper.NUM_BIN_PER_WINDOW)
	last_window_end_row_index = int(response_df.shape[0]) # the last row
	last_window_response_df = response_df.iloc[range(last_window_start_row_index, last_window_end_row_index), :]
	last_window_response_df.to_csv(last_output_fn, header = True, index = False, sep = '\t', compression = 'gzip')
	print ('Done saving file: ' + last_output_fn)
	return 

def predict_segmentation_one_genomic_window(regression_machine, all_ct_chrom_pos_folder, chromosome, predict_outDir, train_cell_types, num_chromHMM_state):
	# based on the machine created through training, predict the segmentation corresponding to one specific window on the genome, specified in segment_fn. And print out, for each position, and for each chromHMM state, the probability that the region fall in to the state.
	# 1. Get the data of predictor cell types into the right format
	predictor_df = get_predictorX_chromHMM_pos_data(all_ct_chrom_pos_folder, chromosome, train_cell_types, num_chromHMM_state) # columns: predictor cell type - state combinations, rows: positions inside a  window on the genome
	# 2. Do the prediction job. Different model has different prediction functions
	response_df = regression_machine.predict_proba(predictor_df) # --> 2D: rows: positions (observations), columns: states (types) --> probability that each obs is of each type
	response_df = pd.DataFrame(response_df) # convert to a dataframe 
	response_df.columns = list(map(lambda x: 'state_' + str(x),regression_machine.classes_))
	# 3. Turn the results into readable format, then write to file. 
	save_predictor_df_into_multiple_files(response_df, chromosome, predict_outDir, num_chromHMM_state)
	print("Done getting files for: " + chromosome)
	return 

def one_job_run_predict_segmentation(regression_machine, all_ct_chrom_pos_folder, genome_pos_list, predict_outDir, train_cell_types, num_chromHMM_state):
	'''
	segment_fn_list and output_fn_list: the orders of regions in these two lists are similar (look at function predict_segmentation).
	Each element corresponds to a region on the genome
	'''
	for (window_index, chromosome) in enumerate(genome_pos_list):
		predict_segmentation_one_genomic_window(regression_machine, all_ct_chrom_pos_folder, chromosome, predict_outDir, train_cell_types, num_chromHMM_state)
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


def predict_segmentation (regression_machine, all_ct_chrom_pos_folder, predict_outDir, train_cell_types, response_ct, num_chromHMM_state):
	# 1. Get list of segmentation files corresponding to different windows on the genome.
	genome_pos_list = list(map(lambda x: 'chr' + x, helper.CHROMOSOME_LIST))
	# 2. partition the list of file names into groups, for later putting into jobs for multiple processes
	# num_cores = 2
	# partition_genome_pos_list = partition_file_list(genome_pos_list, num_cores)#genome_pos_list[:3], num_cores) # the input segment fn are actually by chromosome, because this file is all about chromHMM posterior data
	one_job_run_predict_segmentation(regression_machine, all_ct_chrom_pos_folder, genome_pos_list, predict_outDir, train_cell_types, num_chromHMM_state)
	# processes = [mp.Process(target = one_job_run_predict_segmentation, args = (regression_machine, all_ct_chrom_pos_folder,  partition_genome_pos_list[i], predict_outDir, train_cell_types, num_chromHMM_state)) for i in range(num_cores)] 
	# for p in processes:
	# 	p.start()
	# for i, p in enumerate(processes):
	# 	p.join()
	# 	print("Process " + str(i) + " is finished!")
	return 

def main():
	num_mandatory_args = 7
	if len(sys.argv) != num_mandatory_args:
		usage()
	print("training chromHMM_pos model")
	train_posterior_data_fn = sys.argv[1]
	helper.check_file_exist(train_posterior_data_fn)
	predict_outDir = sys.argv[2]
	helper.make_dir(predict_outDir)
	all_ct_chrom_pos_folder = sys.argv[3] # where we get data for the prediction after getting the regression machine
	helper.check_dir_exist(all_ct_chrom_pos_folder)
	try:
		num_chromHMM_state = int(sys.argv[4])
		assert num_chromHMM_state > 0, "num_chromHMM_state needs to be positive"
	except:
		print("num_chromHMM_state or num_train_ct is not valid")
		usage()
	response_ct = sys.argv[5]
	# remove the response ct from all_ct list
	train_ct_fn = sys.argv[6]
	helper.check_file_exist(train_ct_fn)
	all_ct = [line.strip() for line in open(train_ct_fn, "r").readlines()]
	try:
		all_ct.remove(response_ct)
	except:
		pass
	train_cell_types = all_ct # a mistake in the code, but let's keep it like this for now
	print("Done getting command line arguments train_chromHMM_pos_model.py")
	# 1. get the training data that we want
	Xtrain_cPos_df, Y_df = get_XY_segmentation_and_Cpos_data (train_cell_types, response_ct, num_chromHMM_state, train_posterior_data_fn)
	print("Done get_XY_segmentation_and_Cpos_data")
	# 2. Train the model based on the data that we got
	regression_machine = train_multinomial_logistic_regression(Xtrain_cPos_df, Y_df, num_chromHMM_state)
	print ("Done training")    
	# 2. no need to train model. Predict segmentation at each position 
	predict_segmentation (regression_machine, all_ct_chrom_pos_folder, predict_outDir, train_cell_types, response_ct, num_chromHMM_state)
	print("Done predicting whole genome")

def usage():
	print("train_chromHMM_posterior_model.py")
	print("train_posterior_data_fn: where chromHMM_pos data of all cell types  for this cell groups are stored, for the entire genome, so that we can get data for prediction out.")
	print("predict_outDir: where the beta values obtained from training the data will be stored (e.g. val_E###/average_predictions/")
	print("all_ct_chrom_pos_folder: where segmentation data of all cell types are stored, for the entire genome, so that we can get data for prediction out.")
	print("num_chromHMM_state: Number of chromHMM states that are shared across different cell types")
	print("response_ct : the cell type that we are trying to predict")
	print("train_ct_fn: filename of all ct being used in the model to produce the representative data of all these samples/ct")
	exit(1)
main()
