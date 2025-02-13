#!/usr/bin/env python
'''
This file will train a multi-variate logistic regression model to predict the chromatin state maps for a sample, based on the chromatin state maps in other samples. 

command-line argument: 
python train_multiLog_auto1Hot.py 
train_data_folder: where the state assignment and of training data are stored for all cell types. Each cell type has its own file
all_ct_segment_folder: where segmentation data of all cell types are stored, for the entire genome, so that we can get data for prediction out.
predict_outDir: where output data of the predictions of cell types are stored
response_ct: the cell type that we are trying to predict from the training dataset. This data is the Y value in our model training
num_chromHMM_state: Number of chromHMM states that are shared across different cell types
all_ct_fn: number of cell types that we will train
replace_existing_files: whether or not we would want to replace_existing_ output files 0 (no, only create result files for those that have not been outputted) or 1 (yes, rewrite everything)
seed: random seed for reproducibility
'''
import pandas as pd 
import numpy as np 
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import make_column_transformer
# import multiprocessing as mp
import os
import sys
import glob
import helper
import time
def get_X_colnames (train_cell_types, num_chromHMM_state):
	'''
	train_cell_types = ['E034', 'E037']
	num_chromHMM_state = 2
	--> ['E034_1', 'E034_2', 'E037_1', 'E037_2']
	'''
	results = []
	for ct in train_cell_types:
		results += list(map(lambda x: ct + "_" + str(x + 1), range(num_chromHMM_state)))
	return results

def get_state_number(state_annot):
	# 'E18' --> 18 (integers)
	try:
		return int(state_annot[1:])
	except:
		print('State number format in Y_df is not correct. You should recheck your input data. After preprocessing of data, all states should be of the form \'E<state_index_one_based>\', except for the case in training data where it could be that multiple states get mapped to the same training location. If the multi-state map problem occurs, it may be that your input data is not valid, because our assumption is that each position on the genome should be assigned to one single state. You can check this by looking at the bed file showing the regions  selected for training data, and bedtools intersect between the sampled regions and the input chromatin state map for the response_ct. If there are positins being mapped to multiple states, your data is invalid. If you cannot figure out, you may contact the authors: Ha Vu and Jason Ernst. Exitting ...')
		exit(1)

def transform_oneHot_trainCt_df(X_segment_df, train_cell_types, num_chromHMM_state, n_jobs):
	if not (X_segment_df.columns == train_cell_types).all():
		X_segment_df = X_segment_df[train_cell_types] # to make sure that we order the ct correctly, so that there will be matching columns as input features into the model for both training and for prediction
	all_states = np.array(list(map(lambda x: 'E{}'.format(x+1), range(num_chromHMM_state))))
	transformer = make_column_transformer((OneHotEncoder(categories = [all_states]*len(train_cell_types), handle_unknown = 'ignore'), train_cell_types), n_jobs = n_jobs, remainder='drop')
	X_data = transformer.fit_transform(X_segment_df)
	return X_data.toarray() # do not forget the .toarray function because I found out that the weights change slightly with and without the .toarray() function

def get_XY_segmentation_data (train_cell_types, response_ct, num_chromHMM_state, train_data_folder):
	# given the segmentation data of multiple cell types, we want to extract information from the cell types that we will use as predictor (X) and response (Y). Predictor will be binarized with ct_state combinations as columns. Response will be just state labels 1 --> num_chromHMM_state. Rows of these two data frames correspond to different positions on the genome. Genome positions in all cell types' data are ordered exactly as in /u/home/h/havu73/project-ernst/diff_pete/roadmap/sample_genome_regions.gz
	all_segment_df = pd.DataFrame()
	for ct in train_cell_types + [response_ct]:
		this_ct_fn = os.path.join(train_data_folder, ct + '_train_data.bed.gz') # file correponding to this X_ct
		this_ct_df = pd.read_csv(this_ct_fn, sep = '\t', header = 0) # open that file
		this_ct_df = this_ct_df[ct] # only pick columns that annotates the chromatin state for this cell type at each of those position
		all_segment_df = pd.merge(all_segment_df, this_ct_df, how = 'outer', left_index = True, right_index = True) # join columns, index-based. This is equivalent to a cbind in R
	all_segment_df = all_segment_df[all_segment_df.apply(lambda x: (~x.str.contains('[.,]', regex=True)))].dropna() # drop rows where in at least one cell type the state annnotation is either an empty match (.) or a multiple-state match (,). The multiple state match should not happen if the input data provided by users are directly learned from ChromHMM. However, in some cases, when the input annotations are actually lifted-Over from one ref.genome to another, it can happen that multiple states are maped to the same place. Usually, we want to get rid of those regions, but if the users forgot to do that, we will do that instead here for training data.
	Xtrain_segment_df = all_segment_df[train_cell_types]
	# after getting the state segmentations for all response cell types (E003 --> E127). Now we binarize the data: E003_S1 --> E003_S18 etc.
	print("Xtrain_segment_df before one hot encoding")
	n_jobs = 4
	Xtrain_segment_df = transform_oneHot_trainCt_df(Xtrain_segment_df, train_cell_types, num_chromHMM_state, n_jobs) #now , Xtrain_segment_df is actually a numpy array with 0 and 1 --> one-hot representation of states in training cell types
	print('Done getting binarized data for input sample annotations')
	Y_df = all_segment_df[response_ct]  
	Y_df = Y_df.apply(get_state_number)
	return Xtrain_segment_df, Y_df # X is binarized, Y is just state label 1 --> num_chromHMM_state


def get_predictorX_segmentation_data(train_cell_types, num_chromHMM_state, segment_fn):
	# given the segmentation data of multiple cell types, we want to extract information from the cell types that we will use as predictor (X). Predictor will be binarized with ct_state combinations as columns. Rows of this data frame correspond to different positions on the genome.
	segment_df = pd.read_csv(segment_fn, sep = '\t', header = 0)
	segment_df = segment_df[train_cell_types] # only get the data of the cell types that we need as predictors
	n_jobs = 4
	segment_df = transform_oneHot_trainCt_df(segment_df, train_cell_types, num_chromHMM_state, n_jobs)
	return segment_df

def train_multinomial_logistic_regression(X_df, Y_df, num_chromHMM_state, seed):
	# give the Xtrain_segment_df and Y_df obtained from get_XY_segmentation_data --> train a logistic regression object
	np.random.seed(seed)
	regression_machine = LogisticRegression(random_state = 0, solver = 'lbfgs', multi_class = 'multinomial', max_iter = 10000).fit(X_df, Y_df)
	print(regression_machine.coef_)
	return regression_machine 

def predict_segmentation_one_genomic_window(segment_fn, output_fn, train_cell_types, response_ct, num_chromHMM_state, regression_machine):
	# based on the machine created through training, predict the segmentation corresponding to one specific window on the genome, specified in segment_fn. And print out, for each position, and for each chromHMM state, the probability that the region fall in to the state.
	# 1. Get the data of predictor cell types into the right format
	predictor_df = get_predictorX_segmentation_data(train_cell_types, num_chromHMM_state, segment_fn) # rows: predictor cell type - state combinations, columns: positions inside a  window on the genome
	# 2. Do the prediction job. Different model has different prediction functions
	response_df = regression_machine.predict_proba(predictor_df) # --> 2D: rows: positions (observations), columns: states (types) --> probability that each obs is of each type
	response_df = pd.DataFrame(response_df) # convert to a dataframe 
	response_df.columns = regression_machine.classes_
	# soemtimes, due to the the training data not having observations of certain classes (states), therefore, we do not see that state included in the regression machine. Therefore, we assign the probabilities of missing states from the model to be 0. Usually, this is not an issue when we train using 10% of the genome.
	missing_states = np.setdiff1d(np.arange(1, num_chromHMM_state+1, 1), regression_machine.classes_) # missing states (not in the model)
	for state in missing_states:
		response_df[state] = 0 # fill up missing states with probabilities 0
	response_df = response_df[np.arange(1, num_chromHMM_state + 1, 1)] # rearrange the columns from 1 --> num_chromHMM_state
	# 3. Turn the results into readable format, then write to file. 
	response_df.columns = list(map(lambda x: "state_" + str(x + 1), range(num_chromHMM_state)))
	response_df.to_csv(output_fn, header = True, index = False, compression = 'gzip', sep = '\t')	
	print ("Done producing file: " + output_fn)


def one_job_run_predict_segmentation(segment_fn_list, output_fn_list, train_cell_types, response_ct, num_chromHMM_state, regression_machine):
	'''
	segment_fn_list and output_fn_list: the orders of regions in these two lists are similar (look at function predict_segmentation).  
	Each element corresponds to a region on the genome
	'''
	for (window_index, segment_fn) in enumerate(segment_fn_list):
		output_fn = output_fn_list[window_index]
		predict_segmentation_one_genomic_window(segment_fn, output_fn, train_cell_types, response_ct, num_chromHMM_state, regression_machine)
	return 

def find_uncalculated_gene_regions(predict_outDir, all_ct_segment_folder, replace_existing_files):
	"""
	predict_outDir: the output folder where output files <gene_region>_pred_out.txt.gz are stored
	all_ct_segment_folder: the folder of input data files for the training. 
	output: a list of genomic regions that were in the input folder, but missing in the output folder --> uncalculated regions of the genome
	"""
	calculated_fn_list = glob.glob(predict_outDir + '/*_pred_out.txt.gz')
	calculated_region_list = list(map(lambda x: x.split('/')[-1].split('_pred_out.txt.gz')[0], calculated_fn_list)) # chr<chrom>_<region_index>
	all_input_fn_list = glob.glob(all_ct_segment_folder + '/*_combined_segment.bed.gz')
	all_region_list = list(map(lambda x: x.split('/')[-1].split('_combined_segment.bed.gz')[0], all_input_fn_list))
	if replace_existing_files == 1: # the user wants to replace existing files in the output_folder, we willl return all the regions. If not, we return only regions whose output files are missing and hence need to be recalcualted
		return all_region_list
	not_calculated_region_list = list(np.setdiff1d(all_region_list, calculated_region_list))
	return not_calculated_region_list

def predict_segmentation (all_ct_segment_folder, regression_machine, predict_outDir, train_cell_types, response_ct, num_chromHMM_state, replace_existing_files):
	# 1. Get list of segmentation files corresponding to different windows on the genome.
	uncalculated_region_list = find_uncalculated_gene_regions(predict_outDir, all_ct_segment_folder, replace_existing_files) 
	segment_fn_list = list(map(lambda x: os.path.join(all_ct_segment_folder, x + '_combined_segment.bed.gz'), uncalculated_region_list))
	output_fn_list = list(map(lambda x: os.path.join(predict_outDir, x + "_pred_out.txt.gz"), uncalculated_region_list)) # get the output file names corresponding to different regions on the genome
	# 2. partition the list of file names into groups, for later putting into jobs for multiple processes
	one_job_run_predict_segmentation(segment_fn_list, output_fn_list, train_cell_types, response_ct, num_chromHMM_state, regression_machine)
	# num_cores = 4
	# partition_segment_fn_list = helper.partition_file_list(segment_fn_list, num_cores) # [process_index][file_index]
	# partition_output_fn_list = helper.partition_file_list(output_fn_list, num_cores)
	# processes = [mp.Process(target = one_job_run_predict_segmentation, args = (partition_segment_fn_list[i], partition_output_fn_list[i], train_cell_types, response_ct, num_chromHMM_state, regression_machine)) for i in range(num_cores)]
	# for p in processes:
	# 	p.start()
	# for i, p in enumerate(processes):
	# 	p.join()
	# 	print ("Process " + str(i) + " is finished!")
	
def get_train_cell_types(all_ct_fn, response_ct):
	# given the files that list all cell types of this cell groups, we would like to get the list of train cell types, which is the list of cell types that are not response_ct and are also listed in all_ct_fn
	ct_list =  list(pd.read_csv(all_ct_fn, sep = '\n', header = None)[0]) # -->  a list with each entry being the cell type in this cell group
	# get the index of the validation ct and then remove it from the list, so that we can focus this training process on ct other than the validation ct
	response_ct_index = ct_list.index(response_ct)
	assert response_ct_index != -1, "the validation ct is not present in the list of all cell type of the group that we are trying to train on"
	ct_list = ct_list[:response_ct_index] + ct_list[(response_ct_index + 1):] # skip the validation ct
	return ct_list

def main():
	start_time = time.time()
	num_mandatory_args = 9
	if len(sys.argv)!= num_mandatory_args: 
		usage()
	train_data_folder = sys.argv[1]
	helper.check_dir_exist(train_data_folder)
	all_ct_segment_folder = sys.argv[2] # where the genome-wide segmentation data of all cell types are combined, and stored in files corresponding to different regions in the genome.
	helper.check_dir_exist(all_ct_segment_folder)
	predict_outDir = sys.argv[3]
	helper.make_dir(predict_outDir)
	response_ct = sys.argv[4]
	try: 
		num_chromHMM_state = int(sys.argv[5])
		assert num_chromHMM_state > 0, "num_chromHMM_state needs to be positive"
	except:
		print ("num_chromHMM_state is not valid")
		usage()
	all_ct_fn = sys.argv[6]
	helper.check_file_exist(all_ct_fn)
	replace_existing_files = helper.get_command_line_integer(sys.argv[7])
	assert replace_existing_files in [0,1], 'replace_existing_files should be 0 (no, only create result files for those that have not been outputted) or 1 (yes, rewrite everything)'
	# get the list of train_cell_types as our training features
	seed = helper.get_command_line_integer(sys.argv[8])
	train_cell_types = get_train_cell_types(all_ct_fn, response_ct)
	print ("Done getting command line arguments")
	# 1. Get the data of predictors and response for training
	Xtrain_segment_df, Y_df = get_XY_segmentation_data (train_cell_types, response_ct, num_chromHMM_state, train_data_folder)
	end_time = time.time()
	print ("Done getting one hot data: {}".format(end_time - start_time))
	# 2. Get the regression machine
	regression_machine = train_multinomial_logistic_regression(Xtrain_segment_df, Y_df, num_chromHMM_state, seed)
	print(regression_machine.coef_)
	end_time = time.time()
	print ("Done training: {}".format(end_time - start_time))
	# 3. Based on the machine just created, process training data and then predict the segmentation at each position for the response_ct
	predict_segmentation (all_ct_segment_folder, regression_machine, predict_outDir, train_cell_types, response_ct, num_chromHMM_state,  replace_existing_files)
	end_time = time.time()
	print ("Done predicting whole genome: {}".format(end_time - start_time))
	
def usage():
	print ("python train_multiLog_auto1Hot.py ")
	print ("train_data_folder: where the state assignment and of training data are stored for all cell types. Each cell type has its own file")
	print ("all_ct_segment_folder: where segmentation data of all cell types are stored, for the entire genome, so that we can get data for prediction out.")
	print ("predict_outDir: where output data of the predictions of cell types are stored")
	print ("response_ct: the cell type that we are trying to predict from the training dataset. This data is the Y value in our model training")
	print ("num_chromHMM_state: Number of chromHMM states that are shared across different cell types")
	print ("all_ct_fn: number of cell types that we will train")
	print ("replace_existing_files: whether or not we would want to replace_existing_ output files 0 (no, only create result files for those that have not been outputted) or 1 (yes, rewrite everything)")
	print ('seed: random seed for reproducibility')
	exit(1)

if __name__ == '__main__':
	main()