import pandas as pd
import numpy as np
import os
import sys
import time
import glob
import helper
import multiprocessing as mp

def get_X_colnames (train_cell_types, num_chromHMM_state):
        '''
        train_cell_types = ['E034', 'E037']
        num_chromHMM_state = 2
        --> ['E034_1', 'E034_2', 'E037_1', 'E037_2']
        '''
        results = []
        for ct in train_cell_types:
                results += [ct + "_" + str(x + 1) for x in range(num_chromHMM_state)]
        return results

def get_state_number(state_annot):
        # 'E18' --> 18 (integers)
        return int(state_annot[1:])

def get_binary_state_assignment(state_train_segment, num_chromHMM_state):
        '''
        a series of state assignment with indices ['E034', 'E037'] --> [E1, E2]
        num_chromHMM_state = 2
        --> a series of 0/1 with indices ['E034_1', 'E034_2', 'E037_1', 'E037_2'] --> [1,0,0,1]
        '''
        results = [0] * (len(state_train_segment) * num_chromHMM_state)
        for ct_index, train_segment in enumerate(state_train_segment):
                train_segment = train_segment[1:] # get rid of the E as a prefix
                index_to_put_one = (int(train_segment) - 1) + ct_index * num_chromHMM_state
                results[index_to_put_one] = 1
        return pd.Series(results)


def get_predictorX_segmentation_data(train_cell_types, num_chromHMM_state, segment_fn):
        # given the segmentation data of multiple cell types, we want to extract information from the cell types that we will use as predictor (X). Predictor will be binarized with ct_state combinations as columns. Rows of this data frame correspond to different positions on the genome.
        segment_df = pd.read_csv(segment_fn, sep = '\t', header = 0)
        segment_df = segment_df[train_cell_types] # only get the data of the cell types that we need as predictors
        #segment_df = segment_df.applymap(get_state_number)# get the state train_segmentation from 'E18' --> 18 (integers)
        segment_df = segment_df.apply(lambda x: get_binary_state_assignment(x, num_chromHMM_state), axis = 1) # apply function row-wise, change from the state train_segmentation to binarized of state
        segment_df.columns = get_X_colnames(train_cell_types, num_chromHMM_state)
        return segment_df

def predict_baseline_segmentation(predictor_df, num_chromHMM_state):
        # predictor_df : <ct>_state<state_index> --> each columns corresponds to 0/1 based on segmentation of that specific cell type
        # rows: genomic positions
        # return response_df: rows: genomic pos, columns: state<state_index> --> probability that each pos is in each state
        num_train_ct = int(len(predictor_df.columns) / num_chromHMM_state)
        train_df_list = []
        state_colnames = ["state_" + str(x+1) for x in range(num_chromHMM_state)] # ex: [state_1, state_2, ..., state_18]
        for ct_index in range(num_train_ct): # get data for each ct that we use to train the model
            start_colname_index = ct_index * num_chromHMM_state # starting from <ct>_state1
            end_colname_index = (ct_index + 1) * num_chromHMM_state # till <ct>_state<num_chromHMM_state>
            this_ct_train_df = predictor_df[predictor_df.columns[start_colname_index:end_colname_index]] # get the data corresponding to states in this particular cell type
            this_ct_train_df.columns = state_colnames
            train_df_list.append(this_ct_train_df)
        # Now we have obtained the training data for each of the cell types, we will average out the one-hot encoding of the state assignment across cell types. That's the basis of the baseline training where basically the state assigned to each position is the state that have the highest count of cell types where this position is assigned to.
        response_df = pd.concat(train_df_list).groupby(level = 0).mean() # get the average across all the df. So what we get is a df : rows: genomic positions, columns: states, each entry is the average of the respective cells in all the input dfa
        return response_df

def predict_segmentation_one_genomic_window(segment_fn, output_fn, train_cell_types, num_chromHMM_state, train_mode):
        # based on the machine created through training, predict the segmentation corresponding to one specific window on the genome, specified in segment_fn. And print out, for each position, and for each chromHMM state, the probability that the region fall in to the state.
        # 1. Get the data of predictor cell types into the right format
        predictor_df = get_predictorX_segmentation_data(train_cell_types, num_chromHMM_state, segment_fn) # rows: predictor cell type - state combinations, columns: positions inside a  window on the genome
        # 2. Do the prediction job. Different model has different prediction functions
        response_df = predict_baseline_segmentation(predictor_df, num_chromHMM_state)
        response_df.to_csv(output_fn, header = True, index = False, compression = 'gzip', sep = '\t')
        print("Done producing file: " + output_fn)

def one_job_run_predict_segmentation(segment_fn_list, output_fn_list, train_cell_types, num_chromHMM_state, train_mode):
        '''
        segment_fn_list and output_fn_list: the orders of regions in these two lists are similar (look at function predict_segmentation).
        Each element corresponds to a region on the genome
        '''
        for (window_index, segment_fn) in enumerate(segment_fn_list):
                output_fn = output_fn_list[window_index]
                predict_segmentation_one_genomic_window(segment_fn, output_fn, train_cell_types, num_chromHMM_state, train_mode)
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

def predict_segmentation (all_ct_segment_folder, predict_outDir, train_cell_types, num_chromHMM_state, train_mode):
        # 1. Get list of segmentation files corresponding to different windows on the genome.
        segment_fn_list = glob.glob(all_ct_segment_folder + "/*.bed.gz")
        genome_pos_list = [(x.split('/')[-1]).split('_combined_segment.bed.gz')[0] for x in segment_fn_list] # from /path/to/chr9_14_combined_segment.bed.gz --> chr9_14
        output_fn_list = [os.path.join(predict_outDir, x + "_avg_pred.txt.gz") for x in genome_pos_list]# get the output file names corresponding to different regions on the genome
# 2. partition the list of file names into groups, for later putting into jobs for multiple processes
        num_cores = 1
        partition_segment_fn_list = partition_file_list(segment_fn_list, num_cores)
        partition_output_fn_list = partition_file_list(output_fn_list, num_cores)
        processes = [mp.Process(target = one_job_run_predict_segmentation, args = (partition_segment_fn_list[i], partition_output_fn_list[i], train_cell_types, num_chromHMM_state, train_mode)) for i in range(num_cores)]
        for p in processes:
            p.start()
        for i, p in enumerate(processes):
            p.join()
            print("Process " + str(i) + " is finished!")

def main():
	num_mandatory_args = 6
	if len(sys.argv) != num_mandatory_args:
		usage()
	print("training baseline model")
	train_data_folder = sys.argv[1]
	helper.check_dir_exist(train_data_folder)
	all_ct_segment_folder = sys.argv[2]
	helper.check_dir_exist(all_ct_segment_folder)
	predict_outDir = sys.argv[3]
	helper.make_dir(predict_outDir)
	try:
		num_chromHMM_state = int(sys.argv[4])
		assert num_chromHMM_state > 0, "num_chromHMM_state needs to be positive"
	except:
		print("num_chromHMM_state or num_train_ct is not valid")
		usage()
	# remove the response ct from all_ct list
	all_ct_fn = sys.argv[5]
	print("Done getting command line arguments train_baseline_model.py")
	all_ct = [line.strip() for line in open(all_ct_fn, "r").readlines()]
	# 1. Get the data of predictors and response for training
	train_mode = 'baseline'
	# 2. no need to train model
	# 3. process training data and predict segmentation at each position
	predict_segmentation (all_ct_segment_folder, predict_outDir, all_ct, num_chromHMM_state, train_mode)
	print("Done predicting whole genome")

def usage():
	print("train_baseline_model.py")
	print("train_data_folder: where the state assignment and of training data are stored for all cell types. Each cell type has its own file")
	print("all_ct_segment_folder: where segmentation data of all cell types are stored, for the entire genome, so that we can get data for prediction out.")
	print("predict_outDir: where the beta values obtained from training the data will be stored (e.g. pred_E###")
	print("num_chromHMM_state: Number of chromHMM states that are shared across different cell types")
	print("all_ct_fn: filename of all ct being used in the model to get the baseline representative cell types")
	exit(1)
main()
