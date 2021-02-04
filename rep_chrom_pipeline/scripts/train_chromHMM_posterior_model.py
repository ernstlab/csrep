import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
import os
import sys
import time
import glob
sys.path.append("/u/project/ernst/havu73/source_pete/train_and_evaluate/")
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


def get_predictorX_chromHMM_pos_data(all_ct_chrom_pos_folder, chromosome, train_cell_types, num_chromHMM_state):
    # given the all_ct_chrom_pos_folder where all ct's posterior data are store, we want to out the files corresponding to the train_cell_types, and put them into a list so that we can average them out, which gives us the prediction for the chromatin state assignments for states across the chromsome. 
    # Return: a list of df, each corresponds to a train_cell_type. inside each df: rows: states, columns: regions on the chromosome that we are reading data in. 
    cPos_df_list = []
    one_train_ct_colnames = list(map(lambda x: 'state_' + str(x+1), range(num_chromHMM_state)))
    for ct_index, ct in enumerate(train_cell_types):
        this_ct_cPos_fn = os.path.join(all_ct_chrom_pos_folder, ct, ct + '_18_core_K27ac_' + chromosome + '_posterior.txt.gz') # example: /u/home/h/havu73/project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded/posterior/E003/E003_18_core_K27ac_chr9_posterior.txt.gz
        this_ct_cPos_df = pd.read_csv(this_ct_cPos_fn, header = 0, skiprows = 1, index_col = None, sep = '\t')
        this_ct_cPos_df.columns = one_train_ct_colnames
        cPos_df_list.append(this_ct_cPos_df)
    return cPos_df_list

def save_predictor_df_into_multiple_files(response_df, chromosome, predict_outDir, num_chromHMM_state):
    total_segment = response_df.shape[0] # number of rows is the number of semgents here.
    num_output_fn = int(np.ceil(total_segment / helper.NUM_BIN_PER_WINDOW))
    print(num_output_fn)
    for window_index in range(num_output_fn - 1): #except the last window, which we will save after this
        output_fn = os.path.join(predict_outDir, chromosome + '_' + str(window_index) + '_avg_pred.txt.gz')
        this_window_start_row_index = int(window_index * helper.NUM_BIN_PER_WINDOW)
        this_windw_end_row_index = int((window_index + 1) * helper.NUM_BIN_PER_WINDOW) 
        this_window_response_df = response_df.iloc[range(this_window_start_row_index, this_windw_end_row_index), :]
        this_window_response_df = this_window_response_df.reset_index(drop = True) #  now the index is always from 0 to the end of the window
        this_window_response_df.to_csv(output_fn, header = True, index = False, sep = '\t', compression = 'gzip')
        print("Done saving file: " + output_fn)
    last_output_fn = os.path.join(predict_outDir, chromosome + '_' + str(num_output_fn - 1) + '_avg_pred.txt.gz')
    last_window_start_row_index = int((num_output_fn - 1) * helper.NUM_BIN_PER_WINDOW)
    last_window_end_row_index = int(response_df.shape[0]) # the last row
    last_window_response_df = response_df.iloc[range(last_window_start_row_index, last_window_end_row_index), :]
    last_window_response_df.to_csv(last_output_fn, header = True, index = False, sep = '\t', compression = 'gzip')
    print ('Done saving file: ' + last_output_fn)
    return 

def predict_segmentation_one_genomic_window(all_ct_chrom_pos_folder, chromosome, predict_outDir, train_cell_types, num_chromHMM_state):
    # 1. Get the data of predictor cell types into the right format
    predictor_df_list = get_predictorX_chromHMM_pos_data(all_ct_chrom_pos_folder, chromosome, train_cell_types, num_chromHMM_state) # a list of df. Within each df: rows: states, columns: positions inside a  chromosome on the genome
    # 2. Do the prediction job. For this model, it is simply a averaging over all the train_ct
    response_df = pd.concat(predictor_df_list).groupby(level = 0).mean() # average over all the available dataframes, each dataframe corresponds to a train_ct
    # now save the files such that each file only corresponds to each region on that chromosome. 
    save_predictor_df_into_multiple_files(response_df, chromosome, predict_outDir, num_chromHMM_state)
    return

def one_job_run_predict_segmentation(all_ct_chrom_pos_folder, genome_pos_list, predict_outDir, train_cell_types, num_chromHMM_state):
    '''
    segment_fn_list and output_fn_list: the orders of regions in these two lists are similar (look at function predict_segmentation).
    Each element corresponds to a region on the genome
    '''
    for (window_index, chromosome) in enumerate(genome_pos_list):
        predict_segmentation_one_genomic_window(all_ct_chrom_pos_folder, chromosome, predict_outDir, train_cell_types, num_chromHMM_state)
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

def predict_segmentation (all_ct_chrom_pos_folder, predict_outDir, train_cell_types, num_chromHMM_state):
        # 1. Get list of segmentation files corresponding to different windows on the genome.
    segment_fn_list = glob.glob(all_ct_chrom_pos_folder + '/' +  train_cell_types[0] + "/*.txt.gz")
    genome_pos_list = [((x.split('/')[-1]).split('_posterior.txt.gz')[0]).split('_18_core_K27ac_')[1] for x in segment_fn_list] # from /u/home/h/havu73/project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded/posterior/E003/E003_18_core_K27ac_chr9_posterior.txt.gz --> chr9
    # 2. partition the list of file names into groups, for later putting into jobs for multiple processes
    genome_pos_list.remove('chrY') # we focus our analysis on only the normal chromosome and chromosome X
    num_cores = 2
    partition_genome_pos_list = partition_file_list(genome_pos_list, num_cores)#genome_pos_list[:3], num_cores) # the input segment fn are actually by chromosome, because this file is all about chromHMM posterior data
    processes = [mp.Process(target = one_job_run_predict_segmentation, args = (all_ct_chrom_pos_folder,  partition_genome_pos_list[i], predict_outDir, train_cell_types, num_chromHMM_state)) for i in range(num_cores)]
    for p in processes:
        p.start()
    for i, p in enumerate(processes):
        p.join()
        print("Process " + str(i) + " is finished!")

def main():
	num_mandatory_args = 5
	if len(sys.argv) != num_mandatory_args:
		usage()
	print("training chromHMM_pos model")
	all_ct_chrom_pos_folder = sys.argv[1]
	helper.check_dir_exist(all_ct_chrom_pos_folder)
	predict_outDir = sys.argv[2]
	helper.make_dir(predict_outDir)
	try:
		num_chromHMM_state = int(sys.argv[3])
		assert num_chromHMM_state > 0, "num_chromHMM_state needs to be positive"
	except:
		print("num_chromHMM_state or num_train_ct is not valid")
		usage()
	# remove the response ct from all_ct list
	all_ct_fn = sys.argv[4]
	print("Done getting command line arguments train_chromHMM_pos_model.py")
	all_ct = [line.strip() for line in open(all_ct_fn, "r").readlines()]
	train_cell_types = all_ct # a mistake in the code, but let's keep it like this for now
	# 2. no need to train model
	# 3. process training data and predict segmentation at each position 
	predict_segmentation (all_ct_chrom_pos_folder, predict_outDir, train_cell_types, num_chromHMM_state)
	print("Done predicting whole genome")

def usage():
	print("train_chromHMM_posterior_model.py")
	print("all_ct_chrom_pos_folder: where segmentation data of all cell types are stored, for the entire genome, so that we can get data for prediction out.")
	print("predict_outDir: where the beta values obtained from training the data will be stored (e.g. val_E###/average_predictions/")
	print("num_chromHMM_state: Number of chromHMM states that are shared across different cell types")
	print("all_ct_fn: filename of all ct being used in the model to produce the representative data of all these samples/ct")
	exit(1)
main()
