import pandas as pd 
import numpy as np 
import sys
import os
import glob
import helper

def get_confusion_matrix_data_one_cell_type(fn, num_chromHMM_state):
	df = pd.read_csv(fn, header = 0, index_col = 0, sep = '\t')
	right_columns = list(map(lambda x: 'pred_E' + str(x+1), range(num_chromHMM_state)))
	right_rows = list(map(lambda x: 'true_E' + str(x+1), range(num_chromHMM_state)))
	assert (df.index == right_rows).all(), "The row names are not as expected accoring to num_chromHMM_state"
	assert (df.columns == right_columns).all(), "The column names are not as expected accoring to num_chromHMM_state"
	fraction_true_state_correct_pred = list(map(lambda x: df.iloc[x,x], range(num_chromHMM_state))) # each entry corresponds to a state, the reported value for the state is the fraction of the state's true segmentations that are correctly predicted to that state in the representative data that excludes this state.
	return (fraction_true_state_correct_pred)

def calculate_summary_staistics_across_ct(cg_dir, out_dir, num_chromHMM_state, ct_list):
	ct_confusion_matrix_fn_list = list(map(lambda x: os.path.join(cg_dir, 'val_' + x, 'validation_results', 'confusion_matrix_prediction_performance.txt.gz'), ct_list))
	# get files like /Users/vuthaiha/Desktop/window_hoff/diff_pete/roadmap/blood/baseline/val_E034/roc/tpr_fpr_all_states.txt.gz for all cell types
	map(helper.check_file_exist, ct_confusion_matrix_fn_list) # make sure that we have all the roc data files that we need
	column_names = list(map(lambda x: 'fraction_true_S' + str(x + 1) + "_correct", range(num_chromHMM_state)))
	sum_df = pd.DataFrame(columns = column_names)
	for ct_index, ct in enumerate(ct_list):
		this_ct_diag_list = get_confusion_matrix_data_one_cell_type(ct_confusion_matrix_fn_list[ct_index], num_chromHMM_state) # a list of values for the states for this cell type
		sum_df.loc[ct_index] = this_ct_diag_list
	sum_df.index = ct_list # index of the result dataframe
	output_fn = os.path.join(out_dir, 'summary_confusion_matrix.txt.gz')
	sum_df.to_csv(output_fn, header = True, index = True, sep = '\t', compression = 'gzip') # save to a pickle file


def main():
	if len(sys.argv) != 5:
		usage()
	cg_dir = sys.argv[1]
	helper.check_dir_exist(cg_dir)
	out_dir = sys.argv[2]
	helper.make_dir(out_dir)
	num_chromHMM_state = helper.get_command_line_integer(sys.argv[3])
	cell_type_list_fn = sys.argv[4]
	ct_list = helper.get_list_from_line_seperated_file(cell_type_list_fn)
	helper.check_file_exist(cell_type_list_fn)
	print ("Done getting command line arguments")
	calculate_summary_staistics_across_ct(cg_dir, out_dir, num_chromHMM_state, ct_list)
	print ("Done!")

def usage():
	print ("Given that results of validation the representative chormatin state map for a group of cell types, we want to summarize the results of all the confusion matrices that we got from analyses that leave one cell type out at a time. Each of such analysis will calculate how well the representative data from the non-excluded cell types can help predict the segmentation of the excluded cell types")
	print ("python get_summary_cofusion_matrix_data.py")
	print ("cg_dir: directory where the validation data for all the cell types are stored for this cell group")
	print ("out_dir: where the summazied roc data for across the cell types inside this cell group is stored")
	print ("num_chromHMM_state")
	print ("cell_type_list_fn: fn where each cell types are stored separate lines")
	exit(1)
main()