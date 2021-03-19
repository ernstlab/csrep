# tested 03/08/2021
import os
import sys
import helper
import pandas as pd
import glob
import numpy as np
import multiprocessing as mp

def process_one_ct_one_chrom_segment_fn(fn, ct, chrom):
	df = pd.read_csv(fn, header = None, index_col = None, sep = '\t')
	df.columns = ['chrom', 'start_bp', 'end_bp', 'state']
	df = df[df['chrom'] == chrom] # filter only rows with the chromosome that we are interested in
	df['num_bin'] =  ((df['end_bp'] - df['start_bp']) / helper.NUM_BP_PER_BIN).astype('int')
	segment_series = []
	for index, row in df.iterrows():
		segment_series.extend([row['state']] * row['num_bin'])
	return segment_series

def print_one_chrom_data_to_text(chrom_df, chrom, output_folder):
	# save the data from chromsome into multiple files, each corresponding to a region on the chromosome
	num_files_to_print = int(np.ceil(chrom_df.shape[0] / helper.NUM_BIN_PER_WINDOW))
	for file_index in range(num_files_to_print - 1):
		start_row_index = file_index * helper.NUM_BIN_PER_WINDOW
		end_row_index = start_row_index + helper.NUM_BIN_PER_WINDOW
		save_fn = os.path.join(output_folder, chrom + '_' + str(file_index) + '_combined_segment.bed.gz')
		save_df = chrom_df.loc[start_row_index:end_row_index]
		save_df.to_csv(save_fn, header = True, index = False, sep = '\t', compression = 'gzip')
	last_index = num_files_to_print - 1
	last_save_df = chrom_df.loc[last_index * helper.NUM_BIN_PER_WINDOW:]
	last_save_fn = os.path.join(output_folder, chrom + '_' + str(last_index) + '_combined_segment.bed.gz')
	last_save_df.to_csv(last_save_fn, header = True, index = False, sep = '\t', compression = 'gzip')
	return 

def process_one_chrom_all_ct(input_fn_list, ct_list, chrom, output_folder):
	result_df = pd.DataFrame()
	for ct_index, ct in enumerate(ct_list):
		segment_series = process_one_ct_one_chrom_segment_fn(input_fn_list[ct_index], ct, chrom)
		if len(segment_series) == 0: 
			print ('chromsome data: ' + chrom + ' and cell type: ' + ct + ' show no segmentation data. You may want to check your input!')
			return
		if (ct_index != 0):
			assert len(segment_series) == result_df.shape[0], 'The number of bins in cell type: ' + ct + ' does not match with other cell types in chromosome ' + chrom + ' (' + str(len(segment_series)) + '). Check your input segmentation data, all cell types should have the same chromosome length in the segmentation data.'
		result_df[ct] = segment_series
	result_df['chrom'] = chrom
	result_df['start_bp_this_window'] = helper.NUM_BP_PER_BIN * result_df.index
	result_df['end_bp_this_window'] = result_df['start_bp_this_window'] + helper.NUM_BP_PER_BIN
	result_df = result_df[['chrom', 'start_bp_this_window', 'end_bp_this_window'] + ct_list]
	print_one_chrom_data_to_text(result_df, chrom, output_folder)
	return 

def combine_all_ct_segment_folder_one_process(input_fn_list, ct_list, chrom_list, output_folder):
	for chrom in chrom_list:
		chrom = 'chr' + chrom # currenntly chrom is just 1 --> 22, X. Again, we exclude Y because not every samples have Ys data
		process_one_chrom_all_ct(input_fn_list, ct_list, chrom, output_folder)
	return 

def combine_all_ct_segment_folder(input_folder, output_folder, input_suffix):
	input_fn_list = glob.glob(input_folder + '/*' + input_suffix) # ex: input_folder/E003_chr22_core_K27ac_segments.bed.gz
	print(input_fn_list)
	ct_list = list(map(lambda x: x.split('/')[-1].split(input_suffix)[0], input_fn_list)) # E003
	print (ct_list)
	num_cores = 4
	partition_chrom_list = helper.partition_file_list(helper.CHROMOSOME_LIST, num_cores)
	processes = [mp.Process(target = combine_all_ct_segment_folder_one_process, args = (input_fn_list, ct_list, partition_chrom_list[i], output_folder)) for i in range(num_cores)]
	for p in processes:
		p.start()
	for i, p in enumerate(processes):
		p.join()
		print('Done with process ' + str(i + 1))
	return 

def main():
	if len(sys.argv) != 4:
		usage()
	input_folder = sys.argv[1]
	helper.check_dir_exist(input_folder)
	output_folder = sys.argv[2]
	helper.make_dir(output_folder)
	input_suffix = sys.argv[3]
	print ("Done getting command line arguments")
	combine_all_ct_segment_folder(input_folder, output_folder, input_suffix)
	print ('Done!')
def usage():
	print ("python get_all_ct_segment_folder.py")
	print ("input_folder: where all the raw data for all the ct are stored")
	print ("output_folder: where we store files that correspond to regions on the genome, and within each file, we store the data of segmentation for all the ct that are in the input_folder")
	print ("input_suffix: inside the input_folder, the file names are <ct><input_suffix>. In testdata, it would be chr22_core_K27ac_segments.bed.gz")
	exit(1)
main()
