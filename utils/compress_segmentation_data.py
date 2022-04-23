import pandas as pd 
import helper 
import numpy as np 
import os 
import argparse
parser = argparse.ArgumentParser(description = 'This code takes in a segmentation file and reaarranges the rows such that if consecutive genomic bins are annotated as the same state, they will be combined. At the same time, the output file will be sorted by chrom, start_bp such that it will be readily available to use for downstream analysis with bedtools or pybedtools')
parser.add_argument('--segment_bed_fn', type = str, required = True,
	help = 'input segment_bed_fn')
parser.add_argument('--output_fn', type = str, required = True, help = 'output_fn')
args = parser.parser_args()
print(args)
helper.check_file_exist(args.segment_bed_fn)
helper.create_folder_for_file(args.output_fn)

def get_compressed_state_segment_one_chrom(chrom_segment_df, chrom):
	# combine multiple rows with the same state into a row so that we can save space on the computer
	# chrom_segment_df has columns 'chrom', 'start_bp', 'end_bp', 'state'
	chrom_segment_df = chrom_segment_df.reset_index(drop = True)
	chrom_segment_df = chrom_segment_df.sort_values(['start_bp']) # sort by ascending start_bp
	result_df = pd.DataFrame(columns = ['chrom', 'start_bp', 'end_bp', 'state'])
	current_start_bin_index = 0
	current_start = chrom_segment_df.loc[current_start_bin_index, 'start_bp']
	current_end = chrom_segment_df.loc[current_start_bin_index, 'end_bp']
	current_state = chrom_segment_df.loc[current_start_bin_index, 'state']
	for index in range(1, chrom_segment_df.shape[0]): # skip the first one because it's already reported into the current data
		if (current_state == chrom_segment_df.loc[index, 'state']) and current_end == chrom_segment_df.loc[index, 'start_bp']: 
			# if in the same state and the previous row and the current row are continous segments on the genome (no gap)
			current_end = chrom_segment_df.loc[index, 'end_bp']
		else: # change of state or there is a gap in the segment --> report into the result_df
			add_row = [chrom, current_start, current_end, current_state]
			result_df.loc[result_df.shape[0]] = add_row
			current_state = chrom_segment_df.loc[index, 'state']
			current_start = chrom_segment_df.loc[index, 'start_bp']
			current_end = chrom_segment_df.loc[index, 'end_bp']
	add_row = [chrom, current_start, current_end, current_state]
	result_df.loc[result_df.shape[0]] = add_row
	# now we are done producing
	return result_df

def read_compressed_segment_fn(segment_bed_fn, output_fn):
	segment_df = pd.read_csv(segment_bed_fn, header = None, sep = '\t', index_col = None)
	try:
		segment_df.columns = ['chrom', 'start_bp', 'end_bp', 'state']
	except:
		print('segment data does not have 4 columns correspoidng to chrom, start_bp, end_bp, state. Please check your input {}'.format(segment_bed_fn))
		exit(1)
	group_df = segment_df.groupby('chrom')
	result_df = pd.DataFrame(columns = ['chrom', 'start_bp', 'end_bp', 'state'])
	for chrom, chrom_segment_df in group_df:
		compressed_chrom_df = get_compressed_state_segment_one_chrom(chrom_segment_df, chrom)
		result_df = result_df.append(compressed_chrom_df)
	result_df = result_df.sort_values(['chrom'])
	result_df.to_csv(output_fn, header = False, index = False, sep = '\t', compression = 'gzip')
	return 

read_compressed_segment_fn(args.segment_bed_fn, args.output_fn)