import pandas as pd 
import numpy as np 
import os
import glob
import helper
import argparse
import pybedtools as bed 
parser = argparse.ArgumentParser(description = 'Given a summary chromatin state map and a 1-1 mapping of genomic bins from one assembly to another (with all overlapping mapped regions removed), we will create a mapped summary chromatin state map in the destimation assembly')
parser.add_argument('--sum_fn', type=str, required=True,
	help = 'the summary chromatin state map in original assembly')
parser.add_argument('--map_fn', type=str, required=True,
	help = 'the 1-1 mapping of genomic bins from one region to another. This should be the output of the liftOver pipeline that will also get rid of regions that got mapped from multiple regions in the original assembly')
parser.add_argument('--skip_rows_sum', type=int, required=False, default=1,
	help = 'the first N rows that we will skip when we read sum_fn. This is because if we write the summary data in ucsc genome browser format, the first row should not be read into pandas dataframe because it\'s about the browser setting. But if we had the sum_fn as a normal bed file then this number should be set to 0.')
parser.add_argument('--output_fn', type=str, required=True,
	help = 'output_fn. should be .gz')
args = parser.parse_args()
print(args)
helper.check_file_exist(args.sum_fn)
helper.check_file_exist(args.map_fn)
helper.create_folder_for_file(args.output_fn)

def get_compressed_state_segment_one_chrom(chrom_segment_df, chrom):
	# combine multiple rows with the same state into a row so that we can save space on the computer
	# chrom_segment_df has columns 'chrom', 'start', 'end', 'state'
	chrom_segment_df = chrom_segment_df.reset_index(drop = True)
	chrom_segment_df = chrom_segment_df.sort_values(['start']) # sort by ascending start
	result_df = pd.DataFrame(columns = ['chrom', 'start', 'end', 'state'])
	current_start_bin_index = 0
	current_start = chrom_segment_df.loc[current_start_bin_index, 'start']
	current_end = chrom_segment_df.loc[current_start_bin_index, 'end']
	current_state = chrom_segment_df.loc[current_start_bin_index, 'state']
	for index in range(1, chrom_segment_df.shape[0]): # skip the first one because it's already reported into the current data
		if (current_state == chrom_segment_df.loc[index, 'state']) and current_end == chrom_segment_df.loc[index, 'start']: 
			# if in the same state and the previous row and the current row are continous segments on the genome (no gap)
			current_end = chrom_segment_df.loc[index, 'end']
		else: # change of state or there is a gap in the segment --> report into the result_df
			add_row = [chrom, current_start, current_end, current_state]
			result_df.loc[result_df.shape[0]] = add_row
			current_state = chrom_segment_df.loc[index, 'state']
			current_start = chrom_segment_df.loc[index, 'start']
			current_end = chrom_segment_df.loc[index, 'end']
	add_row = [chrom, current_start, current_end, current_state]
	result_df.loc[result_df.shape[0]] = add_row
	# now we are done producing
	return result_df

def compress_segmentation(segment_df, output_fn):
	# segment_df should have 4 columns: chrom, start, end, state
	group_df = segment_df.groupby('chrom')
	result_df = pd.DataFrame(columns = ['chrom', 'start', 'end', 'state'])
	for chrom, chrom_segment_df in group_df:
		compressed_chrom_df = get_compressed_state_segment_one_chrom(chrom_segment_df, chrom)
		result_df = result_df.append(compressed_chrom_df)
	result_df = result_df.sort_values(['chrom'])
	result_df.to_csv(output_fn, header = False, index = False, sep = '\t', compression = 'gzip')
	return 

def read_summary_chrom_state_fn(sum_fn, skip_rows_sum):
	sum_df = pd.read_csv(sum_fn, skiprows = skip_rows_sum, header = None, index_col = None, sep = '\t')
	assert sum_df.shape[1] >= 4, 'sum_df: {} DOES NOT HAVE 4 COLUMNS. WE ASSUME THAT THE FIRST 4 COLUMNS OF SUM_DF WILL SHOW CHROM, START, END, STATE. Program will exit now'.format(sum_fn)
	sum_df = sum_df[range(4)] # only need the first 4 columns
	sum_df.columns = ['chrom', 'start', 'end', 'state']
	print('NOTE: we assume that sum_fn contains genome coordinate data that has been sorted by the code used to generate this file, so we do not sort it here. If your input data is generated otherwise and it is not sorted, then pybedtools may crash. Here, we do not sort the data to save computation time in this code.')
	sum_bed = bed.BedTool.from_dataframe(sum_df)
	return sum_bed


def convert_map_df_to_orgAssembly(map_fn):
	# map_df: fourth will have the form chrom_start_end in the original assembly, the first 3 columns show the coordinate of bins in the destimation assembly
	# --> transform such that the first three columns will show the coordinates in original assembly, while the last one will show the corresponding coordinate in the destination assembly.
	map_df = pd.read_csv(map_fn, header = None, index_col = None, sep = '\t') 
	assert map_df.shape[1] >= 4, 'sum_df: {} DOES NOT HAVE 4 COLUMNS. WE ASSUME THAT THE FIRST 4 COLUMNS OF MAP_DF WILL SHOW CHROM, START, END, ORG_BIN, where CHROM, START, END correspond to the coordinates in the destination assembly. Program will exit now'.format(map_fn)
	map_df = map_df[range(4)]
	map_df.columns = ['chrom', 'start', 'end', 'org_bin'] # org_bin will have the form chrom_start_end in the original assembly, the first 3 columns show the coordinate of bins in the destimation assembly
	map_df['dest_bin'] = map_df.apply(lambda x: '{c}_{s}_{e}'.format(c=x['chrom'], s=x['start'], e=x['end']), axis = 1)
	map_df.drop(['chrom', 'start', 'end'], axis = 1, inplace = True) # drop columns
	map_df['chrom'] = map_df.apply(lambda x: x['org_bin'].split('_')[0], axis = 1) # get the chrom in the original assembly
	map_df['start'] = map_df.apply(lambda x: int(x['org_bin'].split('_')[1]), axis = 1) # get the start in the original assembly
	map_df['end'] = map_df.apply(lambda x: int(x['org_bin'].split('_')[2]), axis = 1) # get the start in the original assembly
	map_df.drop('org_bin', axis=1, inplace=True) 
	map_df = map_df.sort_values(['chrom', 'start'])
	map_df = map_df[['chrom', 'start', 'end', 'dest_bin']]
	map_bed = bed.BedTool.from_dataframe(map_df)
	return map_bed

def map_summary_chromatin_state_map(sum_fn, map_fn, output_fn, skip_rows_sum):
	sum_bed = read_summary_chrom_state_fn(sum_fn, skip_rows_sum)
	map_bed = convert_map_df_to_orgAssembly(map_fn) # chrom, start, end, dest_bin: first 3 show the coordinates in orignal assembly, last column show the corresponding assembly in destination assembly
	print('Done reading in sum_df and map_df')
	map_bed = map_bed.map(sum_bed, c=4, o='collapse') # chrom, start, end, dest_bin, state
	map_bed = map_bed.to_dataframe() # still give it the variable name map_bed even though it is a dataframe now because that would save a significant amount of memory. map_bed has #rows ~ 15 millions
	map_bed.columns = ['org_chrom', 'org_start', 'org_end', 'dest_bin', 'state']
	map_bed.drop(['org_chrom', 'org_start', 'org_end'], axis = 1, inplace = True)
	map_bed['chrom'] = map_bed.apply(lambda x: x['dest_bin'].split('_')[0], axis = 1) # get the chrom in the destination assembly
	map_bed['start'] = map_bed.apply(lambda x: x['dest_bin'].split('_')[1], axis = 1) # get the chrom in the destination assembly
	map_bed['end'] = map_bed.apply(lambda x: x['dest_bin'].split('_')[2], axis = 1) # get the chrom in the destination assembly
	map_bed.drop('dest_bin', axis=1, inplace=True)
	map_bed = map_bed[['chrom', 'start', 'end', 'state']]
	compress_segmentation(map_bed, output_fn) # compress the results such that consecutive segments of the same state will be combined into one row
	print ('Done getting the transformed data into the destination assembly')
	return 

map_summary_chromatin_state_map(args.sum_fn, args.map_fn, args.output_fn, args.skip_rows_sum)