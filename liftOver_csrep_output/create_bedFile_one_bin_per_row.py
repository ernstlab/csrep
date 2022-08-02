#!/usr/bin/env python

'''
Given the file that shows the length of chromosomes, create a file with 4 columns chrom, start, end, chrom_start_end. The purpose of the output is to eventually find a 1-1 mapping of regions on hg19 to region on hg38, which we will eventually use to map data of summary chromatin state maps from hg19 to hg38. In other words, output from this file will be passed to liftOver software.
Please use python create_bedFile_one_bin_per_row.py --help for full details on how to run this script.
'''
import pandas as pd 
import numpy as np 
import os
import glob
import helper
import argparse
parser = argparse.ArgumentParser(description = 'Given the file that shows the length of chromosomes, create a file with 4 columns chrom, start, end, chrom_start_end. The purpose of the output is to eventually find a 1-1 mapping of regions on hg19 to region on hg38, which we will eventually use to map data of summary chromatin state maps from hg19 to hg38. In other words, output from this file will be passed to liftOver software.')
parser.add_argument('--chrom_length_fn', type = str, required=True,
	help = '3 columns: chrom, start(0), end (chrom_lengths)')
parser.add_argument('--output_fn', type=str, required=True,
	help = 'should be .gz')
parser.add_argument('--NUM_BP_PER_BIN', type=int, required=False, default=200,
	help = 'size of the bin that we are trying to divide the genome into')
args = parser.parse_args()
print(args)
helper.check_file_exist(args.chrom_length_fn)
helper.create_folder_for_file(args.output_fn)

def create_result_df_one_chrom(chrom_length, chrom, NUM_BP_PER_BIN):
	df = pd.DataFrame() # name should have the format: chrom_start_end
	num_bins = int(chrom_length / NUM_BP_PER_BIN) # round down
	df['start'] = np.arange(num_bins) * NUM_BP_PER_BIN
	df['end'] = df['start'] + NUM_BP_PER_BIN
	df['chrom'] = chrom
	df['name'] = df.apply(lambda x: '{c}_{s}_{e}'.format(c=x['chrom'], s=x['start'], e=x['end']), axis = 1)
	df = df[['chrom', 'start', 'end', 'name']]
	return df

def create_bedFile_one_bin_per_row(chrom_length_fn, output_fn, NUM_BP_PER_BIN):
	chrom_len_df = pd.read_csv(chrom_length_fn, header = None, sep = '\t', index_col = None)
	chrom_len_df.columns = ['chrom', 'start', 'end']
	result_df = pd.DataFrame(columns = ['chrom', 'start', 'end', 'name'])
	for index, row in chrom_len_df.iterrows():
		this_chrom_df = create_result_df_one_chrom(row['end'], row['chrom'], NUM_BP_PER_BIN)
		result_df = result_df.append(this_chrom_df)
	result_df = result_df.sort_values(['chrom', 'start'])
	result_df.to_csv(output_fn, header = False, index = False, sep = '\t', compression = 'gzip')
	return 

create_bedFile_one_bin_per_row(args.chrom_length_fn, args.output_fn, args.NUM_BP_PER_BIN)