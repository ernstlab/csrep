#!/usr/bin/env python
'''
If you have a file that is the output of a pipeline that finds 1-1 mappings of genomic bins from one assembly to another, the format of the file would be dest_chrom, dest_start, dest_end, org_bin. org_bin should take the form orgChrom_orgStart_orgEnd. Then, this code will output a reverse file of that, which means the output would have the columns org_chrom, org_start, org_end, destBin. destBin would be the form destChrom_destStart_destEnd
Please use python convert_map_liftOver_to_orgAssembly.py --help for full details on how to run this script.
'''
import pandas as pd 
import numpy as np 
import os
import glob
import helper
import argparse

parser = argparse.ArgumentParser(description = 'If you have a file that is the output of a pipeline that finds 1-1 mappings of genomic bins from one assembly to another, the format of the file would be dest_chrom, dest_start, dest_end, org_bin. org_bin should take the form orgChrom_orgStart_orgEnd. Then, this code will output a reverse file of that, which means the output would have the columns org_chrom, org_start, org_end, destBin. destBin would be the form destChrom_destStart_destEnd')
parser.add_argument('--destOrg_fn', required = True, type=str,
	help = 'the file that shows the mapping of genomic regions from destination assembly to the original assembly')
parser.add_argument('--orgDest_fn', required=True, type=str,
	help = 'the output file that shows the mapping of genomic regions from original assembly to the destination assembly. Format columns: org_chrom, org_start, org_end, destBin. destBin would be the form destChrom_destStart_destEnd')
args = parser.parse_args()
print(args)
helper.check_file_exist(args.destOrg_fn)
helper.create_folder_for_file(args.orgDest_fn)

def convert_map_df_to_orgAssembly(destOrg_fn, orgDest_fn):
	# map_df: fourth will have the form chrom_start_end in the original assembly, the first 3 columns show the coordinate of bins in the destimation assembly
	# --> transform such that the first three columns will show the coordinates in original assembly, while the last one will show the corresponding coordinate in the destination assembly.
	map_df = pd.read_csv(destOrg_fn, header = None, index_col = None, sep = '\t') 
	assert map_df.shape[1] >= 4, 'map_df: {} DOES NOT HAVE 4 COLUMNS. WE ASSUME THAT THE FIRST 4 COLUMNS OF MAP_DF WILL SHOW CHROM, START, END, ORG_BIN, where CHROM, START, END correspond to the coordinates in the destination assembly. Program will exit now'.format(destOrg_fn)
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
	map_df.to_csv(orgDest_fn, header = False, index = False, sep = '\t', compression = 'gzip')
	return

convert_map_df_to_orgAssembly(args.destOrg_fn, args.orgDest_fn)