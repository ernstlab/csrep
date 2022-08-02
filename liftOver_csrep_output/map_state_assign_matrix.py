#!/usr/bin/env python
'''
Given the summary chromatin state assignment probabilities outputted by CSREP/ base_count in one assembly, this approach will map the summary chromatin state maps to another assembly. Please use python map_state_assign_matrix.py --help to see full details on how to run this script.
'''
import pandas as pd 
import numpy as np 
import os
import glob
import helper
import argparse
import multiprocessing as mp

parser = argparse.ArgumentParser(description = 'Given a summary chromatin state map and a 1-1 mapping of genomic bins from one assembly to another (with all overlapping mapped regions removed), we will create a mapped summary chromatin state map in the destimation assembly')
parser.add_argument('--csrep_folder', type=str, required=True,
	help = 'where there are files, each representing a region in the genome. The data should be the summary chromatin state assigment matrix')
parser.add_argument('--orgDest_fn', type=str, required=True,
	help = 'the 1-1 mapping of genomic bins from the orgiginal assembly to destimation assembly. This should be the output of convert_map_liftOver_to_orgAssembly.py. Format columns: org_chrom, org_start, org_end, destBin. destBin would be the form destChrom_destStart_destEnd')
parser.add_argument('--output_folder', type=str, required=True,
	help = 'where there are files, each representing a region in the genome. The data should be the summary chromatin state assigment matrix IN THE DESTINATION ASSEMBLY')
parser.add_argument('--chromhmm_state_num', type=int, required=True,
	help = 'chromhmm_state_num')
parser.add_argument('--rewrite_existing_chrom', action='store_true',
	help = 'if this flag is present, we will rewrite results for all chromosomes, even if the output files associated with some of them are already produced')
parser.set_defaults(rewrite_existing_chrom=False)
args = parser.parse_args()
print(args)
helper.check_dir_exist(args.csrep_folder)
helper.check_file_exist(args.orgDest_fn)
helper.make_dir(args.output_folder)

def read_rep_df(fn, chromhmm_state_num):
	right_colnames = list(map(lambda x: 'state_' + str(x+1), range(chromhmm_state_num)))
	print(fn)
	try:
		df = pd.read_csv(fn, header = 0, index_col = 0, sep = '\t')
		assert list(df.columns) == right_colnames, 'index_col = 0 is wrong'
	except:
		df = pd.read_csv(fn, header = 0, index_col = None, sep = '\t')
		assert list(df.columns) == right_colnames, 'index_col = 0 is wrong'
	return df


def map_assign_matrix_one_orgFn(org_fn, chrom_map_df, region_index, chromhmm_state_num):
	'''
	chrom_map_df: org_chrom, org_start, org_end, destBin. destBin would be the form destChrom_destStart_destEnd. The data here should be restricted to one genomic region on the original assembly
	'''
	org_df = read_rep_df(org_fn, chromhmm_state_num)
	region_start_bp = helper.NUM_BP_PER_WINDOW * region_index
	region_end_bp = region_start_bp + helper.NUM_BP_PER_WINDOW
	region_map_df = chrom_map_df[(chrom_map_df['start'] >= region_start_bp) & (chrom_map_df['end'] < region_end_bp)]
	region_map_df['index_within_region'] = ((region_map_df['start'] - region_start_bp) / helper.NUM_BP_PER_BIN).astype(int)
	region_map_df = region_map_df.merge(org_df, left_on = 'index_within_region', right_index = True, how = 'inner')
	region_map_df.drop(['chrom', 'start', 'end'], axis = 1, inplace=True)
	region_map_df['chrom'] = region_map_df['destBin'].apply(lambda x: x.split('_')[0])
	region_map_df['start'] = region_map_df['destBin'].apply(lambda x: int(x.split('_')[1]))
	region_map_df['end'] = region_map_df['destBin'].apply(lambda x: int(x.split('_')[2]))
	state_columns = region_map_df.columns[pd.Series(region_map_df.columns).str.startswith('state')]
	region_map_df = region_map_df[['chrom', 'start', 'end'] + list(state_columns)] # the chrom, start, end show the coordinates of summary state assignment in the destination assembly
	return region_map_df

def write_result_destChrom_diff_orgChrom(other_chrom_df, orgChrom, output_folder):
	# this is the case where regions from orgChrom in the original assembly got mapped to a different chromosome in the destination assembly. In this case, we will first filter out all regions that are not in the orgChom in this lifted-over summary chromatin state map. Then, we will write them in separate files so that later we can merge back to one file per chromosome
	group_other_chrom_df = other_chrom_df.groupby('chrom')
	for chrom, oChrom_df in group_other_chrom_df:
		save_fn = os.path.join(output_folder, '{o}_to_{d}_prob_state_map.txt.gz'.format(o=orgChrom, d=chrom))
		oChrom_df.to_csv(save_fn, header = True, index = False, sep = '\t', compression = 'gzip')
	return 

def map_assign_matrix_one_chrom(csrep_folder, chrom, orgDest_df, chromhmm_state_num, output_folder):
	chrom = 'chr{}'.format(chrom)
	orgDest_df = orgDest_df[orgDest_df['chrom'] == chrom] # i named this orgDest_df instead of something else because I wanted to save computation space, orgDest_df for the whole genome is itself a very large dataframe
	num_org_fn = len(glob.glob(csrep_folder + '/{}_*_avg_pred.txt.gz'.format(chrom)))
	output_colnames = ['chrom', 'start', 'end'] + list(map(lambda x: 'state_{}'.format(x+1), range(chromhmm_state_num)))
	result_df = pd.DataFrame(columns = output_colnames)
	for region_index in range(num_org_fn):
		org_fn = os.path.join(csrep_folder, '{c}_{i}_avg_pred.txt.gz'.format(c=chrom, i=region_index))
		this_region_df = map_assign_matrix_one_orgFn(org_fn, orgDest_df, region_index, chromhmm_state_num)
		result_df = result_df.append(this_region_df) # columns: chrom, start, end, state_1 --> state_18
	# now that we got all the liftedOver data for regions in the chrom in the original assembly
	# we will first filter out only region that are in this chrom first
	other_chrom_df = result_df[result_df['chrom'] != chrom] 
	if other_chrom_df.shape[0] >0: # if there are regions originally from this chrom but then get mapped to other chromosomes, we will write that data into separate files
		write_result_destChrom_diff_orgChrom(other_chrom_df, chrom, output_folder)
	# after we are done writing to files the regions that got mapped to other chromosomes, now just save the mapping of regions that got mapped to the same chromosome. 
	result_df = result_df[result_df['chrom'] == chrom]
	result_df = result_df.sort_values('start') # sort by ascending start bp
	save_fn = os.path.join(output_folder, '{}_liftOver_probState.txt.gz'.format(chrom))
	result_df.to_csv(save_fn, header=True, index=False, sep='\t', compression='gzip')
	return 

def map_assign_matrix_multiple_chrom(csrep_folder, chrom_list, orgDest_df, chromhmm_state_num, output_folder):
	# for each chrom C, map all regions in C from the original chrom to the destination assembly 
	# this function should be called for each of the processes running in parallel during multiprocessing
	for chrom in chrom_list:
		map_assign_matrix_one_chrom(csrep_folder, chrom, orgDest_df, chromhmm_state_num, output_folder)
	return

def combine_mapped_data_from_multiChrom_to_oneChrom_in_dest_assembly(output_folder, chrom):
	# chrom is just '1', '2', etc. not 'chr1', etc.
	chrom = 'chr{}'.format(chrom)
	fn_list = glob.glob(output_folder + '*_to_{}_prob_state_map.txt.gz'.format(chrom))
	if len(fn_list) ==  0:
		return # do nothing, because there are not regions from other chromosome in org_assembly that got mapped to this chromosome in the dest_assembly
	fn_list += ['{}_liftOver_probState.txt.gz'.format(chrom)]
	df_list = list(map(lambda x: pd.read_csv(x, header=0, index_col = None, sep = '\t'), fn_list))
	result_df = pd.concat(df_list)
	result_df = result_df.sort_values('start')
	save_fn = os.path.join(output_folder, '{}_liftOver_probState.txt.gz'.format(chrom))
	result_df.to_csv(save_fn, header = True, index = False, sep = '\t', compression = 'gzip')
	return 

def combine_mapped_data_one_process(output_folder, chrom_list):
	for chrom in chrom_list:
		combine_mapped_data_from_multiChrom_to_oneChrom_in_dest_assembly(output_folder, chrom)
	return 

def get_chrom_list(output_folder, rewrite_existing_chrom):
	if rewrite_existing_chrom:
		return helper.CHROMOSOME_LIST
	chrom_segment_fn_list = glob.glob(output_folder + '/chr*_liftOver_probState.txt.gz')
	chrom_list = [(x.split('/')[-1]).split('_gen_pos_segment_fn_list')[0].split('chr')[1] for x in chrom_segment_fn_list] # get the list of all genomic positions available: [chr9_11, chr9_12, etc.]
	chrom_list = np.setdiff1d(helper.CHROMOSOME_LIST, chrom_list)
	print("Number of positions that will be calculated in average_pred_results: " + str(len(chrom_list)))
	print(chrom_list)
	return chrom_list

def map_assign_matrix_allChrom_multiProcess(csrep_folder, orgDest_fn, chromhmm_state_num, output_folder, rewrite_existing_chrom):
	orgDest_df = pd.read_csv(orgDest_fn, header = None, index_col=None, sep='\t')
	orgDest_df.columns = ['chrom', 'start', 'end', 'destBin']
	num_cores = 4
	chrom_list = get_chrom_list(output_folder, rewrite_existing_chrom)
	partition_chrom_list = helper.partition_file_list(chrom_list, num_cores)
	processes = [mp.Process(target = map_assign_matrix_multiple_chrom, args =(csrep_folder, partition_chrom_list[i], orgDest_df, chromhmm_state_num, output_folder)) for i in range(num_cores)]
	for p in processes:
		p.start()
	for i, p in enumerate(processes):
		p.join()
		print("Process " + str(i) + " to map individual chromosomes is finished!")
	print('Done first round of mapping regions between assemblies')
	# now  that we wrote down all the data for individual chromosomes. For each chromosome, if this chrom was mapped from multiple chromosomes in the original assembly, we will need to combine that data into one file for the destination chromosome
	processes = [mp.Process(target = combine_mapped_data_one_process, args = (output_folder, partition_chrom_list[i])) for i in range(num_cores)]
	for p in processes:
		p.start()
	for i, p in enumerate(processes):
		p.join()
		print("Process " + str(i) + " to combine regions that were mapped to the same chromosome in dest_assembly is finished!")
	print ('Done!')
	return

map_assign_matrix_allChrom_multiProcess(args.csrep_folder, args.orgDest_fn, args.chromhmm_state_num, args.output_folder, args.rewrite_existing_chrom)
