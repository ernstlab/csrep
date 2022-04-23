import pandas as pd 
import helper 
import numpy as np 
import os 
import pybedtools as bed 
import sys

def read_state_annot(state_annot_fn):
	state_annot_df = pd.read_csv(state_annot_fn, header = 0, index_col = False, sep = '\t')
	assert state_annot_df.columns[0] == 'state' and state_annot_df.columns[1] == 'mnenomic', 'The first two columns names of the state annotation file should be: state, mnenomic. Please refer to the sample state annotation file provided with this package for an example'
	state_annot_dict = dict(zip(state_annot_df.mnenomic, state_annot_df.state))
	print(state_annot_dict)
	return state_annot_dict # keys: mnenomic, values: state index (1 based)

def transform_state_into_desired_form(state, state_annot_dict):
	try: 
		int_state = int(state[1:]) # if it is in the form such as E1, E2, ... E18
	except:
		try:
			int_state = state_annot_dict[state] # if not in the above from, it should be of the form that are consistent with the mneunomic
		except:
			print("state annotation is not in the right format: {}".format(state))
			exit(1)
	return 'E' + str(int_state)


def get_sample_genome_segmentaton(sample_bed_fn, segment_bed_fn, output_fn, state_annot_dict, ct_name):
	sample_df = pd.read_csv(sample_bed_fn, header = None, sep = '\t', index_col = None)
	sample_df = sample_df.sort_values([0,1])
	print('Done reading in sample_df')
	segment_df = pd.read_csv(segment_bed_fn, header = None, sep = '\t', index_col = None)  # a df with columns chrom, start_bp, end_bp, state such that the df is already storted by chrom and start_bp, this is the assumption that the code will make
	segment_df = segment_df.sort_values([0,1])
	print('Done reading in segment_df')
	sample_b = bed.BedTool.from_dataframe(sample_df)
	segment_b = bed.BedTool.from_dataframe(segment_df)
	result_b = sample_b.map(segment_b, c = 4, o = 'collapse')
	result_df = result_b.to_dataframe().rename(columns = {'name' : ct_name})
	result_df[ct_name] = result_df[ct_name].apply(lambda x: transform_state_into_desired_form(x, state_annot_dict))
	result_df.to_csv(output_fn, header = True, index = False, sep = '\t', compression = 'gzip')
	return 

def main():
	if len(sys.argv) != 6:
		usage()
	sample_bed_fn = sys.argv[1]
	helper.check_file_exist(sample_bed_fn)
	segment_bed_fn = sys.argv[2]
	helper.check_file_exist(segment_bed_fn)
	state_annot_fn = sys.argv[3]
	helper.check_file_exist(state_annot_fn)
	state_annot_dict = read_state_annot(state_annot_fn)
	output_fn = sys.argv[4]
	helper.create_folder_for_file(output_fn)
	ct_name = sys.argv[5]
	print('Done getting command line arguments')
	get_sample_genome_segmentaton(sample_bed_fn, segment_bed_fn, output_fn, state_annot_dict, ct_name)
	print('Done!')

def usage():
	print("python get_sample_genome_segmentation.py")
	print('sample_bed_fn: where we save the coordinates of regions that we picked as the training regions')
	print('segment_bed_fn: a biosample\'s segmentation data')
	print('state_annot_fn: information about characteristics of states')
	print('output_fn')
	print('ct_name: the ct code that we will use as column name for the output dataframe')
	exit(1)

main()