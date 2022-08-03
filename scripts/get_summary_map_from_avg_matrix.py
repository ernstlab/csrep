#!/usr/bin/env python
import pandas as pd 
import numpy as np 
import os
import glob
import helper
import argparse
parser = argparse.ArgumentParser(description = 'Given the summary chromatin state assignment matrix the CSREP/base_count calculates (outputted by average_pred_results.py), this code will output the summary chromatin state map for a group of samples. The output will be a bed-format file that can be readable in the ucsc genome browser (if the user sets the flag --igv-format to 1).')
parser.add_argument('--avg_folder', type = str, required = True,
	help = 'Where there are files showing the representative chroamtin state assignment matrices for different regions on the genome')
parser.add_argument('--output_fn', type = str, required = True, 
	help = 'output_fn')
parser.add_argument('--igv_format', type = int, choices = [0,1], required = False, default = 1,
		help = '1 means that output summary chromatin state map will be written in a form that can be read into ucsc genome browser. 0 the the output summary chromatin state map will just be a normal bed file with columns: chrom, start, end, state (ex: E1 --> E18)') # had to set this to either 0 or 1 instead of being a boolean value right  away because then I can run this on my snakemake pipeline
parser.add_argument('--igv_track_name', type=str, required=False,
	default = 'csrep_summary',
	help = 'Name of the track on ucsc genome browser should you want to upload the output file to genome browser for visualization. For example, the name of sample group that you are summarizing.')
parser.add_argument('--state_annot_fn', type = str, required = False,
	default = 'test_data/state_annotation.txt',
	help = 'If you dont provide the state_annot_fn, we will use the default which is test_data/state_annotation.txt, which is from roadmap\'s 18-state annotation. If your annotation is not from the same model, please provide this parameter otherwise the code will stop without producing output.')
args = parser.parse_args()
print(args)
helper.check_dir_exist(args.avg_folder)
helper.create_folder_for_file(args.output_fn)
helper.check_file_exist(args.state_annot_fn)
igv_format = True
if (args.igv_format == 1):
	igv_format = True
elif args.igv_format == 0:
	igv_format = False
else: 
	print ('igv_format needs to be 0 or 1')
	exit(1)


def read_state_annot_fn(state_annotation_fn):
	df = pd.read_csv(state_annotation_fn, header = 0, sep = '\t')
	df['state'] = list(map(lambda x: 'E' + str(x), df['state']))
	return df

def get_rgb_format_right(rgb):
	# convert from (255, 245, 238) to 255,245,238
	# numbers = (rgb[1:-1]).split(',') # get rid of () 
	numbers = (rgb[:]).split(',')
	numbers = list(map(lambda x: str(int(x)), numbers))
	return ",".join(numbers)

def get_standard_ucsc_format_bed_file(segment_df, state_annot_df, output_fn, igv_track_name):
	segment_df = pd.merge(segment_df, state_annot_df
		, how = 'left', left_on = 'state', right_on = 'state', left_index = False, right_index = False)
	segment_df = segment_df [['chrom', 'start_bp', 'end_bp', 'state', 'itemRgb']]
	(nrow, ncol) = segment_df.shape
	segment_df['itemRgb'] = (segment_df['itemRgb']).apply(get_rgb_format_right)
	segment_df['score'] = ['1'] * nrow
	segment_df ['strand'] = ['.'] * nrow
	segment_df['thickStart'] = segment_df['start_bp']
	segment_df['thickEnd'] = segment_df['start_bp']
	segment_df = segment_df.rename(columns = {'start_bp' : 'chromStart', 'end_bp' : 'chromEnd', 'state' : 'name'})
	segment_df = segment_df[['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb']]
	outF = open(output_fn, 'a')
	header_comment = "track name=\"" + igv_track_name + "_" + "\" description=\"\" visibility=1 itemRgb=\"On\"\n"
	outF.write(header_comment) # write the comment first so that genome browser can read the file
	segment_df.to_csv(outF, sep = '\t', header = False, index = False)
	outF.close()


def get_max_prob_state_segmentation_one_region(state_rep_prob_fn, genome_pos):
	print (state_rep_prob_fn)
	segment_df = pd.read_csv(state_rep_prob_fn, header = 0, index_col = 0, sep = '\t')
	max_state = segment_df.idxmax(axis = 1) # get the column max, but return the column names
	max_state = list(map(lambda x: 'E' + str(x.split('_')[1]), max_state)) # chagne from state_18 to E18 for all elements in the max_state series
	chr_data = genome_pos.split('_')[0]
	offset_bp = int(genome_pos.split('_')[1])
	offset_bp = offset_bp * helper.NUM_BP_PER_WINDOW # each window contains data in terms of NUM_BP_PER_WINDOW
	result_df = pd.DataFrame(columns = ['chrom', 'start_bp', 'end_bp', 'state'])
	current_start_bin_index = 0
	current_end_bin_index = 1 
	current_state = max_state[current_start_bin_index]
	for index in range(1, len(max_state)): # skip the first one because it's already reported into the current data
		if (current_state == max_state[index]):
			current_end_bin_index = index + 1 # add one more bin to the stread of bins that are segmented here
		else: # change of state --> report into the result_df
			start_bp_report = offset_bp +  current_start_bin_index * helper.NUM_BP_PER_BIN
			end_bp_report = offset_bp + current_end_bin_index * helper.NUM_BP_PER_BIN
			add_row = [chr_data, start_bp_report, end_bp_report, current_state]
			result_df.loc[result_df.shape[0]] = add_row
			current_state = max_state[index]
			current_start_bin_index = index
			current_end_bin_index = index + 1
	start_bp_report = offset_bp + current_start_bin_index * helper.NUM_BP_PER_BIN
	end_bp_report = offset_bp + current_end_bin_index * helper.NUM_BP_PER_BIN
	add_row = [chr_data, start_bp_report, end_bp_report, current_state]
	result_df.loc[result_df.shape[0]] = add_row
	# now we are done producing
	return result_df


def create_igv_format_bed(avg_folder, state_annot_fn, output_fn, igv_format, igv_track_name):
	rep_prob_fn_list = glob.glob(avg_folder + '/*_avg_pred.txt.gz')
	genome_pos_list = list(map(lambda x: x.split('/')[-1].split('_avg_pred.txt.gz')[0], rep_prob_fn_list))
	state_annot_df = read_state_annot_fn(state_annot_fn)
	result_df = pd.DataFrame(columns = ['chrom', 'start_bp', 'end_bp', 'state'])
	for index, rep_prob_fn in enumerate(rep_prob_fn_list):
		this_region_max_state_df = get_max_prob_state_segmentation_one_region(rep_prob_fn, genome_pos_list[index])
		result_df = result_df.append(this_region_max_state_df)
	result_df.sort_values(['chrom', 'start_bp'], inplace = True)
	if igv_format: # if user specified they wanted output that can then be read into ucsc genome browser
		get_standard_ucsc_format_bed_file(result_df, state_annot_df, output_fn, igv_track_name)
	else:
		result_df.to_csv(output_fn, header = False, index = False, sep = '\t')
	print('Done!')
	return 


create_igv_format_bed(args.avg_folder, args.state_annot_fn, args.output_fn, igv_format, args.igv_track_name)