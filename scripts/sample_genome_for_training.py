#!/usr/bin/env python

'''
This script will take in the file specifying the length of chromosomes, and user-specified sample_fraction. It will then sample the corresponding fraction of the genome, which will be used as training regions for CSREP. The output of this file is a bed-format file: chrom, start, end --> training genomic bins. Each line will correspond to a 200-bp bin for which we will sample. 

command-line arguments:
python sample_genome.py
chrom_length_fn: a bed file with 3 columns, no headers. Columns should correspond to: chromsoome (chr1, chr2, etc.), start_bp (0 in all chromsomes), end_bp (the length of the chromosome, which will be a multiple of 200 because it's the resolution of the chromatin state annotation)
sample_fraction: the fraction of the genome that we want to sample. Remember fractions are not percentages.
output_fn: where the data of sampled regions and state segementation will be stored for all the cell types that we chose
seed: random seed, for reproducibility
The result should give us around 1518140 200-bp bins
'''

# tested 02/25/2021
import pandas as pd 
import numpy as np 
import sys
import os
import helper 
NUM_BP_PER_BIN = 200


def get_chrom_length(chrom_length_fn):
	chr_length_df = pd.read_table(chrom_length_fn, sep = '\t', header = None)
	chr_length_df.columns = ['chrom', 'start_bp', 'end_bp'] # this file shows all the available segments in different chromosome that we can sample from.  
	chr_length_df['length'] = chr_length_df['end_bp'] - chr_length_df['start_bp']
	# chr_length_df has the following columns: chrom, length, end_bp, start_bp, start_bin, end_bin --> data about the chromosome
	chr_length_df = chr_length_df.groupby('chrom')['length'].sum() # a series with index: chr1 --> chrX/chrY, etc. 
	return chr_length_df.to_dict()

def sample_genome_positions(chrom_length_fn, sample_fraction, output_fn, seed):
	np.random.seed(seed) # set the random seed for reproducibility
	chrom_length_dict = get_chrom_length(chrom_length_fn) # keys: chr1--> chrX/chrY, values: length of the chromosomes
	sample_df = pd.DataFrame(columns = ['chrom', 'start_bp', 'end_bp']) # all bp index are zero-based --> a dataframe of all the positions that we will sample
	for chrom_index in helper.CHROMOSOME_LIST:
		try: 
			this_chrom_length = chrom_length_dict['chr{}'.format(chrom_index)]
			assert this_chrom_length > 0
		except:
			print ('Chromosome {} is skipped because we did not want to include this chrom'.format(chrom_index))
			continue
		num_bin_this_chrom = int(this_chrom_length / NUM_BP_PER_BIN )
		num_bin_to_sample = int(num_bin_this_chrom * sample_fraction) 
		sample_bins_this_chrom = np.random.choice(num_bin_this_chrom, size = num_bin_to_sample, replace = False)
		sample_bins_this_chrom.sort()
		sample_this_chrom_df = pd.DataFrame() 
		sample_this_chrom_df['chrom'] = ['chr' + chrom_index] * len(sample_bins_this_chrom)
		sample_this_chrom_df['start_bp'] = sample_bins_this_chrom * NUM_BP_PER_BIN
		sample_this_chrom_df['end_bp'] = sample_this_chrom_df['start_bp'] + NUM_BP_PER_BIN
		sample_df = sample_df.append(sample_this_chrom_df)
		print ("Done with chromosome: " + chrom_index)
	sample_df = sample_df.sort_values(['chrom', 'start_bp'])
	# save to file 
	sample_df.to_csv(output_fn, sep = '\t', header = False, index = False, compression = 'gzip')
	return sample_df


def main():
	if len(sys.argv) != 5:
		usage()
	chrom_length_fn = sys.argv[1]
	helper.check_file_exist(chrom_length_fn)
	sample_fraction = sys.argv[2]
	try: 
		sample_fraction = float(sample_fraction)
	except: 
		print("sample_fraction should be a float number")
		usage()
	assert sample_fraction > 0 and sample_fraction < 1.0, "sample_fraction should be greater than 0 and smaller than 1"
	output_fn = sys.argv[3]
	seed = helper.get_command_line_integer(sys.argv[4])
	helper.create_folder_for_file(output_fn)
	print ("Done getting command line arguments")
	# select regions on the genome that we will sample from
	genome_sample_df = sample_genome_positions(chrom_length_fn, sample_fraction, output_fn, seed) # --> a dataframe of 3 columns: "chromosome", "start_bp", 'end_bp'

	
def usage():
	print ("python sample_genome.py")
	print ("chrom_length_fn: a bed file with 3 columns, no headers. Columns should correspond to: chromsoome (chr1, chr2, etc.), start_bp (0 in all chromsomes), end_bp (the length of the chromosome, which will be a multiple of 200 because it's the resolution of the chromatin state annotation)")
	print ("sampling fraction: the fraction of the genome that we want to sample. Remember fractions are not percentages.")
	print ("output_fn: where the data of sampled regions and state segementation will be stored for all the cell types that we chose")
	print ('seed: random seed, for reproducibility')
	print ("The result should give us around 1518140 200-bp bins")
	exit(1)
main()