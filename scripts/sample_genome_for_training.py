# tested 02/25/2021
import pandas as pd 
import numpy as np 
import sys
import os
import helper 


def sample_genome_positions(chrom_length_fn, sample_fraction, output_fn):
	NUM_BP_PER_BIN = 200
	chr_length_df = pd.read_table(chrom_length_fn, sep = '\t', header = None)
	chr_length_df = chr_length_df[[0,2]] # get rid of the second column because they are all 1's
	chr_length_df.columns = ['chrom', 'length'] 
	chr_length_df.drop(len(helper.CHROMOSOME_LIST), axis = 0) # drop the last row because it's chromosome Y
	# chr_length_df has the following columns: chrom, length, end_bp, start_bp, start_bin, end_bin --> data about the chromosome
	sample_df = pd.DataFrame(columns = ['chrom', 'start_bp', 'end_bp']) # all bp index are zero-based --> a dataframe of all the positions that we will sample
	for chrom_index in helper.CHROMOSOME_LIST:
		this_chrom_length = (chr_length_df[chr_length_df['chrom'] == 'chr' + chrom_index])['length']
		num_bin_this_chrom = int(this_chrom_length / NUM_BP_PER_BIN )
		num_bin_to_sample = int(num_bin_this_chrom * sample_fraction) 
		sample_bins_this_chrom = np.random.choice(range(num_bin_this_chrom), size = num_bin_to_sample, replace = False)
		sample_bins_this_chrom.sort()
		sample_this_chrom_df = pd.DataFrame() 
		sample_this_chrom_df['chrom'] = ['chr' + chrom_index] * len(sample_bins_this_chrom)
		sample_this_chrom_df['start_bp'] = sample_bins_this_chrom * NUM_BP_PER_BIN
		sample_this_chrom_df['end_bp'] = sample_this_chrom_df['start_bp'] + NUM_BP_PER_BIN
		sample_df = sample_df.append(sample_this_chrom_df)
		print ("Done with chromosome: " + chrom_index)
	# save to file 
	sample_df.to_csv(output_fn, sep = '\t', header = None, index = False, compression = 'gzip')
	return sample_df


def main():
	if len(sys.argv) != 4:
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
	helper.create_folder_for_file(output_fn)
	print ("Done getting command line arguments")
	# select regions on the genome that we will sample from
	genome_sample_df = sample_genome_positions(chrom_length_fn, sample_fraction, output_fn) # --> a dataframe of 3 columns: "chromosome", "start_bp", 'end_bp'

	
def usage():
	print ("python sample_genome.py")
	print ("chrom_length_fn: a bed file with 3 columns, no headers. Columns should correspond to: chromsoome (chr1, chr2, etc.), start_bp (0 in all chromsomes), end_bp (the length of the chromosome, which will be a multiple of 200 because it's the resolution of the chromatin state annotation)")
	print ("sampling fraction: the fraction of the genome that we want to sample. Remember fractions are not percentages.")
	print ("output_fn: where the data of sampled regions and state segementation will be stored for all the cell types that we chose")
	print ("The result should give us around 1518140 200-bp bins")
	exit(1)
main()