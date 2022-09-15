#!/usr/bin/env python

'''
This file contains some useful functions for file-managements, reading input data, etc. that can be used by other scripts
List of variables: 
CHROMOSOME_LIST: 1--22, X --> For CSREP, we only gets summary chromatin state maps for chromosomes 1--> 22, X, because not all our input samples have Y chromosome. 
NUM_BP_PER_BIN: 200 (bp)
NUM_BIN_PER_WINDOW: 50,000 (windows, each of length NUM_BP_PER_BIN)
NUM_BP_PER_WINDOW: NUM_BP_PER_BIN * NUM_BIN_PER_WINDOW

List of functions:
- make_dir(dir_path) --> create a directory (recursively) if it does not exist yet
- check_file_exist(fn) --> if not, exit the program
- create_folder_for_file(fn) --> usually used when fn is an output file path. This function will create the folder that contains the file fn
- check_dir_exist(dir_path) --> if not, exit the program
- get_command_line_integer(argument) --> try to convert to integer, if not succesfully then exit the program. This function is not entirely useful anymore given that we use argparse for all our scripts now. 
- get_list_from_line_seperated_file(fn) --> read in a file such that each line is an item in a list, using pandas series
- partition_file_list (file_list, num_cores) --> partition the file_list into a list of num_cores lists, as evely distributed as possible. Useful when we want to partition the list of outptu files into lists that can then be divided for different cores  to produce in parallel. 
'''
import os
import numpy as np
import pandas as pd

CHROMOSOME_LIST = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X'] # we exclude chromosome Y because it's not available in all cell types

NUM_BP_PER_BIN = 200
NUM_BIN_PER_WINDOW = 50000 # each window is 10 000 000 base pairs
NUM_BP_PER_WINDOW = NUM_BP_PER_BIN * NUM_BIN_PER_WINDOW 

def make_dir(directory):
	try:
		os.makedirs(directory)
	except:
		print('Folder' + directory + ' is already created')



def check_file_exist(fn):
	if not os.path.isfile(fn):
		print("File: " + fn + " DOES NOT EXIST")
		exit(1)
	return True
	
def create_folder_for_file(fn):
	last_slash_index = fn.rfind('/')
	if last_slash_index != -1: # path contains folder
		make_dir(fn[:last_slash_index])
	return 

def check_dir_exist(dirName):
	if not os.path.isdir(dirName):
		print("Directory: " + dirName + " DOES NOT EXIT")
		exit(1)
	return

def get_command_line_integer(arg):
	try: 
		arg = int(arg)
		return arg
	except:
		print("Integer: " + str(arg) + " IS NOT VALID")
		exit(1)

def get_list_from_line_seperated_file(fn):
	# from a text file where each line contains an item, we will get a list of the items
	result =  list(pd.read_csv(fn, sep = '\n', header = None)[0]) # -->  a list with each entry being an element in a list. Note the [0] is necessary for us to get the first column
	return result

def partition_file_list(file_list, num_cores):
	results = [] # list of lists of file names
	num_files_per_core = int(np.around(len(file_list) / num_cores)) # round to the nearest decimal point, either up or down
	if num_files_per_core == 0:
		for file in file_list:
			results.append([file])
		for core_i in range(len(file_list), num_cores):
			results.append([])
		return results
	# else, the number of files per core is greater than 0, meaning only
	for core_i in range(num_cores):
		if core_i < (num_cores - 1):
			this_core_files = file_list[core_i * num_files_per_core : (core_i + 1) * num_files_per_core]
		elif core_i == (num_cores - 1):
			this_core_files = file_list[core_i * num_files_per_core :]
		results.append(this_core_files)
	return results
