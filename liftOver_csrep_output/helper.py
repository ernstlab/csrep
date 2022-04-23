import string
import os
import sys
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
	num_files_per_core = int(len(file_list) / num_cores)
	for core_i in range(num_cores):
		if core_i < (num_cores - 1):
			this_core_files = file_list[core_i * num_files_per_core : (core_i + 1) * num_files_per_core]
		elif core_i == (num_cores - 1):
			this_core_files = file_list[core_i * num_files_per_core :]
		results.append(this_core_files)
	return results
