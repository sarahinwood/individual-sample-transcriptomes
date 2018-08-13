!#/bin/env python3

import pathlib2
import pandas
import os

#############
# FUNCTIONS #
#############

def find_read_files(read_dir):
#Make list of files
	path_generator = os.walk(read_dir, followlinks = True)
	my_files = list((dirpath, filenames)
		for (dirpath, dirname, filenames)
		in path_generator)
#Make new empty dictionary & populate with files (flowcell = key)
	my_fastq_files = {}
	for dirpath, filenames in my_files:
		for filename in filenames:
			if filename.endswith('.fastq.gz'):
				my_flowcell = pathlib2.Path(dirpath).name
				my_fastq = str(pathlib2.Path(dirpath,filename))
				if my_flowcell in my_fastq_files:
					my_fastq_files[my_flowcell].append(my_fastq)
				else:
					my_fastq_files[my_flowcell]= []
					my_fastq_files[my_flowcell].append(my_fastq)
	return(my_fastq_files)

###########
# GLOBALS #
###########

read_dir = 'data/reads'

sample_key_file = 'data/sample_key.csv'

#########
# SETUP #
#########

sample key = pandas.read_csv(sample_key_file)

# generate name to filename dictionary
all_fastq = find_read_files(read_dir)

#########
# RULES #
#########

