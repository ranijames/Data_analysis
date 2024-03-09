import os
import h5py    
import numpy as np
import pandas as pd
from itertools import chain
import argparse
import sys
from argparse import ArgumentParser


"""
Usage : python samples_df_generator.py sample_ids
Author: Alva James, September 17 2019
"""
parser = argparse.ArgumentParser()
parser.add_argument("sample_ids", help="The file with absolute path to the new line seperated list of sample_ids. Each sample_id should start with the BSSE string")
parser.add_argument("read_pos", help="The read 1 and 2 position in sample id string, for eg: if sample id is MA_1_CTG_TTG_R1.fastq.gz the read position is 4")
parser.add_argument("sample_pos", help="The sample name position in sample id string, for eg: if sample id is MA_1_CTG_TTG_R1.fastq.gz the sample name position is 0")
parser.add_argument("output", help="The path where you need your output samples.tsv file")

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(0)

args = parser.parse_args()
samples            = args.sample_ids
read_pos           = args.read_pos
sample_pos         = args.sample_pos
output             =  args.output

#read_pos=10
#sample_pos=0

def opening_samples_as_list(samples,l1):
    """
    Function to open the given set of samples as a list
    """
    sams=open("samples","r")
    for lines in sams.readlines():
        for samples in lines.strip().split(" "):
            l1.append(samples)

l1=[]
opening_samples_as_list(samples,l1)

# extracting the string for sample_id, and fatsqs from reads one and two
col1     = [x.split('_')[sample_pos] for x in l1]
col3     = [x for x in l1 if x.split('_')[read_pos] == "R1"]
col4     = [x for x in l1 if x.split('_')[read_pos] == "R2"]

# Converting the above list to dictionaries with sample_id as there key

dictionary_fq1    = dict((x,[]) for x in col1)
dictionary_fq2    = dict((x,[]) for x in col1)
dictionary_lane   = { y : [1,2] for y in col1}

# Populating teh dictionaries for each fastqs
for id2 in dictionary_fq2:
    for names in col4:
        if id2 in names:
            dictionary_fq2[id2].append(names)
            
for id1 in dictionary_fq1:
    for names in col3:
        if id1 in names:
            dictionary_fq1[id1].append(names)
            
# Converting the dictionaries into dataframe  of interest
df = pd.DataFrame({'lane': pd.DataFrame(dictionary_lane).unstack(),
              'fq1': pd.DataFrame(dictionary_fq1).unstack(),
              'fq2': pd.DataFrame(dictionary_fq2).unstack()}).reset_index(level=0)
df.rename(columns={'level_0' :'sample'}, inplace=True)

# saving the above generated datafarme as samples.tsv file which can be used for SNAKEMAKE worfkflow
df.to_csv(os.path.join(output,"samples.tsv"),sep="\t",index=True)
