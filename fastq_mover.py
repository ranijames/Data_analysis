import os
import pandas as pd
import fnmatch
import shutil, sys
from os import listdir
from os.path import isfile, join
import argparse
import sys
import warnings
warnings.filterwarnings("ignore")


"""
Useage: python fastq_mover.py <fastq_dir> <dropbox>
Example: python fastq_mover.py /cluster/work/grlab/projects/tumor_profiler/data/fastq /cluster/work/grlab/projects/tumor_profiler/bkRNA_samples/dropbox
Author: Alva James, September 24, 2019
Purpose:
    The program converts mutiple dataframes from drug matrix and mutation file into the MOFA acceptable format
"""

parser = argparse.ArgumentParser()
parser.add_argument("fastq_dir", help=" The Sample_ids")
parser.add_argument("dropbox", help="The path to drug matries")
args   = parser.parse_args()

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(0)

fastqpath   = args.fastq_dir
dropbox     = args.dropbox
fastq_files = [os.path.join(fastqpath,i).split('/')[-1] for i in os.listdir(fastqpath) if i.startswith("BSSE")]
#print(fastq_files)

# Making list of samples and converting the list to dataframe and save it in the fastq folder
sam_lis      =list()
for sam in fastq_files:
    sam_list =sam.split('_')[5]
    sam_lis.append(sam_list)
    sam_lis  =pd.unique(sam_lis).tolist()
sample_id    = pd.DataFrame({'sample_id':sam_lis })

# Using the above list
#sample_id   = sam_lis
sample_id   =['MODUDOL']
print(sample_id)
# Function for moving the raw fastqs from the fastq directory to the dropbox
def fastqmover(sample_id,fastq_files,dropbox):
    """
    Function to move fastq files from the common place to the destination folder
    """
    for samples in sample_id:
        for fastqs in fastq_files:
            if samples in fastqs:
                desination = dropbox + "/bkRNA/" + samples + "/raw/"
                if not os.path.isdir(desination):
                    os.makedirs(desination)
        for rawfiles in fnmatch.filter(fastq_files, pat="*"):
            if samples in rawfiles:
                shutil.move(os.path.join(fastqpath,rawfiles), os.path.join(desination,rawfiles))
                

fastqmover(sample_id,fastq_files,dropbox)


