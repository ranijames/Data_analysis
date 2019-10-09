import os
import pandas as pd
import numpy as np
import glob
import fnmatch
import shutil, sys
import argparse
from argparse import ArgumentParser

"""
Usage : python data_management.py <sampl_lab> <results> <dropbox>
Author: Alva James

"""
parser = argparse.ArgumentParser()
parser.add_argument("sampl_lab", help="The samples and labkey as tab seperated dataframe")
parser.add_argument("results", help="The absolute path to the SNAKEMAKE generated the outputs")
parser.add_argument("dropbox", help="The path to the tumor profiler dropbox")

if len(sys.argv) < 3:
    parser.print_help()
    sys.exit(0)

args               = parser.parse_args()
sampl_lab          = pd.read_csv(args.sampl_lab)
results            = args.results
dropbox            = args.dropbox

sample = sampl_lab['sample'].tolist()
samp_dic                = dict((samples,[]) for samples in sample)
#print(samp_dic)
tp_results              = ['alignments','expression','hypoxia','tcga_boxplot/png_images','qc/fastqc']
snakemake_folders       = [dirs for dirs in os.listdir(results)]
snakemake_subfolders    = [os.path.join(results, name) for name, sub, files in os.walk(results) if os.path.isdir(os.path.join(results,name))]
snakemake_subfolders_filtered = [s for s in snakemake_subfolders if s.endswith(tuple(tp_results))]
snakemake_files         = [os.path.join(results,name) for root, sub, files in os.walk(results) for name in files]
extensions              = ['all.bam','fastqc.zip','.expression.tsv','tcga_cohort.png','hypoxia_score.tsv']

file_paths = []
for folder, subs, files in os.walk(results):
    for filename in files:
        file_paths.append(os.path.abspath(os.path.join(folder, filename)))
files_filtered = [f for f in file_paths if f.endswith(tuple(extensions))]
for filetred in snakemake_subfolders_filtered:
    for files in files_filtered:
        if filetred in files:
            for sample,lab_key in samp_dic.items():
                if sample in files:
                    desination = dropbox + "/bkRNA/" + sample + "/derived/"
                    if not os.path.isdir(desination):
                            os.makedirs(desination)
                    for desination_files in fnmatch.filter(files_filtered, pat="*"):
                         if sample in desination_files:
                                try:
                                    #print(os.path.join(results,desination_files),os.path.join(desination,os.path.basename(os.path.normpath(desination_files)))
                                    os.symlink(os.path.join(results,desination_files),os.path.join(desination,os.path.basename(os.path.normpath(desination_files))))
                                    break
                                except FileExistsError:
                                    pass
