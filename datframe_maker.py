import os
import pandas as pd
import numpy as np
import argparse
import sys
import warnings
warnings.filterwarnings("ignore")


"""
Useage: python <sample_id> <path> <mutation_file>
# Example: python drug_df_translation.py /Users/alvajames/Mount_temp/tupro_data_frames/data/sample_ds /Users/alvajames/Mount_temp/tupro_data_frames/data /Users/alvajames/Mount_temp/tupro_data_frames/data/fmidata_summary_formatted.txt
Author: Alva James, July 25, 2019
Purpose:
    The program converts mutiple dataframes from drug matrix and mutation file into the MOFA acceptable format
"""

parser = argparse.ArgumentParser()
parser.add_argument("sample_id", help=" The Sample_ids")
parser.add_argument("path", help="The path to drug matries")
parser.add_argument("mutation_file", help="Your input muation file")

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(0)

args = parser.parse_args()
# The inputs
sample_id     = args.sample_id
path          = args.path
mutation_file = pd.read_csv(args.mutation_file,sep="\t",header=0)

# Reading the sample_list  as a list and removeing newline
sample_list= [line.rsplit('\n') for line in open(sample_id)]
all_files  = [os.path.join(path,i) for i in os.listdir(path) if i.endswith('__DrugResponse_matrix.txt')]

# Defining an empty datframe
main       = []

# Defining the column of interest in this case its is All melonoma samples for both drug and mutation datasets
cols       = ['AllMel']

# parsing through files and appending the dataframes into previously defined empty dataframe main
# Renaming AllMel to the filename as header name
# Spliting the headers
# Making sure all samples in sample_id list are there in the new DataFrame
# Adding NA for empty values
# Changing the index name to drugs
# Saving the file

for filename in all_files:
    name=os.path.basename(os.path.normpath(filename))
    dfs = pd.DataFrame.from_csv(filename,sep='\t')
    dfs =dfs[[col for col in cols if col in dfs.columns]]
    dfs.rename(columns={'AllMel':name},inplace=True)
    dfs.rename(columns=lambda x:x.split('-F')[0],inplace=True)
    main.append(dfs)
merged_df = pd.concat(main,axis=1)
merged_df1= merged_df.replace(np.nan,'NA',regex=True)
merged_df1= merged_df1.reindex(columns=sample_list,fill_value='NA')
#merged_df1= merged_df1[sample_list]
merged_df1.index.name ="Drug_id"
merged_df1.to_csv(os.path.join(path,"Drug_dataframe"),sep='\t',index=True)

## Converting the mutation table into the dataframe as we need using pivot tables
mutation_dfs=mutation_file[["SAMPLE_ID","GENE"]]
mutation_dfs=mutation_dfs.pivot_table(index='GENE',columns='SAMPLE_ID',fill_value=0,aggfunc='size')
mutation_dfs=mutation_dfs.reindex(columns=sample_list,fill_value=0)
mutation_dfs.to_csv(os.path.join(path,"Mutation_dataframe"),sep="\t",index=True)
print("  Successfully finished making dataframes!!  ")
