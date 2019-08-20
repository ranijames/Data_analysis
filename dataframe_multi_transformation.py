# For Protein datasets
import os
import pandas as pd
import numpy as np
import glob
import argparse
import sys

# Hide the warning from the user
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

"""
Useage: python protein_exp_dataframe_transformation.py <path> <sample_id>
Author: Alva James, August 20, 2019
Multi-Purpose script:
    1. The current program removes or strips characters like: square brackets, digits, white spaces from the
    column names of a given set of dataframe(s)

    2. Merge a given set of dataframes with unequal diamensions, based on condition:
       If the one of the col is equal to the file name
    3. Setting index for a dataframe using commonly occuring column, here Protein ID
    4. Protein ID was not able to set as index because it was not unique accross the samples, at the same time
    we needed all protein expression values and fill the missing samples with NA.
    5. Reshape the such a dataframe into long format
    6. The decimals where seperated using comma not by period, had to address this using read_csv
    7. Finall, adding all sample_ids from the provided list and filling in the missing with NA

"""

parser = argparse.ArgumentParser()
parser.add_argument("path", help="The path to drug matries")
parser.add_argument("sample_id", help=" The Sample_ids")


if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(0)

args = parser.parse_args()


path           = args.path
sample_ids     = args.sample_id

sample_list = [line.rstrip('\n') for line in open(sample_ids)]
pc_norm     = [os.path.join(path,i) for i in os.listdir(path) if i.endswith('Normalized_Protein_Report.tsv')]

for filename in pc_norm:
    name       = os.path.basename(os.path.normpath(filename))
    name_list  = name.split('-')[0]
    dfs        = pd.DataFrame.from_csv(filename, sep="\t")
    dfs.rename(columns = lambda x:x.split('.PG.Quantity')[0],inplace=True)
    dfs        = dfs.rename(columns=lambda x: x.replace('[','').replace(']',''))
    dfs.columns= dfs.columns.str.replace('\d+ ', '')
    dfs.columns= dfs.columns.str.replace(' ', '')
    dfs.rename(columns=lambda x:x.split('-')[0],inplace=True)
    dfs.drop_duplicates(inplace=True)
    dfs.to_csv(os.path.join(path,name_list+"_"+"protein_modified"),sep='\t',index=True)


pc_modified            = [os.path.join(path,i) for i in os.listdir(path) if i.endswith('protein_modified')]
pc_modified_files      = [os.path.basename(os.path.normpath(files)) for files in pc_modified]
filename_short         = [filenam.split('-')[0] for filenam in pc_modified_files]
main           = []
colnames       = []
main           = []

# Choosing the column with which has same name as the filename, in order to make sure we are capturing Protein ID and
# its expression column for each sample_id

for files in pc_modified:
    name = os.path.basename(os.path.normpath(files))
    filename_short  = name.split('_')[0]
    #print(filename_short,files)
    try:
        df = pd.read_csv(files, usecols=[filename_short] + ['PG.ProteinAccessions'], index_col=['PG.ProteinAccessions'], sep='\t')
        df.reset_index(drop=False, inplace=True)
        main.append(df)
    except:
        pass

merged = pd.concat(main,axis=1,sort=True)

merged.drop_duplicates(subset='PG.ProteinAccessions', keep='first',inplace=True)

df2 = pd.DataFrame()

# Converting the above dataframe to long format
for i in range(len(merged.columns)//2):
    key          = merged.columns[2*i+1]
    dfx          = pd.DataFrame()
    dfx['ID']    = merged.iloc[:,2*i]
    dfx['Key']   = key
    dfx['Value'] = merged.iloc[:,2*i+1]
    df2          = pd.concat([df2,dfx], sort=False)

df2.to_csv(os.path.join(path,"protein_modified.tsv"),sep='\t')

pc          = pd.read_csv(os.path.join(path,"protein_modified.tsv"), sep = '\t', decimal=',')
pc['Value'] = pc['Value'].convert_objects(convert_numeric = True)
pc_final    = pc.reset_index().pivot_table(values='Value',index='ID',columns='Key')
pc_final    = pc_final.reindex(columns=sample_list,fill_value='NA')
pc_final    = pc_final[[s for s in sample_list if s in pc_final.columns ]]
pc_final.index.name = "Protein_ID"
pc_final.fillna('NA',inplace=True)

pc_final.to_csv(os.path.join(path,"Protein_dataframe.tsv"),sep='\t',index=True)

# deleting the intermediate files generated
for fname in os.listdir(path):
    if fname.endswith("protein_modified"):
        os.remove(os.path.join(path, fname))

print("Successfully finished transforming the protein expression files into MOFA format !!")
