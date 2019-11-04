
# coding: utf-8

# In[ ]:


#!/usr/bin/env
import os
import h5py    
import numpy as np
import pandas as pd
import sys
import argparse

__author__ = 'Alva James'

"""
Purpose:  
To convert the hd formated from single cell seq files to RNA-seq expression matrix
Conversion is mainly by taking column-wise sum of the cells expression from single-cell seq datasets in a given list of h5 files in the input
Finally, convert that into a dataframe or matrix in tsv format and save the files
"""

parser = argparse.ArgumentParser()
parser.add_argument("input",help="The path to the hd5 files")
parser.add_argument("output",help="The path to where you wanna write your input files")

if len(sys.argv)<2:
    parser.print_help()
    sys.exit(0)
    
args   = parser.parse_args()
path   = args.input
output = args.output

# Making the list of files for the h5 files in the input location
hdfiles = [os.path.join(path,i) for i in os.listdir(path) if i.endswith('__raw.h5')]

def main(hdfiles):
    """
    The function does the conversion and column-wise sum of cell expression values from single cell datasets
    """
    for files in hdfiles:
        name    = os.path.basename(os.path.normpath(files))
        sampleid=name.split('-')[0]
        raw     =h5py.File(files,'r')
        rawcounts=raw['raw_counts']
        cells =raw['cell_attrs']
        genes =raw['gene_attrs']
        rawcountsA =np.array(rawcounts)
        rawcountsT=rawcountsA.transpose()
        gene_ids=np.array(genes['gene_ids'])
        mydataframe=pd.DataFrame(data=rawcountsT,index=gene_ids,columns=cells['cell_names'])
        mydataframe["cell_sum"]=mydataframe.sum(axis=1,skipna=True)
        newdf = mydataframe[["cell_sum"]]
        newdf.reset_index(level=0,inplace=True)
        newdf.rename(columns={'index':'gene'},inplace=True)
        newdf['gene']=newdf['gene'].str.decode('utf-8')
        newdf['cell_sum'] =newdf['cell_sum'].astype(int)
        newdf.to_csv(os.path.join(output,sampleid+".tsv"),sep='\t',index=False,header=None)
    return 0


if __name__ ="__main__"
    main(hdfiles)
        
        

