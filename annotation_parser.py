#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
from argparse import ArgumentParser, SUPPRESS
import sys
import os


"""
Purpose:    This program parses the GENCODE annoation file
Author:     Alva James
Date:       October 8, 2019
Usage:      python annotation_parser.py <-p paralog/optional> <gtf> <option> <output> <organism>
example for mouse: python annotation_parser.py --paralog mouse_paralogs.GRCm38.p6.biomart.tsv Mouse_vM22_annotation_altered2.gtf CDS /Users/alvajames/Personal_offcial/Olga_task_Sep mouse
"""

parser = argparse.ArgumentParser()
parser.add_argument('--paralog',help='You can provide the paralog file, if exists')
parser.add_argument('GTF', help='Can be human or mouse')
parser.add_argument('type', help='String, eg. exon or CDS')
parser.add_argument('output', help='Path you need the output')
parser.add_argument('organism', help='which organism are you using, options: human or mouse')


args = parser.parse_args()
 
if len(sys.argv) < 4:
    parser.print_help()
    sys.exit(0)
    
option     = args.type
organism   = args.organism
outputpath = args.output
version    = args.GTF.split('_')[1]
gtffile    = pd.read_csv(args.GTF,sep='\t',header=0,skipinitialspace=True)

def main(df,paralogs=None, option=None, organism=None,*args):
    """
    The function to convert GFT file to tsv file based on exon or CDS positions
    The exon or CDS positions are collected as comma seperated list or values in correposnding new column 
    in the tsv file
    The approach here is making list of intronic start and end postions based on exonic start and end positions
    The same process will be repeated for the CDS
    The function will take care of if or not the paralog is given, currently the main function assumes only mouse has a paralog
    """
    NaN ='NaN'
    df['Gene']          = df.apply(
                            lambda row: NaN if row['gene'] !="transcript" else row['gene_symbol'],axis=1
                                  )
    df                  = df.replace('NaN',np.NaN)
    def filter_transcripts_with_multi_exons(df):
        """
        Filter the Genes which has transcripts with less than 2 exons
        """
        df["isagene"]           = np.where(df['gene'].isin(['transcript','exon']), df['Gene'].ffill(), 'no')
        exon_gene_count         = df.groupby('isagene').size()
        exon_gene_count1        = pd.DataFrame(exon_gene_count)
        exon_gene_count1.reset_index(inplace=True)
        exon_gene_count1.columns=["Gene","exon_number"]
        exon_gene_count1        = exon_gene_count1.query('exon_number <=2')
        Gene_listwithsingle_exon= exon_gene_count1.Gene.tolist()
        return Gene_listwithsingle_exon
    
    def filtering(df):
        """
        Remove all genes without a gene symbol based on the defined list genes_exclusion
        """
        genes_exclusion   = ['gene_status KNOWN','NaN','nan','gene_status PUTATIVE','gene_status NOVEL']
        gene_list         = df.Gene.drop_duplicates().tolist()
        gene_list_cleaned = [x for x in gene_list if str(x)  != tuple(genes_exclusion)]
        return gene_list_cleaned
    
    Gene_listwithsingle_exon = filter_transcripts_with_multi_exons(df)
    gene_list_cleaned        = filtering(df)
    
    if paralogs is not None:
        mouse_with_paralogs    = pd.merge(df[["gene_id","Gene"]],paralogs[["gene_id"]],on='gene_id')
        mouse_with_paralogs.drop_duplicates(inplace=True)
        paralog_mouse_gene_list= mouse_with_paralogs.Gene.tolist()
        def paralogs_func(row,paralog_mouse_gene_list):
            if row['Gene'] in (paralog_mouse_gene_list):
                return 1
            else:
                return 0
    if option =="exon":
        df         = df.sort_values(['chr','start','end','gene','transcript_symbol','gene_symbol'], ascending=True)
        df3        = df[df['gene']=='exon'].copy()
        df3['s']   = df3.groupby(['chr','gene_symbol', 'transcript_symbol'])['end'].shift()

        df3        = df3.dropna(subset=['s']).assign(exon_start = lambda x: x['s'].astype(int).astype(str),
                                      exon_end = lambda x: x['start'].astype(str))
        df4        = df3.groupby(['chr','gene_symbol', 'transcript_symbol'])['exon_start','exon_end'].agg(','.join).reset_index()
        df5        = df[df['gene']=='transcript'].drop_duplicates(['chr','gene_symbol', 'transcript_symbol']).drop('gene', 1)
        result     = df5.merge(df4)
        result["exon_start"] = result['exon_start'].apply(
                       lambda x: ','.join(x.split(',')))
        result["exon_end"]   = result['exon_end'].apply(
                       lambda x: ','.join(x.split(',')))
        result                 = result[result['gene_symbol'].isin(gene_list_cleaned)]
        result                 = result[~result['Gene'].isin(Gene_listwithsingle_exon)]
        result['exon_start']   = result['exon_start'].apply(lambda x: str(x)+',')
        result['exon_end']     = result['exon_end'].apply(lambda x: str(x)+',')
        if  paralogs is not None:
            print("Given annotation is mouse, position position looking for is",option,"----")
            result["paralog"]       = result.apply(lambda row: paralogs_func(row,paralog_mouse_gene_list), axis=1)
            colnames                = ['Gene','paralog','chr','strand','start','end','exon_end','exon_start','transcript_symbol','transript_id']
            result.rename(columns   = {'exon_end':'exon_start','exon_start':'exon_end'},inplace=True)
            result                  = result[colnames]
        else:
            result["paralog"]       = 0
            colnames                = ['Gene','paralog','chr','strand','start','end','exon_end','exon_start','transcript_symbol']
            result.rename(columns   = {'exon_end':'exon_start','exon_start':'exon_end'},inplace=True)
            result                  = result[colnames]
    elif option =="CDS":
        df         = df.sort_values(['chr','start', 'end','gene','transcript_symbol','gene_symbol'], ascending=True)
        df3        = df[df['gene']=='CDS'].copy()
        df3['s']   = df3.groupby(['chr','gene_symbol', 'transcript_symbol'])['end'].shift()
        df3        = df3.dropna(subset=['s']).assign(CDS_start = lambda x: x['s'].astype(int).astype(str),
                                      CDS_end = lambda x: x['start'].astype(str))
        df4        = df3.groupby(['chr','gene_symbol', 'transcript_symbol'])['CDS_start','CDS_end'].agg(','.join).reset_index()
        df_group   = df[df['gene']=='CDS'].groupby(['chr','gene_symbol', 'transcript_symbol','strand','transript_id'])
        df_new     = pd.DataFrame(columns =['start', 'end'])
        df_new[['start', 'end']]   = df_group.agg({'start':'first', 'end': 'last'})
        df_new.reset_index(inplace = True)
        df_group.head()
        df_new     = pd.DataFrame(columns =['start', 'end'])
        df_new[['start', 'end']]   = df_group.agg({'start':'first', 'end': 'last'})
        df_new.reset_index(inplace = True)
        result     = df_new.merge(df4)
        result.rename(columns={'gene_symbol':'Gene'},inplace=True)
        result["CDS_start"] = result['CDS_start'].apply(
                       lambda x: ','.join(x.split(',')))
        result["CDS_end"]   = result['CDS_end'].apply(
                       lambda x: ','.join(x.split(',')))
        result                 = result[result['Gene'].isin(gene_list_cleaned)]
        result                 = result[~result['Gene'].isin(Gene_listwithsingle_exon)]
        result['CDS_start']    = result['CDS_start'].apply(lambda x: str(x)+',')
        result['CDS_end']      = result['CDS_end'].apply(lambda x: str(x)+',')
        if  paralogs is not None:
            result["paralog"]      = result.apply(lambda row: paralogs_func(row,paralog_mouse_gene_list), axis=1)
            colnames                = ['Gene','paralog','chr','strand','start','end','CDS_end','CDS_start','transcript_symbol','transript_id']
            result.rename(columns   = {'CDS_end':'CDS_start','CDS_start':'CDS_end'},inplace=True)
            result                  = result[colnames]
        else:
            result["paralog"]       = 0
            colnames                = ['Gene','paralog','chr','strand','start','end','CDS_end','CDS_start','transcript_symbol']
            result.rename(columns   = {'CDS_end':'CDS_start','CDS_start':'CDS_end'},inplace=True)
            result                  = result[colnames]

    return(result)

if __name__ == "__main__":
    if organism    == "mouse":
        print("Input organism is: ",organism,  "searching for", option)
        paralogs   = pd.read_csv(args.paralog,sep='\t',header=0,skipinitialspace=True)
        result     = main(gtffile,paralogs=paralogs,option=option)
        result.to_csv(os.path.join(outputpath,organism+"_"+version+"_"+option+".tsv"),sep="\t",index=False)
    elif organism == "human":
        print("Input organism is: ",organism, "searching for", option)
        result     = main(gtffile,option=option)
        result.to_csv(os.path.join(outputpath,organism+"_"+version+"_"+option+".tsv"),sep="\t",index=False)