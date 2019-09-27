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
Date:       September, 2019
Usage:      python GTF_converter.py <-p paralog/optional> <gtf> <option> <output> <organism>
"""

parser = argparse.ArgumentParser()
parser.add_argument('GTF', help='Can be human or mouse')
parser.add_argument('option', help='String, eg. exon or CDS')
parser.add_argument('output', help='Path you need the output')
parser.add_argument('organism', help='which organism are you using, options: human or mouse')
parser.add_argument('-p','--paralog',help='You can provide the paralog file, if exists')

args = parser.parse_args()
 
if len(sys.argv) < 3:
    parser.print_help()
    sys.exit(0)
    
colnames   = ['chr','gene', 'start','end','strand','gene_id','gene_symbol','transcript_symbol','transript_id'] 
gtffile    = pd.read_csv(args.GTF,sep='\t',names=colnames,skipinitialspace=True)
option     = args.option
paralog    = pd.read_csv(args.paralog,sep='\t',skipinitialspace=True)
organism   = args.organism
outputpath = args.output


def main(df,option=None, paralogs=None, *args):
    """
    The function to convert GFT file to tsv file based on exon or CDS positions
    The exon or CDS positions are collected as comma seperated list or values in correposnding new column 
    in the tsv file
    The approach here is making list of intronic start and end postions based on exonic start and end positions
    So as for the CDS positions
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
        Remove all genes without a gene symbol based on list gene exclusion
        """
        genes_exclusion   = ['gene_status KNOWN','NaN','nan','gene_status PUTATIVE','gene_status NOVEL']
        gene_list         = df.Gene.drop_duplicates().tolist()
        gene_list_cleaned = [x for x in gene_list if str(x)  != tuple(genes_exclusion)]
        return gene_list_cleaned
    
    Gene_listwithsingle_exon = filter_transcripts_with_multi_exons(df)
    gene_list_cleaned        = filtering(df)
    
    if paralogs is not None:
        #print("----- input is a mouse paralog, with",option," ----")
        mouse_paralogs         = paralogs
        mouse_with_paralogs    = pd.merge(df[["gene_id","Gene"]],mouse_paralogs[["gene_id"]],on='gene_id')
        mouse_with_paralogs.drop_duplicates(inplace=True)
        paralog_mouse_gene_list= mouse_with_paralogs.Gene.tolist()
        def paralogs(row,paralog_mouse_gene_list):
            if row['Gene'] in (paralog_mouse_gene_list):
                return 1
            else:
                return 0
    if option =="exon":
        #print("----- sequence type is", option," ----, begining to collect intron positions")
        df         = df.sort_values(['chr','start', 'end','gene','transcript_symbol','gene_symbol'], ascending=True)
        df3        = df[df['gene']=='exon'].copy()
        df3['s']   = df3.groupby(['chr','gene_symbol', 'transcript_symbol'])['end'].shift()

        df3        = df3.dropna(subset=['s']).assign(exon_start = lambda x: x['s'].astype(int).astype(str),
                                      exon_end = lambda x: x['start'].astype(str))
        df4        = df3.groupby(['chr','gene_symbol', 'transcript_symbol'])['exon_start','exon_end'].agg(','.join).reset_index()
        df5        = df[df['gene']=='transcript'].drop_duplicates(['chr','gene_symbol', 'transcript_symbol']).drop('gene', 1)
        result     = df5.merge(df4)
        result["exon_start"] = result['exon_start'].apply(
                       lambda x: ','.join(sorted(x.split(','))))
        result["exon_end"]   = result['exon_end'].apply(
                       lambda x: ','.join(sorted(x.split(','))))
        result                  = result[result['gene_symbol'].isin(gene_list_cleaned)]
        result                  = result[~result['Gene'].isin(Gene_listwithsingle_exon)]
        result['exon_start']   = result['exon_start'].apply(lambda x: str(x)+',')
        result['exon_end']     = result['exon_end'].apply(lambda x: str(x)+',')
        if  paralogs is not None:
            #print("Given annotation is mouse, position position looking for is",option,"----")
            result["paralog"]      = result.apply(lambda row: paralogs(row,paralog_mouse_gene_list), axis=1)
            colnames                = ['Gene','paralog','chr','strand','start','end','exon_end','exon_start','transcript_symbol','transript_id']
            result.rename(columns   = {'exon_end':'exon_start','exon_start':'exon_end'},inplace=True)
            result                  = result[colnames]
            print("--------Successfully completed generating tsv file for",organism,option, "--------")
        else:
            result["paralog"]       = 0
            colnames                = ['Gene','paralog','chr','strand','start','end','exon_end','exon_start','transcript_symbol']
            result.rename(columns   = {'exon_end':'exon_start','exon_start':'exon_end'},inplace=True)
            result                  = result[colnames]
            print("--------Successfully completed generating tsv file for",organism,option, "--------")
    elif option =="CDS":
        #print("----- sequence type is", option," ---- begining to collect intron positions")
        df         = df.sort_values(['chr','start', 'end','gene','transcript_symbol','gene_symbol'], ascending=True)
        df3        = df[df['gene']=='CDS'].copy()
        df3['s']   = df3.groupby(['chr','gene_symbol', 'transcript_symbol'])['end'].shift()
        df3        = df3.dropna(subset=['s']).assign(CDS_start = lambda x: x['s'].astype(int).astype(str),
                                      CDS_end = lambda x: x['start'].astype(str))
        df4        = df3.groupby(['chr','gene_symbol', 'transcript_symbol'])['CDS_start','CDS_end'].agg(','.join).reset_index()
        df5        = df[df['gene']=='transcript'].drop_duplicates(['chr','gene_symbol', 'transcript_symbol']).drop('gene', 1)
        result     = df5.merge(df4)
        result["CDS_start"] = result['CDS_start'].apply(
                       lambda x: ','.join(sorted(x.split(','))))
        result["CDS_end"]   = result['CDS_end'].apply(
                       lambda x: ','.join(sorted(x.split(','))))
        result                 = result[result['gene_symbol'].isin(gene_list_cleaned)]
        result                 = result[~result['Gene'].isin(Gene_listwithsingle_exon)]
        result['CDS_start']    = result['CDS_start'].apply(lambda x: str(x)+',')
        result['CDS_end']      = result['CDS_end'].apply(lambda x: str(x)+',')
        if  paralogs is not None:
            #print("Given annotation is mouse, position looking for is",option,"----")
            result["paralog"]      = result.apply(lambda row: paralogs(row,paralog_mouse_gene_list), axis=1)
            colnames                = ['Gene','paralog','chr','strand','start','end','CDS_end','CDS_start','transcript_symbol','transript_id']
            result.rename(columns   = {'CDS_end':'CDS_start','CDS_start':'CDS_end'},inplace=True)
            result                  = result[colnames]
            print("--------Successfully completed generating tsv file for",organism,option, "--------")
        else:
            result["paralog"]       = 0
            colnames                = ['Gene','paralog','chr','strand','start','end','CDS_end','CDS_start','transcript_symbol']
            result.rename(columns   = {'CDS_end':'CDS_start','CDS_start':'CDS_end'},inplace=True)
            result                  = result[colnames]
            print("--------Successfully completed generating tsv file for",organism,option, "--------")
    return(result)

if __name__ == "__main__":
    result     = main(gtffile,option,paralog,outputpath)
    result.to_csv(os.path.join(outputpath,organism+"_"+option+"_cannonical.tsv"),sep="\t",index=False)