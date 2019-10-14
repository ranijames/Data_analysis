import numpy as np
import argparse 
import pandas as pd
import sys
import os

"""
Purpose:    This program parses the GENCODE annoation file
Author:     Alva James
Date:       September 19, 2019
Usage:      python annotation_parser.py <-p paralog/optional> <gtf> <option> <output> <organism>
example for mouse: python annotation_parser.py --paralog mouse_paralogs.GRCm38.p6.biomart.tsv Mouse_vM22_annotation_altered2.gtf CDS /Users/alvajames/Personal_offcial/Olga_task_Sep mouse
"""

parser = argparse.ArgumentParser()
parser.add_argument('Mouse_annotation',help='The converted mouse annotation file: gene_transcript_mouse_annotation_part1.txt')
parser.add_argument('metafeature', help='The metafeature file downloaded from GENCODE: gencode.vM22.metadata.Transcript_supporting_feature')
parser.add_argument('outputpath',  help='Path where you need to write the output file')


args = parser.parse_args()
 
if len(sys.argv) < 3:
    parser.print_help()
    sys.exit(0)
    
output                  = args.outputpath    
gene_transcript         = pd.read_csv(args.Mouse_annotation,sep='\t',header=0,skipinitialspace=True)
gene_transcript.columns = gene_transcript.columns.str.strip()
meta_feature            = pd.read_csv(args.metafeature,sep='\t',skipinitialspace=True)
meta_feature.columns    = meta_feature.columns.str.strip()

# Merge the mouse transcripts and the meta feature file on the transcript ids
mouse_meta              = pd.merge(gene_transcript, meta_feature[["transcript_id","evidence"]],
                         on="transcript_id")
# Generate the comma seperated list of evidence based on the evidence column
df_mouse_meta = (mouse_meta['evidence'].str.split(', ')
                    .groupby(mouse_meta['transcript_id'])
                    .agg(lambda x: ', '.join(set(y for z in x for y in z)))
                    .reset_index())

df_final= pd.merge(df_mouse_meta,gene_transcript, on="transcript_id")
df_final[["gene_id","transcript_id","evidence"]].to_csv(os.path.join(output+"/"+"mouse_metainfo_gene_transcript_final.txt"),sep='\t',index=False)