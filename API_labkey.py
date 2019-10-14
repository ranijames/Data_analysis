#!/usr/bin/env python

from labkey.utils import create_server_context
from labkey.query import select_rows, QueryFilter
import pandas as pd
import re
import argparse
import os
import sys

"""
Usage : python API_labkey.py <output>
Author: Alva James, September 2019
"""
parser = argparse.ArgumentParser()
parser.add_argument("output", help="The path where you need your output labkeyprefix")

if len(sys.argv) < 1:
    parser.print_help()
    sys.exit(0)

args               = parser.parse_args()
output             =  args.output

# Configs to be added to ~/.netrc file as documented in `https://github.com/LabKey/labkey-api-python`

labkey_server = 'tp-labkey.ethz.ch:443'
project_name  = "Tumor Profiler - Melanoma"
contextPath   = 'labkey'

schema = 'study'
#table  = 'bkRNA_Sequencing'
table  = 'bkRNA_Analysis'
# 
# Create server context
server_context = create_server_context(labkey_server,
                                       project_name,
                                       contextPath,
                                       use_ssl=True)
# filter columns 
columns = "Name,FifthSequencingRunId"

# Load the entire dataset
data    = select_rows(server_context, schema, table,columns=columns)
main    = {}
data    = select_rows(server_context, schema, table)

# convert the dictionary format to the dataframe
main = data['rows']
df   = pd.DataFrame(main)
#df.rename(columns={'BkRnaLibrary':'labkey_prefix'},inplace=True)
print(df)
#df[["labkey_prefix","Name"]].to_csv(os.path.join(output+"/"+"labkeyprefix"),sep='\t',index=False)
