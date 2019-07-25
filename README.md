# Data_analysis


The python snippets for manipulation of big data datasets of multiple data types

The repository consists of codes/python snippets useful for big data manipulations, including converting the files, decoding, 
summing up across the columns, visualizing, etc.

1. `artificial_bulk.ipynb`
This notebook consits of python code to transform a single cell sequencing dataset in `hdf5` format into a bulk RNA-seq expression matrix (`tsv` file)

2. 	`merging_multiple_files_intone.ipynb`
This notebook shows how to merge multiple tsv files into a single dataframe

3. `data_report_visual_analysis.ipynb`
 This notebook represents the visual analysis of a given dataset as boxplots, category plot , and swarmplots using python seaborn library
 
 4. `list_difference_remove_unwanted_characters.ipynb`
 This notebook contains scripts for comparing two lists, for there difference. In addition to that extarcting rows within a dataframe based on which column you need, and removing unwanted characters (`.`) from the gene symbols.
 
 5. `merge_dictionary_values.ipynb`
 This notebook contains scripts for conditional merging of values in a dictionary (based on their keys)

6. `datframe_maker.py`
Helps to concatenate different matrices of different shape into a single dataframe based on conitions, and add additional columns to it fill the empty spaces NA . 
In addition to that, the script helps to aggregate a dataframe based on the frequency of which the column and rows occur together, aka "cross tabulation" using pandas `pivot_table`
