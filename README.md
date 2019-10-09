# Data Analysis


The python snippets for manipulation of big data datasets of multiple data types, from human diseases

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

7. `dataframe_multi_transformation.py`
 Parse through multiple files of different diamensions, select columns if it is equal to the filename and concatenate them into a single dataframe. In addition to that, make long formatted dataframe based on commonly occuring columns names. Striping multiple characters including, square brackets, white space, digits etc.
 
 8. `list_to_dataframe.py`
The script takes in list of fastq files from mutiple lanes and both reads. Then converts the input into a dataframe with 4 columns. Column with `sample_id`, `lane` , `fastq1` and `fastq2` reads.

9. `GT_converter.py`
The script to convert the columns from a GTF file to a tsv file. The contents used for convertion are the exon or CDS positions for each transcripts of a gene. 

 10. `annotation_parser.py`
 Script collects and make a lists of positions for each exon and CDS for all trancripts in a altered GTF file. This conversions can be done on both mouse and human annoation file. 
 
  11. `dropbox_symlink_maker.py`
  Scripts generates symlinks for files from a list of diretcories, for given list of file name extensions.
  
  12. `fastq_mover.py`
  Scriot moves the fastq files from source to desination based on conditions. Additionally, the script makes destination directory based on the given list of sample_ids. The sample_id's are extarcted from the list of fastq file names. The input is path to fastq files and path where you need to generate the desination directories and finally move the `fastq` files
  
  
