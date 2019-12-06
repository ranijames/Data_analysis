#!/bin/bash

set -e

dataset_align=20180802
dataset=20180802_T1012_splDev_ippDev

basedir=/cluster/work/grlab/projects/coukos_immuno
annotation=${basedir}/annotation/gencode.v19.annotation.hs37d5_chr.gtf ### was used for version 1
#annotation=${basedir}/annotation/gencode.v28lift37.annotation.gtf 
outdir=${basedir}/results_${dataset}/expression
mkdir -p $outdir

for file in ${basedir}/alignment_${dataset_align}/*.bam
do
    filebase=$(basename $file)
    cntfile=${outdir}/${filebase%bam}tsv
    logfile=${outdir}/${filebase%bam}log
    if [ ! -f ${cntfile} ]
    then
        echo "cd $(pwd); python count_expression.py -m -B -v -A $file -a $annotation -o ${cntfile} > ${logfile} 2>&1" | bsub -J cntPHRT -n 1 -W 24:00 -R "rusage[mem=15000]" -o /dev/null 
    fi
done
