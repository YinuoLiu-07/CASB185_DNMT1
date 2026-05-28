#!/bin/bash
#$ -V
#$ -l h_data=8G,h_rt=1:00:00
#$ -cwd
#$ -m a
#$ -o motif.out
#$ -e motif.error

### USAGE: qsub find_motif_CCG.sh

### LOAD MODULES ###
. /u/local/Modules/default/init/modules.sh
module load homer/4.11.1

### SET PATHS ###
work_dir=/u/project/gxxiao/yinuoliu/CASB185
bed_dir=${work_dir}/eCLIP_peak/GSE290473_DNMT1_reproducible_peaks_in_293T.bed
unique_bed_dir=${work_dir}/eCLIP_peak/GSE290473_DNMT1_unique_peaks.bed
genome="hg38"
out_dir=${work_dir}/find_motif
pre_parse=${work_dir}/homer_preparse

### LABEL PEAKS ###
# This creates a new BED file where column 4 is given a unique name: peak_1, peak_2, etc.
awk 'BEGIN{OFS="\t"} {$4="peak_"NR; print}' ${bed_dir} > ${unique_bed_dir}

### RUN HOMER ###
findMotifsGenome.pl ${unique_bed_dir} ${genome} ${out_dir} -size 200 -preparsedDir ${pre_parse}

