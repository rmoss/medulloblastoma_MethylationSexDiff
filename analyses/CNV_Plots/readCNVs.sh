#!/bin/bash -l

set -euo pipefail


cd /home/spectorl/moss0134/projects/research/medulloblastoma_MethylationSexDiff/analyses/CNV_Plots

module load R/4.0.4
module load bedtools

#Rscript --vanilla readCNVs.R 

cd /home/spectorl/moss0134/projects/research/MBSexDiff/results/tables/
d=`ls *_CNV.bed`

for f in $d 
do
    name=${f/.bed}
    sort -k1,1 -k2,2n $f > $name.sorted.bed
    name2=${name/_CNV}
    # sort -k1,1 -k2,2n ${name2}_dmp.bed > $name2.dmp.sorted.bed
    # bedtools intersect -a $name2.dmp.sorted.bed -b $name.sorted.bed -sorted -wao > $name2.CNVfeature_intersect.txt
    # bedtools intersect -a $name2.dmp.sorted.bed -b $name.sorted.bed -sorted -C > $name2.CNVfeature_intersect_count.txt
    sort -k1,1 -k2,2n genearm_loc.bed > genearm_loc.sorted.bed
    bedtools intersect -a genearm_loc.sorted.bed -b $name.sorted.bed -sorted -wao > $name2.CNVfeature_intersect_genearm.txt
    bedtools intersect -a genearm_loc.sorted.bed -b $name.sorted.bed -sorted -C > $name2.CNVfeature_intersect_genearm_count.txt  
done