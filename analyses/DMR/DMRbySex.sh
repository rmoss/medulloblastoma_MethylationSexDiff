#!/bin/bash -l
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=moss0134@umn.edu 
#SBATCH --job-name=DMRbysex
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G

set -euo pipefail

CavDir="/home/spectorl/moss0134/projects/research/MBSexDiff/data/Cavalli"
NewDir="/home/spectorl/moss0134/projects/research/MBSexDiff/data/Newcastle"
HeidDir="/home/spectorl/moss0134/projects/research/MBSexDiff/data/Hovestadt"

cd /home/spectorl/moss0134/projects/research/MBSexDiff/scripts


module load R/4.0.4

#Rscript --vanilla DMRbySex.R -g $CavDir/GRset.flt.rds -f $CavDir/phenoCav.csv -p Cavalli -n 500 -o $CavDir

Rscript --vanilla DMRbySex.R -g $NewDir/GRSet_Newcastle_GSE93646_filtMBonly.RDS -f $NewDir/Newcastle_pheno_v3.csv -p Newcastle -n 500 -o $NewDir
