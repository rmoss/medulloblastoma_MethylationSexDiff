#!/bin/bash -l
#SBATCH --time=6:00:00
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=moss0134@umn.edu 
#SBATCH --job-name=DMPbysex
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G

set -euo pipefail

CavDir="/home/spectorl/moss0134/projects/research/MBSexDiff/data/Cavalli"
NewDir="/home/spectorl/moss0134/projects/research/MBSexDiff/data/Newcastle"
HeidDir="/home/spectorl/moss0134/projects/research/MBSexDiff/data/Hovestadt"

cd /home/spectorl/moss0134/projects/research/medulloblastoma_MethylationSexDiff/analyses/DMP


module load R/4.0.4

echo "starting new DMP script"
Rscript --vanilla subtype_DMPbySex.R -g $CavDir/GRset.flt.rds -f $CavDir/pheno_CavalliCombined_subgroupandtype_filtered.csv  -p Cavalli_SHH -n 500 -o $CavDir/subtype_DMP

#Rscript --vanilla DMPbySex.R -g $NewDir/GRSet_Newcastle_GSE93646_filtMBonly.RDS -f $NewDir/Newcastle_pheno_v3.csv -p Newcastle -n 500 -o $NewDir/subtype_DMP
