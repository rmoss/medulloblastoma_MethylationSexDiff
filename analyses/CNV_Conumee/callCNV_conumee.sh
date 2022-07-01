#!/bin/bash -l
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=moss0134@umn.edu 
#SBATCH --job-name=callCNV
#SBATCH --ntasks=1
#SBATCH --mem=100g
#SBATCH --tmp=100g

set -euo pipefail

CavDir="/home/spectorl/moss0134/projects/research/MBSexDiff/data/Cavalli"
NewDir="/home/spectorl/moss0134/projects/research/MBSexDiff/data/Newcastle"
HeidDir="/home/spectorl/moss0134/projects/research/MBSexDiff/data/Hovestadt"
controlDir = "/home/spectorl/moss0134/projects/research/MBSexDiff/data/control"
resultsDir = "/home/spectorl/moss0134/projects/research/MBSexDiff/results"

cd /home/spectorl/moss0134/projects/research/medulloblastoma_MethylationSexDiff/analyses/CNV_Conumee

module load R/4.0.4

Rscript --vanilla callCNV_conumee.R -a $CavDir/MSet_Cavalli_filt.RDS -b $CavDir/phenoCav.csv -c $controlDir/MSet_control_CBL_GSE134379.RDS -d $controlDir/Control_CBL_GSE134379_pheno_Sample_Sheet.csv -p Cavalli -s SHH -o $resultsDir

# Rscript --vanilla callCNV_conumee.R -a $CavDir/MSet_Cavalli_filt.RDS -b $CavDir/phenoCav.csv -c $controlDir/MSet_control_CBL_GSE134379.RDS -d $controlDir/Control_CBL_GSE134379_pheno_Sample_Sheet.csv -p Cavalli -s WNT -o $resultsDir

# Rscript --vanilla callCNV_conumee.R -a $CavDir/MSet_Cavalli_filt.RDS -b $CavDir/phenoCav.csv -c $controlDir/MSet_control_CBL_GSE134379.RDS -d $controlDir/Control_CBL_GSE134379_pheno_Sample_Sheet.csv -p Cavalli -s Group3 -o $resultsDir

# Rscript --vanilla callCNV_conumee.R -a $CavDir/MSet_Cavalli_filt.RDS -b $CavDir/phenoCav.csv -c $controlDir/MSet_control_CBL_GSE134379.RDS -d $controlDir/Control_CBL_GSE134379_pheno_Sample_Sheet.csv -p Cavalli -s Group4 -o $resultsDir

# Rscript --vanilla callCNV_conumee.R -a $NewDir/GRSet_Newcastle_GSE93646_filtMBonly.RDS -f $NewDir/Newcastle_pheno_v3.csv -p Newcastle -n 500 -o $NewDir
