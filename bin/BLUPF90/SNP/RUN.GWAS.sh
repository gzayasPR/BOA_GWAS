#!/bin/bash
#SBATCH --account=mateescu
#SBATCH --job-name=SNP.GWAS
#SBATCH --output=SNP.GWAS_%j.out
#SBATCH --error=SNP.GWASS.errorfile_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8Gb
#SBATCH --time=72:00:00

my_bin=/blue/mateescu/gzayas97/Gabe_Thesis/4.BO.GWAS.Heterosis/Analysis/bin/BLUPF90/SNP/
cd ${my_bin}
for TRAIT in $(echo "MARB" "HCW")
do
    echo "Running Trait ${TRAIT}"
    sbatch BLUPf90.ssGWAS.all.sh ${TRAIT} 10 0.01 "NONLINEAR"
    sbatch BLUPf90.ssGWAS.all.sh ${TRAIT} 200 0.2 "NONLINEAR"
    sbatch BLUPf90.ssGWAS.all.sh ${TRAIT} 1000 1 "NONLINEAR"
    done