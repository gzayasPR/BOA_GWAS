#!/bin/sh
#SBATCH --job-name=Angus.HCW
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=160Gb
#SBATCH --time=72:00:00
date

results=/blue/mateescu/gzayas97/Gabe_Thesis/4.BO.GWAS.Heterosis/Analysis/results
cd ${results}
mkdir BLUPF90.SNP.GWAS
mkdir BLUPF90.BOGWAS
