#!/bin/bash
#SBATCH --account=mateescu
#SBATCH --job-name=Test.it
#SBATCH --output=Test.it_%j.out
#SBATCH --error=errorfile_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=7000mb

results=/blue/mateescu/gzayas97/Gabe_Thesis/3.BOA.GWAS/Analysis/results/BO.BLUPF90/

cd ${results}
mkdir ./Geno.Add
cd ./Geno.Add
cp ${results}/Genotypes/SNP/Seminole.Nutri.PED.add .
#cp ${results}/Genotypes/SNP/HCW.csv  .
#cp ${results}/AIREML.parameter .
#cp ${results}/GBLUP.parameter .
cp ${results}/Genotypes/SNP/Seminole.Nutri.map .
cp ${results}/Genotypes/SNP/Seminole.Nutri.geno.add .

sed -e 's/\t/ /g' -e 's/,/ /g' -e '1d' HCW.csv > pheno.txt

#Edit gblup and aireml files
head HCW.csv
ml blupf90
renumf90 AIREML.parameter
airemlf90 renf90.par
#H2 = 0.31
remlf90 renf90.par
#A2 = 1591
#r2 = 35433543
renumf90 GBLUP.parameter
 blupf90test renf90.par

sed "s/saveGInverse/readGInverse/g" renf90.par | sed "s/saveA22Inverse/readA22Inverse/g" | sed "s/OPTION missing -999/#OPTION missing -999/g" | 
  sed -e '$aOPTION Manhattan_plot_R'  > POSTGS.Parameter


#sed "s/saveGInverse/readGInverse/g" renf90.par | sed "s/saveA22Inverse/readA22Inverse/g" | sed "s/OPTION missing -999/#OPTION missing -999/g" | 
#  sed "s/solution mean.*/windows_variance_mbp 0.01/g" | sed "s/sol se/which_weight nonlinearA 1.05/g" | sed -e '$aOPTION SNP_variance_limit 1.4775' | 
#  sed -e '$aOPTION Manhattan_plot_R'| sed "s/OPTION solv_method FSPAK/#OPTION solv_method FSPAK/g" |
#  sed "s/OPTION omit_ainv/#OPTION omit_ainv/g" | sed "s/OPTION TauOmega 1.0 0.0/#OPTION TauOmega 1.0 0.0/g" |
#  sed "s/OPTION AlphaBeta 0.95 0.05/#OPTION AlphaBeta 0.95 0.05/g"    > POSTGS.Parameter        #to create plots add: sed -e '$aOPTION Manhattan_plot_R'

ml R 
 printf POSTGS.Parameter | postGSf90
 
library(GWASTools)
pval <- read.table("chrsnp_pval")[,c(4,5,6,3)]
names(pval) <- c("SNP","CHR","BP","-log(p-value)")
pvalue <- 10^(pval[,4]*(-1))
qqPlot(pvalue)
