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
date
data=/blue/mateescu/gzayas97/Gabe_Thesis/3.BOA.GWAS/Data
results=/blue/mateescu/gzayas97/Gabe_Thesis/3.BOA.GWAS/Analysis/results

cd ${results}
mkdir BO.GWAS
cd ${results}/BO.GWAS
mkdir Genotypes
cd ${results}/BO.GWAS/Genotypes
mkdir BO
cd ${results}/BO.GWAS/Genotypes/BO
ml plink/1.90b3.39

plink -cow -bfile ${data}/SeminNutri -maf 0.05 -geno 0.1 -mind 0.1 --recode12 --output-missing-genotype 0 --transpose --out SeminNutri


#IBS matrix 
emmax -kin -intel64 -v -s -d 10 SeminNutri

BN (Balding-Nichols) matrix ? % emmax-kin-intel64 -v -d 10 [tped_prefix] (will generate [tped_prefix].aBN.kinf)

awk -F "\t" '{print $1,$2,$6}' FS='\t' $data/SeminNutri.BO.fam > pheno.txt
plink -cow -bfile ${data}/SeminNutri.BO  --linear hide-covar dominant --covar ${data}/mycov.txt --covar-number 1,5 --pheno  pheno.txt --allow-no-sex -out Dom.BO

plink -cow -bfile ${data}/SeminNutri.BO  --linear hide-covar recessive --covar ${data}/mycov.txt --covar-number 1,5 --pheno  pheno.txt --allow-no-sex -out Rec.BO
plink -cow -bfile ${data}/SeminNutri.BO  --linear hide-covar --covar ${data}/mycov.txt --covar-number 1,5 --pheno  pheno.txt --allow-no-sex -out ADD.BO
library(CMplot)
dom <- read.table("Dom.BO.assoc.linear",header=T)
dom.p <- dom[,c(2,1,3,9)]
CMplot(dom.p)

rec <- read.table("Rec.BO.assoc.linear",header=T)
rec.p <- rec[,c(2,1,3,9)]
CMplot(rec.p)


ADD <- read.table("ADD.BO.assoc.linear",header=T)
ADD.p <- ADD[,c(2,1,3,9)]

unique.p <- ADD.p[!duplicated(ADD.p[,c('CHR','P')]),]
write.table(unique.p[1],"Unique.BO.txt",quote=F,row.names=F)
CMplot(dom.p)
