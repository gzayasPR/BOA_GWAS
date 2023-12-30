#!/bin/bash
#SBATCH --account=mateescu
#SBATCH --job-name=BO.GWAS
#SBATCH --output=BO.GWAS_%j.out
#SBATCH --error=errorfile_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12Gb
#SBATCH --time=72:00:00
results=/blue/mateescu/gzayas97/Gabe_Thesis/4.BO.GWAS.Heterosis/Analysis/results
data=/blue/mateescu/gzayas97/Gabe_Thesis/4.BO.GWAS.Heterosis/Analysis/data
my_bin=/blue/mateescu/gzayas97/Gabe_Thesis/4.BO.GWAS.Heterosis/Analysis/bin
resultsBO=${results}/BO.BLUPF90/
resultsSNP=${results}/SNP.BLUPF90/
ALL_GWAS=${results}/ALL.GWAS/
newname=Seminole.Nutri
cd ${ALL_GWAS}
mkdir BO_SNP_GLM
cd BO_SNP_GLM
ml plink
for TRAIT in $(echo "MARB" "HCW" )
do
    GenoBO=${resultsBO}/Genotypes/${TRAIT}/Seminole.Nutri.BO
    GenoSNP=${resultsSNP}/Genotypes/${TRAIT}/Seminole.Nutri.SNP
    cd ${ALL_GWAS}/BO_SNP_GLM
    mkdir ${TRAIT}
    cd ${ALL_GWAS}/BO_SNP_GLM/${TRAIT}
    cp ${ALL_GWAS}/${TRAIT}/SIG.${TRAIT}.BO_SNP.txt .
    pheno=${resultsBO}/Genotypes/${TRAIT}/Cleaned.${TRAIT}.csv
    cp ${pheno} . 
    plink -cow -bfile ${GenoSNP} --recode A -extract SIG.${TRAIT}.BO_SNP.txt -out ${newname}.${TRAIT}.SNP
    plink -cow -bfile ${GenoBO} --ref-allele ${resultsBO}/Genotypes/${TRAIT}/ref.Brahman.txt -extract SIG.${TRAIT}.BO_SNP.txt --recode A -out ${newname}.${TRAIT}.BO
    done
