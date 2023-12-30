#!/bin/bash
#SBATCH --account=mateescu
#SBATCH --job-name=AIREML.SNP
#SBATCH --output=AIREML.SNP_%j.out
#SBATCH --error=AIREML.SNP_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=7000mb
TRAIT=$1
results=/blue/mateescu/gzayas97/Gabe_Thesis/4.BO.GWAS.Heterosis/Analysis/results
data=/blue/mateescu/gzayas97/Gabe_Thesis/4.BO.GWAS.Heterosis/Analysis/data
my_bin=/blue/mateescu/gzayas97/Gabe_Thesis/4.BO.GWAS.Heterosis/Analysis/bin
resultsSNP=${results}/SNP.BLUPF90/
GenoSNP=${resultsSNP}/Genotypes/${TRAIT}/Seminole.Nutri.SNP
AIREML=${resultsSNP}/Genotypes/${TRAIT}/AIREML.parameter
pheno=${resultsSNP}/Genotypes/${TRAIT}/Cleaned.${TRAIT}.csv
cd ${resultsSNP}
mkdir ${TRAIT}
cd ${resultsSNP}/${TRAIT}
echo "Starting ${TRAIT}"
for MODEL in $(echo "add" "dom" "rec" "dom_dev" )
do
    echo "Starting ${MODEL} Model"
    cd ${resultsSNP}/${TRAIT}
    mkdir  SNP.GWAS.${MODEL}
    cd ${resultsSNP}/${TRAIT}/SNP.GWAS.${MODEL}
    mkdir AIREML 
    cd ${resultsSNP}/${TRAIT}/SNP.GWAS.${MODEL}/AIREML 
    cp ${GenoSNP}.geno.${MODEL} .
    cp ${GenoSNP}.PED.${MODEL} .
    cp ${GenoSNP}.map .
    sed -e 's/\t/ /g' -e 's/,/ /g' -e '1d' ${pheno} > pheno.txt
    ml blupf90
    echo "Starting AIREML for ${MODEL} Model"
    renumf90 ${AIREML}.${MODEL} > renum.out.txt
    airemlf90 renf90.par > ${MODEL}.airemlf90.txt
    grep "Final Estimates" -A 40 ${MODEL}.airemlf90.txt
    echo "AIREML ${MODEL}"
    remlf90 renf90.par > ${MODEL}.remlf90.txt
    echo "REML ${MODEL}"
    grep "Final Estimates" -A 40  ${MODEL}.remlf90.txt
done
