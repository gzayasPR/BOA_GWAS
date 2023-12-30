#!/bin/bash
#SBATCH --account=mateescu
#SBATCH --job-name=AIREML.BO
#SBATCH --output=AIREML.BO_%j.out
#SBATCH --error=AIREML.BO_%j.err
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
resultsBO=${results}/BO.BLUPF90/
GenoBO=${resultsBO}/Genotypes/${TRAIT}/Seminole.Nutri.BO 
AIREML=${resultsBO}/Genotypes/${TRAIT}/AIREML.parameter
pheno=${resultsBO}/Genotypes/${TRAIT}/Cleaned.${TRAIT}.csv
cd ${resultsBO}
mkdir ${TRAIT}
cd ${resultsBO}/${TRAIT}
echo "AIREML FOR ${TRAIT}"
for MODEL in $(echo "add" "dom" "rec" "dom.dev" )
do
    echo "Starting ${MODEL} Model"
    cd ${resultsBO}/${TRAIT}
    mkdir  BO.GWAS.${MODEL}
    cd ${resultsBO}/${TRAIT}/BO.GWAS.${MODEL}
    mkdir AIREML 
    cd ${resultsBO}/${TRAIT}/BO.GWAS.${MODEL}/AIREML 
    cp ${GenoBO}.geno.${MODEL} .
    cp ${GenoBO}.PED.${MODEL} .
    cp ${GenoBO}.map .
    sed -e 's/\t/ /g' -e 's/,/ /g' -e '1d' ${pheno} > pheno.txt
    ml blupf90
    echo "Starting AIREML for ${MODEL} Model"
    renumf90 ${AIREML}.${MODEL} > renum.out.txt
    airemlf90 renf90.par > ${MODEL}.airemlf90.txt
    grep "Final Estimates" -A 40  ${MODEL}.airemlf90.txt
    echo "AIREML ${MODEL}"
    remlf90 renf90.par > ${MODEL}.remlf90.txt
    echo "REML ${MODEL}"
    grep "Final Estimates" -A 40  ${MODEL}.remlf90.txt
done
