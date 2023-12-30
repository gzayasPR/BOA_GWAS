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

## Bash code to submit multiple single trait single-step GWAS analysis using BLUPf90 ##
#First round not weighted ssGWAS
#Second and Third rounds weighted ssGWAS (updating G matrix)
#Saving output into different folders (by trait)
#FMR
#Oct 28, 2020
TRAIT=$1
results=/blue/mateescu/gzayas97/Gabe_Thesis/4.BO.GWAS.Heterosis/Analysis/results
data=/blue/mateescu/gzayas97/Gabe_Thesis/4.BO.GWAS.Heterosis/Analysis/data
my_bin=/blue/mateescu/gzayas97/Gabe_Thesis/4.BO.GWAS.Heterosis/Analysis/bin
SIZE_KB=$2
windows_variance_mbp=$3
resultsSNP=${results}/SNP.BLUPF90/
GenoSNP=${resultsSNP}/Genotypes/${TRAIT}/Seminole.Nutri.SNP
GBLUP=${resultsSNP}/Genotypes/${TRAIT}/GBLUP.parameter
pheno=${resultsSNP}/Genotypes/${TRAIT}/Cleaned.${TRAIT}.csv
ALL_GWAS=${results}/ALL.GWAS/
WEIGHTS=$4
ml blupf90
ml R
echo "Starting ${TRAIT}"
cd ${resultsSNP}/${TRAIT}
for MODEL in $(echo "add" "dom" "rec" "dom_dev")
do
    echo "Starting ${MODEL} Model"
    cd ${resultsSNP}/${TRAIT}/SNP.GWAS.${MODEL}
    mkdir GBLUP_${SIZE_KB}
    cd ${resultsSNP}/${TRAIT}/SNP.GWAS.${MODEL}/GBLUP_${SIZE_KB}
    cp ${GenoSNP}.geno.${MODEL} .
    cp ${GenoSNP}.PED.${MODEL} .
    cp ${GenoSNP}.map .
    sed -e 's/\t/ /g' -e 's/,/ /g' -e '1d' ${pheno} > pheno.txt
  weight_default="NONLINEAR"
  if [[ "$WEIGHTS" == "$weight_default" ]]
  then
    round=(1 2 3)	
  else
    round=(1 2)	
  fi
  for i in "${round[@]}"
  do
    cd ${resultsSNP}/${TRAIT}/SNP.GWAS.${MODEL}/GBLUP_${SIZE_KB}/
    mkdir -p "round_"$i                #create one folder per round
    cd "round_"$i
    echo "Start round_"$i
    cp ${GenoSNP}.geno.${MODEL} .
    cp ${GenoSNP}.PED.${MODEL} .
    cp ${GenoSNP}.map .
    sed -e 's/\t/ /g' -e 's/,/ /g' -e '1d' ${pheno} > pheno.txt 
      if (($i == 1))                                  #saving weights to calculate G matrix by blupf90
      then
        awk 'NR>1 {print "1",$1}' ${GenoSNP}.map | awk '{print $1}' > wei.txt    #weights=1 for all SNPs for round_1
      else
        awk 'NR>1 {print $7}' ../round_$((i - 1))/snp_sol > wei.txt           #saving weights from snp_sol of previous round for round_2 and round_3
      fi
    cp ${GBLUP}.${MODEL} .
    dos2unix ${GBLUP}.${MODEL}
    sed -i GBLUP.parameter.${MODEL} -e '$aOPTION weightedG wei.txt'
    printf GBLUP.parameter.${MODEL} | renumf90 > renum.${i}.${MODEL}.txt
    echo "renf90.par" | blupf90 > GBLUP.${i}.${MODEL}.txt
    
      weight_default="NONLINEAR"
      if [[ "$WEIGHTS" == "$weight_default" ]]                              #saving weights to calculate G matrix by blupf90
      then
        sed "s/saveGInverse/readGInverse/g" renf90.par | sed "s/saveA22Inverse/readA22Inverse/g" | sed "s/OPTION missing -999/#OPTION missing -999/g" | sed "s/solution mean.*/windows_variance_mbp $windows_variance_mbp/g" | sed "s/sol se/which_weight nonlinearA 1.05/g" | sed -e '$aOPTION Manhattan_plot_R'  > POSTGS.Parameter        #to create plots add: sed -e '$aOPTION Manhattan_plot_R'    #weights=1 for all SNPs for round_1
      else
        sed "s/saveGInverse/readGInverse/g" renf90.par | sed "s/saveA22Inverse/readA22Inverse/g" | sed "s/OPTION missing -999/#OPTION missing -999/g" | sed "s/solution mean.*/windows_variance_mbp $windows_variance_mbp/g" | sed "s/OPTION sol se/#OPTION sol se/g"  | sed -e '$aOPTION Manhattan_plot_R'  > POSTGS.Parameter        #to create plots add: sed -e '$aOPTION Manhattan_plot_R'
      fi
    sed "s/saveGInverse/readGInverse/g" renf90.par | sed "s/saveA22Inverse/readA22Inverse/g" | sed "s/OPTION solution mean/#OPTION solution mean/g" | sed "s/OPTION missing -999/#OPTION missing -999/g"| sed "s/OPTION sol se/#OPTION sol se/g"  | sed -e '$aOPTION saveHinv' > PREGS.Parameter 
    printf PREGS.Parameter | preGSf90 > preGSf90.${i}.${MODEL}.txt
    printf POSTGS.Parameter | postGSf90 > postGSf90.${i}.${MODEL}.txt
    rm ${name}.geno.${MODEL}
    Rscript Sft*.R
    Rscript Vft*.R
    Rscript Pft*.R
    echo "Finished round_"$i   
    cd ${resultsSNP}/${TRAIT}/SNP.GWAS.${MODEL}/GBLUP_${SIZE_KB}
  done
  cd ${ALL_GWAS}
  mkdir ${TRAIT}
  cd ${ALL_GWAS}/${TRAIT}
  mkdir SNP
  cd ${ALL_GWAS}/${TRAIT}/SNP
  weight_default="NONLINEAR"
  if [[ "$WEIGHTS" == "$weight_default" ]]                              #saving weights to calculate G matrix by blupf90
    then
    cp ${resultsSNP}/${TRAIT}/SNP.GWAS.${MODEL}/GBLUP_${SIZE_KB}/round_3/Vft*_manplot.pdf VSNP.${SIZE_KB}.${TRAIT}.${MODEL}.${WEIGHTS}.pdf
    cp ${resultsSNP}/${TRAIT}/SNP.GWAS.${MODEL}/GBLUP_${SIZE_KB}/round_1/Pft*_manplot.pdf PSNP.${TRAIT}.${MODEL}.${WEIGHTS}.pdf
    cp ${resultsSNP}/${TRAIT}/SNP.GWAS.${MODEL}/GBLUP_${SIZE_KB}/round_3/chrsnp SNP.Effects.${TRAIT}.${MODEL}.${WEIGHTS}.txt
    cp ${resultsSNP}/${TRAIT}/SNP.GWAS.${MODEL}/GBLUP_${SIZE_KB}/round_1/chrsnp_pval SNP.p_val.${TRAIT}.${MODEL}.${WEIGHTS}.txt
    cp ${resultsSNP}/${TRAIT}/SNP.GWAS.${MODEL}/GBLUP_${SIZE_KB}/round_3/chrsnpvar SNP.var.${SIZE_KB}.${TRAIT}.${MODEL}.${WEIGHTS}.txt
    cp $GenoSNP.map .
    echo "Finished ${MODEL} Model"
  else
    cp ${resultsSNP}/${TRAIT}/SNP.GWAS.${MODEL}/GBLUP_${SIZE_KB}/round_1/Pft*_manplot.pdf PSNP.${TRAIT}.${MODEL}.${WEIGHTS}.pdf
    cp ${resultsSNP}/${TRAIT}/SNP.GWAS.${MODEL}/GBLUP_${SIZE_KB}/round_2/chrsnp SNP.Effects.${TRAIT}.${MODEL}.${WEIGHTS}.txt
    cp ${resultsSNP}/${TRAIT}/SNP.GWAS.${MODEL}/GBLUP_${SIZE_KB}/round_1/chrsnp_pval SNP.p_val.${TRAIT}.${MODEL}.${WEIGHTS}.txt
    cp ${resultsSNP}/${TRAIT}/SNP.GWAS.${MODEL}/GBLUP_${SIZE_KB}/round_2/chrsnpvar SNP.var.${SIZE_KB}.${TRAIT}.${MODEL}.${WEIGHTS}.txt
    cp $GenoSNP.map .
    echo "Finished ${MODEL} Model"
  fi
done
