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
resultsBO=${results}/BO.BLUPF90/
GenoBO=${resultsBO}/Genotypes/${TRAIT}/Seminole.Nutri.BO
GBLUP=${resultsBO}/Genotypes/${TRAIT}/GBLUP.parameter
pheno=${resultsBO}/Genotypes/${TRAIT}/Cleaned.${TRAIT}.csv
ALL_GWAS=${results}/ALL.GWAS/
WEIGHTS=$4
ml blupf90
ml R
echo "Starting ${TRAIT}"
cd ${resultsBO}
cd ${resultsBO}/${TRAIT}
for MODEL in $(echo "add" "dom" "rec" "dom.dev")
do
    echo "Starting ${MODEL} Model"
    cd ${resultsBO}/${TRAIT}/BO.GWAS.${MODEL}
    mkdir GBLUP_${SIZE_KB}
    cd ${resultsBO}/${TRAIT}/BO.GWAS.${MODEL}/GBLUP_${SIZE_KB}
    cp ${GenoBO}.geno.${MODEL} .
    cp ${GenoBO}.PED.${MODEL} .
    cp ${GenoBO}.map .
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
    cd ${resultsBO}/${TRAIT}/BO.GWAS.${MODEL}/GBLUP_${SIZE_KB}/
    mkdir -p "round_"$i                #create one folder per round
    cd "round_"$i
    echo "Start round_"$i
    cp ${GenoBO}.geno.${MODEL} .
    cp ${GenoBO}.PED.${MODEL} .
    cp ${GenoBO}.map .
    sed -e 's/\t/ /g' -e 's/,/ /g' -e '1d' ${pheno} > pheno.txt 
      if (($i == 1))                                  #saving weights to calculate G matrix by blupf90
      then
        awk 'NR>1 {print "1",$1}' ${GenoBO}.map | awk '{print $1}' > wei.txt    #weights=1 for all SNPs for round_1
      else
        awk 'NR>1 {print $7}' ../round_$((i - 1))/snp_sol > wei.txt           #saving weights from snp_sol of previous round for round_2 and round_3
      fi
    cp ${GBLUP}.${MODEL} .
    dos2unix ${GBLUP}.${MODEL}
    sed -i GBLUP.parameter.${MODEL} -e '$aOPTION weightedG wei.txt'
    printf GBLUP.parameter.${MODEL} | renumf90 > renum.${i}.${MODEL}.txt
    echo "renf90.par" | blupf90 > GBLUP.${i}.${MODEL}.txt
    
      weight_default="NONLINEAR"
      if [[ "$WEIGHTS" == "$weight_default" ]]                                   #saving weights to calculate G matrix by blupf90
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
    cd ${resultsBO}/${TRAIT}/BO.GWAS.${MODEL}/GBLUP_${SIZE_KB}
  done
  cd ${ALL_GWAS}
  mkdir ${TRAIT}
  cd ${ALL_GWAS}/${TRAIT}
  mkdir BO
  weight_default="NONLINEAR"
  if [[ "$WEIGHTS" == "$weight_default" ]]
    then
    cp ${resultsBO}/${TRAIT}/BO.GWAS.${MODEL}/GBLUP_${SIZE_KB}/round_3/Vft*_manplot.pdf VBO.${SIZE_KB}.${TRAIT}.${MODEL}.${WEIGHTS}.pdf
    cp ${resultsBO}/${TRAIT}/BO.GWAS.${MODEL}/GBLUP_${SIZE_KB}/round_1/Pft*_manplot.pdf PBO.${TRAIT}.${MODEL}.${WEIGHTS}.pdf
    cp ${resultsBO}/${TRAIT}/BO.GWAS.${MODEL}/GBLUP_${SIZE_KB}/round_3/chrsnp BO.Effects.${TRAIT}.${MODEL}.${WEIGHTS}.txt
    cp ${resultsBO}/${TRAIT}/BO.GWAS.${MODEL}/GBLUP_${SIZE_KB}/round_1/chrsnp_pval BO.p_val.${TRAIT}.${MODEL}.${WEIGHTS}.txt
    cp ${resultsBO}/${TRAIT}/BO.GWAS.${MODEL}/GBLUP_${SIZE_KB}/round_3/chrsnpvar BO.var.${SIZE_KB}.${TRAIT}.${MODEL}.${WEIGHTS}.txt
    cp $GenoBO.map .
    echo "Finished ${MODEL} Model"
    else
    cp ${resultsBO}/${TRAIT}/BO.GWAS.${MODEL}/GBLUP_${SIZE_KB}/round_2/Vft*_manplot.pdf VBO.${SIZE_KB}.${TRAIT}.${MODEL}.${WEIGHTS}.pdf
    cp ${resultsBO}/${TRAIT}/BO.GWAS.${MODEL}/GBLUP_${SIZE_KB}/round_1/Pft*_manplot.pdf PBO.${TRAIT}.${MODEL}.${WEIGHTS}.pdf
    cp ${resultsBO}/${TRAIT}/BO.GWAS.${MODEL}/GBLUP_${SIZE_KB}/round_2/chrsnp BO.Effects.${TRAIT}.${MODEL}.${WEIGHTS}.txt
    cp ${resultsBO}/${TRAIT}/BO.GWAS.${MODEL}/GBLUP_${SIZE_KB}/round_1/chrsnp_pval BO.p_val.${TRAIT}.${MODEL}.${WEIGHTS}.txt
    cp ${resultsBO}/${TRAIT}/BO.GWAS.${MODEL}/GBLUP_${SIZE_KB}/round_2/chrsnpvar BO.var.${SIZE_KB}.${TRAIT}.${MODEL}.${WEIGHTS}.txt
    cp $GenoBO.map .
    echo "Finished ${MODEL} Model"
    fi
done
