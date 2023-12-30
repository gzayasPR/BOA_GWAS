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
results=/blue/mateescu/gzayas97/Gabe_Thesis/4.BO.GWAS.Heterosis/Analysis/results
data=/blue/mateescu/gzayas97/Gabe_Thesis/4.BO.GWAS.Heterosis/Analysis/data
my_bin=/blue/mateescu/gzayas97/Gabe_Thesis/4.BO.GWAS.Heterosis/Analysis/bin
BOGeno=${data}/Genotypes/BOA/Seminole.Nutritional.BO
newname=Seminole.Nutri.BO
phenotxt=${data}/Phenotypes/
AIREML=${my_bin}/BLUPF90/AIREML.parameter
GBLUP=${my_bin}/BLUPF90/GBLUP.parameter
cd ${results}
mkdir BO.BLUPF90
cd ${results}/BO.BLUPF90
mkdir Genotypes
cd ${results}/BO.BLUPF90/Genotypes
ml plink
for TRAIT in $(echo)
do
    cd ${results}/BO.BLUPF90/Genotypes
    mkdir ${TRAIT}
    cd ${results}/BO.BLUPF90/Genotypes/${TRAIT}
    cp ${phenotxt}/Cleaned.${TRAIT}.csv .
    cp ${AIREML} .
    TRAITS_HEAD=$(echo "#`head -1 Cleaned.${TRAIT}.csv`")
    sed -e "s/#YOUR_TRAITS/${TRAITS_HEAD}/g" \
    -e "s/YOUR_PED/${newname}.PED.add/g" \
    -e "s/YOUR_GENO/${newname}.geno.add/g" \
    -e "s/YOUR_MAP/${newname}.map/g" \
    AIREML.parameter > AIREML.parameter.add
    sed -e "s/.add/.dom/g" AIREML.parameter.add > AIREML.parameter.dom
    sed -e "s/.add/.rec/g" AIREML.parameter.add > AIREML.parameter.rec
    sed -e "s/.add/.dom.dev/g" AIREML.parameter.add > AIREML.parameter.dom.dev
    rm AIREML.parameter
    cp ${GBLUP} .
    sed -e "s/#YOUR_TRAITS/${TRAITS_HEAD}/g" \
    -e "s/YOUR_PED/${newname}.PED.add/g" \
    -e "s/YOUR_GENO/${newname}.geno.add/g" \
    -e "s/YOUR_MAP/${newname}.map/g" \
    GBLUP.parameter > GBLUP.parameter.add
    sed -e "s/.add/.dom/g" GBLUP.parameter.add > GBLUP.parameter.dom
    sed -e "s/.add/.rec/g" GBLUP.parameter.add > GBLUP.parameter.rec
    sed -e "s/.add/.dom.dev/g" GBLUP.parameter.add > GBLUP.parameter.dom.dev
    rm GBLUP.parameter
    cd ${results}/BO.BLUPF90/Genotypes
done

# Edit ${TRAIT}.csv to only get individuals with ${TRAIT}
#HCW.csv removed uncessary columns, removed missing hcw records, and records with missing year,cg,pca etc. 

for TRAIT in $(echo "MARB" "HCW" )
do
    cd ${results}/BO.BLUPF90/Genotypes/${TRAIT}
    #HCW.csv removed uncessary columns, removed missing hcw records, and records with missing year,cg,pca etc. 
    awk -F "," '{print $1}' Cleaned.${TRAIT}.csv > id.temp
    cp ${BOGeno}.fam FAM.file
    grep -wrf id.temp FAM.file | awk '{print $1,$2}' > id2.temp
    # keep 
    plink -cow -bfile ${BOGeno} -chr 1-29  -keep id2.temp -maf 0.01 -geno 0.1 -mind 0.1 -make-bed --out data.temp
    awk -F " " '{print $2, "A"}' data.temp.bim > ref.Brahman.txt
    plink -cow -bfile data.temp --ref-allele ref.Brahman.txt --recode A  -make-bed --out ${newname}
    plink -cow -bfile ${newname}  --freq --out  ${newname}
    plink -cow -bfile ${newname} --geno-counts --out   ${newname}
    #awk '{print $2,0,0}' ${newname}.fam > ped.txt
    cut -c1-150 ${newname}.raw |head -11

    #BLUPF90 needs
    ## ID genotypes_no_Spaces(aa=0,Aa=1,AA=2,Missing=5)
    #PLINK.raw
    ## FID IID PAT MAT SEX PHENOTYPE genotypes_with_spaces(aa=0,Aa=1,AA=2,Missing=NA)
    ## All with headers
    #First we need to remove the header
    ## Need to take out column 2 (IDs) and all genotypes
    ## In the genotypes we need to replace NA with 5 and remove spaces between genotypes
    #Take Column two out and make new file with IDs while removing first row
    awk '{print $2}'   ${newname}.raw | sed '1d' > Ids.temp
    awk '{ print length }' Ids.temp |sort | uniq
    awk -v m=6 '{printf("%-6s\n", $0)}' Ids.temp  > Ids2.temp
    #Cut first 6 columns , sed -e "removes first colum" -e "replaces NA with5" -e "Remove Blanks spaces between genotypes" 
    # This part takes some time, if you ever want to check what animal you are on do wc -l geno.temp
    cut -f7- ${newname}.raw | sed -e '1d' -e 's/NA/5/g'  -e 's/[[:blank:]]//g' > geno.temp.add
    cut  -f7- ${newname}.raw | sed -e '1d' -e 's/NA/5/g'  -e 's/2/1/g' -e 's/[[:blank:]]//g' > geno.temp.dom
    cut -f7- ${newname}.raw | sed -e '1d' -e 's/NA/5/g'  -e 's/0/1/g' -e 's/2/0/g' -e 's/[[:blank:]]//g' > geno.temp.rec
    cut -f7- ${newname}.raw | sed -e '1d' -e 's/NA/5/g'  -e 's/2/0/g' -e 's/[[:blank:]]//g' > geno.temp.dom.dev
    # count number of characters in each line. which should be your number of markers in this case 221115
    awk '{ print length }' .geno.add |uniq
    # Pasting IDs and Genotypes together
    paste Ids2.temp geno.temp.dom -d " " > ${newname}.geno.dom
    paste Ids2.temp geno.temp.rec -d " " > ${newname}.geno.rec
    paste Ids2.temp geno.temp.add -d " " > ${newname}.geno.add
    paste Ids2.temp geno.temp.dom.dev -d " " > ${newname}.geno.dom.dev
    # List of SNP name from header where the header is <SNPID>_<counted allele>. 
    head -1 ${newname}.raw | cut -d " " -f7- | sed -e 's/[[:space:]]/\n/g' > Header.map
    # Create a map file from original dataset that should be in the same
    awk '{print $2,$1,$4}' ${newname}.bim | sed -e '1i\SNPID CHR POS' > ${newname}.map
    # Craete a log file with both plink log files
    cat data.temp.log ${newname}.log > log.file
    awk -NF " " '{print $1,"0","0"}' ${newname}.geno.dom > ${newname}.PED.dom
    awk -NF " " '{print $1,"0","0"}' ${newname}.geno.dom.dev > ${newname}.PED.dom.dev
    awk -NF " " '{print $1,"0","0"}' ${newname}.geno.rec > ${newname}.PED.rec
    awk -NF " " '{print $1,"0","0"}' ${newname}.geno.add > ${newname}.PED.add
    #remove extra data 
    rm *.nosex
    rm *.log
    rm *.temp*
    rm data.temp.*
    rm Header.map
    #rm *.raw
done
