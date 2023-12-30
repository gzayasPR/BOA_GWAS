#!/bin/bash
#SBATCH --account=mateescu
#SBATCH --job-name=Data
#SBATCH --output=Data_%j.out
#SBATCH --error=errorfile_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=7000mb
data=/blue/mateescu/gzayas97/Gabe_Thesis/4.BO.GWAS.Heterosis/Analysis/data
geno=/blue/mateescu/gzayas97/UF_ILUMINA.ID_OCT23
geno_name=UFID_UF_ILUMINA.ID_OCT23_Illumina.ID_ARS
BO=/blue/mateescu/gzayas97/Gabe_Thesis/2.Pop.Structure.GBC/Analysis/results/Breed.of.Origin/Breed_of_Origin.files/
BO_name=UF.Geno_Jan23.BO
cd $data
mkdir IDs
cd $data
mkdir Genotypes
cd ./Genotypes
mkdir SNP
mkdir BOA
cd $data/Genotypes/SNP
cp $geno/${geno_name}.* .
ml plink
plink -cow -bfile ${geno_name} -chr 1-29  -keep $data/IDs/Seminole.Steers.ID -make-bed --keep-allele-order --out Seminole.Nutritional
cd $data/Genotypes/BOA
cp $BO/${BO_name}.* .
ml plink
plink -cow -bfile ${BO_name} -chr 1-29  -keep $data/IDs/Seminole.Steers.ID -make-bed --keep-allele-order --out Seminole.Nutritional.BO
plink -cow -bfile ${BO_name}  -chr 1-29  -keep $data/IDs/Seminole.Steers.ID --geno-counts -out SemiNutri.BO
#Setting up phenotype files
cd $data/
mkdir Phenotypes
cd $data/Phenotypes
cp $data/IDs/SemiNutri.Phentoypes.csv .
geno_name=Seminole.Nutritional

awk -F "," '{print $2}' SemiNutri.Phentoypes.csv   > id.temp
cp  $data/Genotypes/SNP/${geno_name}.fam FAM.file
grep -wrf id.temp FAM.file | awk '{print $1,$2}' > id2.temp
ml plink
#Quality Controls
plink -cow -bfile $data/Genotypes/SNP/${geno_name}  -chr 1-29  -keep id2.temp -maf 0.01 -geno 0.1 -mind 0.4  -pca 2
plink -cow -bfile $data/Genotypes/SNP/${geno_name}  -chr 1-29  -keep id2.temp --geno-counts -out SemiNutri
rm FAM.file
rm *.log 
rm *.temp
rm *.irem
rm *.nosex
cp SemiNutri.Phentoypes.csv Cleaned.SemiNutri.Phentoypes.csv
# Manually add PCA components and clean data empty spaces,0,#N/A --> -999 (USE VLOOKUP AND REPLACE FUNCTIONS IN EXCEL)
cp Cleaned.SemiNutri.Phentoypes.csv Cleaned.HCW.csv
cp Cleaned.SemiNutri.Phentoypes.csv Cleaned.MARB.csv
#Delete rows with missing values
