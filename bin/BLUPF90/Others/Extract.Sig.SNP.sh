TRAIT=HCW
GENO=/blue/mateescu/gzayas97/Gabe_Thesis/3.BOA.GWAS/Analysis/results/SNP.BLUPF90/Genotypes/${TRAIT}/Seminole.Nutri.SNP 
BO=/blue/mateescu/gzayas97/Gabe_Thesis/3.BOA.GWAS/Analysis/results/BO.BLUPF90/Genotypes/${TRAIT}/Seminole.Nutri.BO
pheno=/blue/mateescu/gzayas97/Gabe_Thesis/3.BOA.GWAS/Analysis/results/SNP.BLUPF90/Genotypes/${TRAIT}/Cleaned.${TRAIT}.csv
GENE_PATH=/blue/mateescu/gzayas97/Gabe_Thesis/3.BOA.GWAS/Analysis/results/ALL.GWAS/HCW/Genes/BO
cd ${GENE_PATH}

cp ${pheno} .
ml plink

plink -cow -bfile ${GENO} -extract HCW.SNP.txt --make-bed --recode A -out Sig.BO.geno
plink -cow -bfile ${BO} -extract HCW.SNP.txt --make-bed --recode A  -out Sig.BO.bo


TRAIT=HCW
GENO=/blue/mateescu/gzayas97/Gabe_Thesis/3.BOA.GWAS/Analysis/results/SNP.BLUPF90/Genotypes/${TRAIT}/Seminole.Nutri.SNP 
BO=/blue/mateescu/gzayas97/Gabe_Thesis/3.BOA.GWAS/Analysis/results/BO.BLUPF90/Genotypes/${TRAIT}/Seminole.Nutri.BO
pheno=/blue/mateescu/gzayas97/Gabe_Thesis/3.BOA.GWAS/Analysis/results/SNP.BLUPF90/Genotypes/${TRAIT}/Cleaned.${TRAIT}.csv
GENE_PATH=/blue/mateescu/gzayas97/Gabe_Thesis/3.BOA.GWAS/Analysis/results/ALL.GWAS/HCW/Genes/SNP
cd ${GENE_PATH}
ml plink
cp ${pheno} .
plink -cow -bfile ${GENO} -extract HCW.SNP.txt --make-bed --recode A -out Sig.SNP.geno
plink -cow -bfile ${BO} -extract HCW.SNP.txt --make-bed --recode A  -out Sig.SNP.bo