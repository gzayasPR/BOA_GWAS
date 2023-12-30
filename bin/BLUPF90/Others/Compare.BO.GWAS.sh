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
mkdir ./Compare.BO.GWAS
cd ./Compare.BO.GWAS

cp ${results}/BO.Dom/chrsnp_pval chrsnp_pval.Dom
cp ${results}/BO.Rec/chrsnp_pval chrsnp_pval.Rec
cp ${results}/BO.Dom.Dev/chrsnp_pval chrsnp_pval.Dom.Dev
cp ${results}/Genotypes/BO/Seminole.Nutri.BO.map .

cp ${results}/Geno.Add/chrsnp_pval chrsnp_pval.Geno.Add
cp ${results}/Geno.Add/Seminole.Nutri.map .

ml plink/1.90b3
plink -cow -bfile ${results}/Genotypes/BO/Seminole.Nutri.BO -extract SNP.Sig.txt --recodeA  -out Sig.SNPs
cp ${results}/Genotypes/BO/HCW.csv .




Dom.p <- read.table("chrsnp_pval.Dom")[,c(3,5,6)]
names(Dom.p) <- c("Dominance Brahman","CHR","POS")
Dom.p <- Dom.p[!duplicated(Dom.p[,c('CHR','Dominance Brahman')]),]

Rec.p <- read.table("chrsnp_pval.Rec")[,c(3,5,6)]
names(Rec.p) <- c("Dominance Angus","CHR","POS")
Rec.p <- Rec.p[!duplicated(Rec.p[,c('CHR','Dominance Angus')]),]

Dom.Dev.p <- read.table("chrsnp_pval.Dom.Dev")[,c(3,5,6)]
names(Dom.Dev.p) <- c("Dominance Deviation","CHR","POS")

Dom.Dev.p <- Dom.Dev.p[!duplicated(Dom.Dev.p[,c('CHR','Dominance Deviation')]),]
MAP <-  read.table("Seminole.Nutri.BO.map",header=T)
names(MAP)[1] <- "SNP"
#put all data frames into list
df_list <- list( Rec.p,Dom.p, Dom.Dev.p,MAP)

#merge all data frames in list
all <- Reduce(function(x, y) merge(x, y,by = c("CHR","POS"),all.x=T), df_list)

all.df <- all[,c(6,1:5)]


write.csv(all.df,"Dominance.pvalues.csv",quote=F,row.names=F)

library(CMplot)

#CMplot(all.df,col=c("orange","blue","green"),multracks=T,LOG10=F)
CMplot(all.df,plot.type = "m",col=c("orange","blue"),multracks=T,LOG10=F,height=5,width=22)

Geno.Add <- read.table("chrsnp_pval.Geno.Add")[,c(3,5,6)]
MAP.geno <-  read.table("Seminole.Nutri.map",header=T)
names(MAP.geno)[1] <- "SNP"
names(Geno.Add ) <- c("Additive","CHR","POS")
add.df <- merge(Geno.Add,MAP.geno,by = c("CHR","POS"))
library(CMplot)
CMplot(add.df ,plot.type = "m",col=c("orange","blue"),LOG10=F)

Sig.SNPs <- read.table("Sig.SNPs.raw",header=T)
HCW <- read.csv("HCW.csv",header=T)
names(HCW)[1] <- "IID" 




df <- merge(Sig.SNPs,HCW,by="IID")
library(emmeans)
lm.hcw <- lm(hcw ~ as.factor(X2.95687817.C.T.rs210039354_B) + as.factor(CG) + as.factor(RanchLOT),data = na.omit(df))
lsmeans(lm.hcw , pairwise ~as.factor(X2.95687817.C.T.rs210039354_B))
table(df$X2.95687817.C.T.rs210039354_B)

lm.hcw <- lm(hcw ~ as.factor(X12.12480612.C.T.rs110854900_B) + as.factor(CG) + as.factor(RanchLOT),data = na.omit(df))
lsmeans(lm.hcw , pairwise ~as.factor(X12.12480612.C.T.rs110854900_B))
table(df$X12.12480612.C.T.rs110854900_B)


lm.hcw <- lm(hcw ~ as.factor(X23.43024105.C.T.rs110862272_B) + as.factor(CG) + as.factor(RanchLOT),data = na.omit(df))
lsmeans(lm.hcw , pairwise ~as.factor(X23.43024105.C.T.rs110862272_B))
table(df$X23.43024105.C.T.rs110862272_B)



lm.hcw <- lm(hcw ~ as.factor(X23.44015434.G.T.rs211088453_B) + as.factor(CG) + as.factor(RanchLOT),data = na.omit(df))
lsmeans(lm.hcw , pairwise ~as.factor(X23.44015434.G.T.rs211088453_B))
table(df$X23.44015434.G.T.rs211088453_B)

lm.hcw <- lm(hcw ~ as.factor(X19.39013681.T.C.rs41920953_B) + as.factor(CG) + as.factor(RanchLOT),data = na.omit(df))
lsmeans(lm.hcw , pairwise ~as.factor(X19.39013681.T.C.rs41920953_B))
table(df$X19.39013681.T.C.rs41920953_B)

lm.hcw <- lm(hcw ~ as.factor(BovineHD0100025146_B) + as.factor(CG) + as.factor(RanchLOT),data = na.omit(df))
lsmeans(lm.hcw , pairwise ~as.factor(BovineHD0100025146_B))
table(df$BovineHD0100025146_B)

lm.hcw <- lm(hcw ~ as.factor(X29.28578769.C.T.rs42179569_B) + as.factor(CG) + as.factor(RanchLOT),data = na.omit(df))
lsmeans(lm.hcw , pairwise ~as.factor(X29.28578769.C.T.rs42179569_B))
table(df$X29.28578769.C.T.rs42179569_B)