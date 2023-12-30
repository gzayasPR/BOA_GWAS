setwd("/blue/mateescu/gzayas97/Gabe_Thesis/3.BOA.GWAS/Analysis/results/ALL.GWAS/MARB/SNP")
MAP <-  read.table("Seminole.Nutri.SNP.map",header=T)
names(MAP)[1] <- "SNP"
SNP.add.var <- read.table("SNP.var.200.MARB.add.txt")[,c(3,5,6)]
names(SNP.add.var ) <- c("Additive","CHR","POS")
SNP.add.var <- SNP.add.var[!duplicated(SNP.add.var[,c('CHR','Additive')]),]
SNP.Dom.var <- read.table("SNP.var.200.MARB.dom.txt")[,c(3,5,6)]
names(SNP.Dom.var) <- c("Dominance Major Allele","CHR","POS")
SNP.Dom.var <- SNP.Dom.var[!duplicated(SNP.Dom.var[,c('CHR','Dominance Major Allele')]),]
SNP.rec.var <- read.table("SNP.var.200.MARB.rec.txt")[,c(3,5,6)]
names(SNP.rec.var) <- c("Dominance Minor Allele","CHR","POS")
SNP.rec.var <- SNP.rec.var[!duplicated(SNP.rec.var[,c('CHR','Dominance Minor Allele')]),]
#put all data frames into list
SNP.df_list <- list( SNP.add.var ,SNP.Dom.var,SNP.rec.var,MAP)
#merge all data frames in list
SNP.all <- Reduce(function(x, y) merge(x, y,by = c("CHR","POS"),all.x=T), SNP.df_list)
SNP.all.df <- SNP.all[,c(6,1:5)]
write.csv(SNP.all.df ,"SNP.MARB.var200.csv",quote=F,row.names=F)

library(CMplot)

CMplot(SNP.all.df,col=c("orange","blue","green"),multracks=T,LOG10=F)

CMplot(SNP.all.df,plot.type = "m",col=c("orange","blue"),multracks=T,LOG10=F,height=5,width=22)
