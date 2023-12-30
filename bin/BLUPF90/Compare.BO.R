setwd("/blue/mateescu/gzayas97/Gabe_Thesis/3.BOA.GWAS/Analysis/results/ALL.GWAS/HCW/SNP/")
####SNP
SNP.Dom.p <- read.table("SNP.p_val.HCW.dom.txt")[,c(3,5,6)]
names(SNP.Dom.p) <- c("Dominance Major Allele","CHR","POS")
SNP.Dom.p <- SNP.Dom.p[!duplicated(SNP.Dom.p[,c('CHR','Dominance Major Allele')]),]
SNP.Rec.p <- read.table("SNP.p_val.HCW.rec.txt")[,c(3,5,6)]
names(SNP.Rec.p) <- c("Dominance Minor Allele","CHR","POS")
SNP.Rec.p <- SNP.Rec.p[!duplicated(SNP.Rec.p[,c('CHR','Dominance Minor Allele')]),]
MAP <-  read.table("Seminole.Nutri.SNP.map",header=T)
SNP.add.p <- read.table("SNP.p_val.HCW.add.txt")[,c(3,5,6)]
names(SNP.add.p ) <- c("Additive","CHR","POS")
SNP.add.p <- SNP.add.p[!duplicated(SNP.add.p[,c('CHR','Additive')]),]
names(MAP)[1] <- "SNP"
#put all data frames into list
SNP.df_list <- list( SNP.add.p ,SNP.Rec.p,SNP.Dom.p,MAP)
#merge all data frames in list
SNP.all <- Reduce(function(x, y) merge(x, y,by = c("CHR","POS"),all.x=T), SNP.df_list)
SNP.all.df <- SNP.all[,c(6,1:5)]
write.csv(SNP.all.df ,"SNP.HCW.pvalues.csv",quote=F,row.names=F)



SNP.Dom.var <- read.table("SNP.var.10.HCW.dom.txt")[,c(3,5,6)]
names(SNP.Dom.var) <- c("Dominance Major Allele","CHR","POS")
SNP.Dom.var <- SNP.Dom.var[!duplicated(SNP.Dom.var[,c('CHR','Dominance Major Allele')]),]
SNP.rec.var <- read.table("SNP.var.10.HCW.rec.txt")[,c(3,5,6)]
names(SNP.rec.var) <- c("Dominance Minor Allele","CHR","POS")
SNP.rec.var <- SNP.rec.var[!duplicated(SNP.rec.var[,c('CHR','Dominance Minor Allele')]),]
MAP <-  read.table("Seminole.Nutri.SNP.map",header=T)
SNP.add.var <- read.table("SNP.var.10.HCW.add.txt")[,c(3,5,6)]
names(SNP.add.var ) <- c("Additive","CHR","POS")
SNP.add.var <- SNP.add.var[!duplicated(SNP.add.var[,c('CHR','Additive')]),]
names(MAP)[1] <- "SNP"
#put all data frames into list
SNP.df_list <- list( SNP.add.var ,SNP.rec.var,SNP.Dom.var,MAP)
#merge all data frames in list
SNP.all <- Reduce(function(x, y) merge(x, y,by = c("CHR","POS"),all.x=T), SNP.df_list)
SNP.all.df <- SNP.all[,c(6,1:5)]
write.csv(SNP.all.df ,"SNP.HCW.var10.csv",quote=F,row.names=F)



setwd("/blue/mateescu/gzayas97/Gabe_Thesis/3.BOA.GWAS/Analysis/results/ALL.GWAS/HCW/BO/")
####BO
BO.Dom.p <- read.table("BO.p_val.HCW.dom.txt")[,c(3,5,6)]
names(BO.Dom.p) <- c("Dominance Brahman","CHR","POS")
BO.Dom.p <- BO.Dom.p[!duplicated(BO.Dom.p[,c('CHR','Dominance Brahman')]),]
BO.Rec.p <- read.table("BO.p_val.HCW.rec.txt")[,c(3,5,6)]
names(BO.Rec.p) <- c("Dominance Angus","CHR","POS")
BO.Rec.p <- BO.Rec.p[!duplicated(BO.Rec.p[,c('CHR','Dominance Angus')]),]
MAP <-  read.table("Seminole.Nutri.BO.map",header=T)
BO.add.p <- read.table("BO.p_val.HCW.add.txt")[,c(3,5,6)]
names(BO.add.p ) <- c("Additive","CHR","POS")
BO.add.p <- BO.add.p[!duplicated(BO.add.p[,c('CHR','Additive')]),]
names(MAP)[1] <- "BO"
#put all data frames into list
BO.df_list <- list( BO.add.p ,BO.Rec.p,BO.Dom.p,MAP)
#merge all data frames in list
BO.all <- Reduce(function(x, y) merge(x, y,by = c("CHR","POS"),all.x=T), BO.df_list)
BO.all.df <- BO.all[,c(6,1:5)]
write.csv(BO.all.df ,"BO.HCW.pvalues.csv",quote=F,row.names=F)

BO.Dom.var <- read.table("BO.var.10.HCW.dom.txt")[,c(3,5,6)]
names(BO.Dom.var) <- c("Dominance Brahman","CHR","POS")
BO.Dom.var <- BO.Dom.var[!duplicated(BO.Dom.var[,c('CHR','Dominance Brahman')]),]
BO.rec.var <- read.table("BO.var.10.HCW.rec.txt")[,c(3,5,6)]
names(BO.rec.var) <- c("Dominance Angus","CHR","POS")
BO.rec.var <- BO.rec.var[!duplicated(BO.rec.var[,c('CHR','Dominance Angus')]),]
MAP <-  read.table("Seminole.Nutri.BO.map",header=T)
BO.add.var <- read.table("BO.var.10.HCW.add.txt")[,c(3,5,6)]
names(BO.add.var ) <- c("Additive","CHR","POS")
BO.add.var <- BO.add.var[!duplicated(BO.add.var[,c('CHR','Additive')]),]
names(MAP)[1] <- "SNP"
#put all data frames into list
BO.df_list <- list( BO.add.var ,BO.rec.var,BO.Dom.var,MAP)
#merge all data frames in list
BO.all <- Reduce(function(x, y) merge(x, y,by = c("CHR","POS"),all.x=T), BO.df_list)
BO.all.df <- BO.all[,c(6,1:5)]
write.csv(BO.all.df ,"BO.HCW.var10.csv",quote=F,row.names=F)
