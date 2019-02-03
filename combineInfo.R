#!/usr/bin/R

setwd("~/Desktop/")

### 读入annovar注释文件
tcga<-read.table("./tcga.crc.exonic_variant_function",sep="\t")
icgc<-read.table("./icgc.crc.exonic_variant_function",sep="\t")
crc111<-read.table("./111crc.crc.exonic_variant_function",sep="\t")
tcganew<-read.table("./tcgaNew.exonic_variant_function",sep="\t")


###读入已有的信息文件
tcga.info<-read.table("./tcga.info",header=T,row.names=1)
icgc.info<-read.table("./icgc.info",header=T,row.names=1)
crc111.info<-read.table("./111crc.info",header=T,row.names=1)
tcganew.info<-read.table("./tcganew.info",header=T,row.names=1)



colnames(tcga)<-c("Mut.Type","INFO","Chr","Start","End","Ref","Alt","Sample")
colnames(icgc)<-colnames(tcga)
colnames(crc111)<-colnames(tcga)
colnames(tcganew)<-colnames(tcga)

###新增一列值Eff，分为Truncated,Missense和synonymous三种
tcga$Eff[tcga$Mut.Type%in%c("frameshift deletion","stopgain","frameshift insertion","frameshift substitution","splicing")]="Truncated"
tcga$Eff[tcga$Mut.Type%in%c("nonsynonymous SNV","nonframeshift deletion","nonframeshift substitution")]="Missense"
tcga$Eff[tcga$Mut.Type%in%c("synonymous SNV")]="synonymous"

icgc$Eff[icgc$Mut.Type%in%c("frameshift deletion","stopgain","frameshift insertion","frameshift substitution","splicing")]="Truncated"
icgc$Eff[icgc$Mut.Type%in%c("nonsynonymous SNV","nonframeshift deletion","nonframeshift substitution")]="Missense"
icgc$Eff[icgc$Mut.Type%in%c("synonymous SNV")]="synonymous"

crc111$Eff[crc111$Mut.Type%in%c("frameshift deletion","stopgain","frameshift insertion","frameshift substitution","splicing")]="Truncated"
crc111$Eff[crc111$Mut.Type%in%c("nonsynonymous SNV","nonframeshift deletion","nonframeshift substitution")]="Missense"
crc111$Eff[crc111$Mut.Type%in%c("synonymous SNV")]="synonymous"


tcganew$Eff[tcganew$Mut.Type%in%c("frameshift deletion","stopgain","frameshift insertion","frameshift substitution","splicing")]="Truncated"
tcganew$Eff[tcganew$Mut.Type%in%c("nonsynonymous SNV","nonframeshift deletion","nonframeshift substitution")]="Missense"
tcganew$Eff[tcganew$Mut.Type%in%c("synonymous SNV")]="synonymous"

### 将INFO分成多列("Gene","Transcript","Exon","ntChange","aaChange")
newd.tcga<-matrix(unlist(lapply(strsplit(as.character(tcga$INFO),"[:,]"),function(x)x[1:5])),nrow=nrow(tcga),ncol=5,byrow=T)
colnames(newd.tcga)<-c("Gene","Transcript","Exon","ntChange","aaChange")
tcga<-cbind(tcga,newd.tcga)

newd.icgc<-matrix(unlist(lapply(strsplit(as.character(icgc$INFO),"[:,]"),function(x)x[1:5])),nrow=nrow(icgc),ncol=5,byrow=T)
colnames(newd.icgc)<-colnames(newd.tcga)
icgc<-cbind(icgc,newd.icgc)

newd.crc111<-matrix(unlist(lapply(strsplit(as.character(crc111$INFO),"[:,]"),function(x)x[1:5])),nrow=nrow(crc111),ncol=5,byrow=T)
colnames(newd.crc111)<-colnames(newd.tcga)
crc111<-cbind(crc111,newd.crc111)

newd.tcganew<-matrix(unlist(lapply(strsplit(as.character(tcganew$INFO),"[:,]"),function(x)x[1:5])),nrow=nrow(tcganew),ncol=5,byrow=T)
colnames(newd.tcganew)<-c("Gene","Transcript","Exon","ntChange","aaChange")
tcganew<-cbind(tcganew,newd.tcganew)

genes<-c("PMS2","MSH2","MLH1","MSH6","TP53","APC","SMAD4","KRAS","PIK3CA","RNF43","BRAF","TCF7L2","PTEN","FAM189A1","TFAP2D","NGF","ARHGAP15","CHST1","BAX")

genes<-c("IRF2BPL","TMEM67","CEACAM6","ZFP36L2","ALPI","FRG2C","CSMD1","MAG")

### smg mutations (whether mutated, aaChange, Effect)
#SMGs中某个基因突变所在行
smg.row.tcga<-sapply(genes,function(x){which(tcga$Gene==x)})
#SMGs突变表格(1:突变,0:野生型)
smg.mut.tcga<-data.frame(lapply(smg.row.tcga,function(x){ifelse(rownames(tcga.info)%in%tcga$Sample[x],1,0)}))
#SMGs中某个基因的突变表格(突变样本；氨基酸改变；突变影响)
smg.aaChange.tcga<-lapply(smg.row.tcga,function(x){data.frame(tcga$Sample[x],tcga$aaChange[x],tcga$Eff[x])})
#SMGs中某个基因在一个人里面发生多种突变，合并这些突变成一行
smg.aaChangeMod.tcga<-lapply(smg.aaChange.tcga,function(x){
  t( sapply(unique(x$tcga.Sample.x.),function(x1){row<-which(x$tcga.Sample.x.%in%x1);c(as.character(x1),
  paste(as.character(x$tcga.aaChange.x.[row]),collapse =";"),paste(as.character(x$tcga.Eff.x.[row]),collapse=";"))}))
   })
smg.tcga<- data.frame(lapply(smg.aaChangeMod.tcga,function(x){ t(sapply(rownames(tcga.info),function(x1){row<-which(x[,1]==x1); if(length(row)>0){x[row,2:3] }else{rep(NA,2)} }) )}))

colnames(smg.tcga)<-paste(rep(genes,each=2),c("aaChange","effect"),sep = "-")

tcga.info.smg<-cbind(tcga.info,smg.mut.tcga,smg.tcga)



smg.row.icgc<-sapply(genes,function(x){which(icgc$Gene==x)})
smg.mut.icgc<-data.frame(lapply(smg.row.icgc,function(x){ifelse(rownames(icgc.info)%in%icgc$Sample[x],1,0)}))
smg.aaChange.icgc<-lapply(smg.row.icgc,function(x){data.frame(icgc$Sample[x],icgc$aaChange[x],icgc$Eff[x])})
smg.aaChangeMod.icgc<-lapply(smg.aaChange.icgc,function(x){
  t( sapply(unique(x$icgc.Sample.x.),function(x1){row<-which(x$icgc.Sample.x.%in%x1);c(as.character(x1),paste(as.character(x$icgc.aaChange.x.[row]),collapse =";"),paste(as.character(x$icgc.Eff.x.[row]),collapse=";"))}))
})
smg.icgc<-data.frame( lapply(smg.aaChangeMod.icgc,function(x){ t(sapply(rownames(icgc.info),function(x1){row<-which(x[,1]==x1); if(length(row)>0){x[row,] }else{rep(NA,3)} }) )}))

colnames(smg.icgc)<-colnames(smg.tcga)

icgc.info.smg<-cbind(icgc.info,smg.mut.icgc,smg.icgc)



smg.row.crc111<-sapply(genes,function(x){which(crc111$Gene==x)})
smg.mut.crc111<-data.frame(lapply(smg.row.crc111,function(x){ifelse(rownames(crc111.info)%in%crc111$Sample[x],1,0)}))
smg.aaChange.crc111<-lapply(smg.row.crc111,function(x){data.frame(crc111$Sample[x],crc111$aaChange[x],crc111$Eff[x])})
smg.aaChangeMod.crc111<-lapply(smg.aaChange.crc111,function(x){
  t( sapply(unique(x$crc111.Sample.x.),function(x1){row<-which(x$crc111.Sample.x.%in%x1);c(as.character(x1),paste(as.character(x$crc111.aaChange.x.[row]),collapse =";"),paste(as.character(x$crc111.Eff.x.[row]),collapse=";"))}))
})
smg.crc111<-data.frame( lapply(smg.aaChangeMod.crc111,function(x){ t(sapply(rownames(crc111.info),function(x1){row<-which(x[,1]==x1); if(length(row)>0){x[row,] }else{rep(NA,3)} }) )}))

smg.crc111<- data.frame(lapply(smg.aaChangeMod.crc111,function(x){ t(sapply(rownames(crc111.info),function(x1){row<-which(x[,1]==x1); if(length(row)>0){x[row,2:3] }else{rep(NA,2)} }) )}))

colnames(smg.crc111)<-colnames(smg.tcga)

crc111.info.smg<-cbind(crc111.info,smg.mut.crc111,smg.crc111)


smg.row.tcganew<-sapply(genes,function(x){which(tcganew$Gene==x)})
#SMGs突变表格(1:突变,0:野生型)
smg.mut.tcganew<-data.frame(lapply(smg.row.tcganew,function(x){ifelse(rownames(tcganew.info)%in%tcganew$Sample[x],1,0)}))
#SMGs中某个基因的突变表格(突变样本；氨基酸改变；突变影响)
smg.aaChange.tcganew<-lapply(smg.row.tcganew,function(x){data.frame(tcganew$Sample[x],tcganew$aaChange[x],tcganew$Eff[x])})
#SMGs中某个基因在一个人里面发生多种突变，合并这些突变成一行
smg.aaChangeMod.tcganew<-lapply(smg.aaChange.tcganew,function(x){
  t( sapply(unique(x$tcganew.Sample.x.),function(x1){row<-which(x$tcganew.Sample.x.%in%x1);c(as.character(x1),
                                                                                       paste(as.character(x$tcganew.aaChange.x.[row]),collapse =";"),paste(as.character(x$tcganew.Eff.x.[row]),collapse=";"))}))
})
smg.tcganew<- data.frame(lapply(smg.aaChangeMod.tcganew,function(x){ t(sapply(rownames(tcganew.info),function(x1){row<-which(x[,1]==x1); if(length(row)>0){x[row,2:3] }else{rep(NA,2)} }) )}))

colnames(smg.tcganew)<-paste(rep(genes,each=2),c("aaChange","effect"),sep = "-")

tcganew.info.smg<-cbind(tcganew.info,smg.mut.tcganew,smg.tcganew)





write.xlsx2(crc111.info.smg,"crc.combined.xlsx",showNA=T,sheetName = "111crc")
write.xlsx2(tcga.info.smg,"crc.combined.xlsx",showNA=T,sheetName = "tcga",append = T)
write.xlsx2(icgc.info.smg,"crc.combined.xlsx",showNA=T,sheetName = "icgc",append = T)
write.xlsx2(tcganew.info.smg,"crc.combined.xlsx",showNA=T,sheetName = "tcga_new",append = T)


write.xlsx2(crc111.info.smg,"crc.combined.picked.xlsx",showNA=T,sheetName = "111crc")
write.xlsx2(tcga.info.smg,"crc.combined.picked.xlsx",showNA=T,sheetName = "tcga",append = T)
write.xlsx2(icgc.info.smg,"crc.combined.picked.xlsx",showNA=T,sheetName = "icgc",append = T)
