##############
#################
## mutation in Amplicon/primer (A1: A2) : A1, A2 according to the Genomic postion of Reference SARS-CoV-2 genome

remove(list=ls())
options(stringsAsFactors=F)
library(tidyverse)
library(ggplot2)

load("nucmer_RMD.Rdata")

assays<-read.csv("Assays.csv")
mutate(assays,Mutation_Ratio=0)

for (i in 1:length(assays$Assay)){
  F1=assays[i,"F1"]
  F2=assays[i,"F2"]
  R1=assays[i,"R1"]
  R2=assays[i,"R2"]
  P1=assays[i,"P1"]
  P2=assays[i,"P2"]
  
  sub_nucmer<-subset(nucmer, nucmer$rpos %in% c(F1:F2,P1:P2,R1:R2))
  #write.csv(sub_nucmer,file=paste0(assays$Assay[i],'_SNP.csv'))
  
  TMN<- table(table(unique(sub_nucmer$sample)))[[1]]
  Total=11100  # Total Cleared GISAID fasta data
  Mutation_Ratio<- round(TMN/Total*100,5)
  assays[i,"Mutation_Ratio"]<- Mutation_Ratio

  
  p<-ggplot(data=sub_nucmer,aes(x=rpos, y=sample,color=M_type))+
    geom_point(size=2)+
    theme_bw()+
    scale_x_continuous(breaks=seq(F1,R2,2),limits =c(F1,R2))+
    geom_vline(aes(xintercept=F1))+
    geom_vline(aes(xintercept=F2))+
    geom_vline(aes(xintercept=R1))+
    geom_vline(aes(xintercept=R2))+
    geom_vline(aes(xintercept=P1))+
    geom_vline(aes(xintercept=P2))+
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+
    labs(x="SARS-CoV-2 Genomic position",
         title=paste0(assays$Assay[i],"-Total Mutant Samples:", TMN,"/Mutation_Ratio:",Mutation_Ratio,"%"))
  
   ggsave(p,filename=paste0(assays$Assay[i],'.png'),width = 12, height = 8, dpi=300)

}
######
write.csv(assays, file="Assays.csv", row.names = F)
par(las=3,mar=c(8,5,5,2))
barplot(data=assays, log2(Mutation_Ratio) ~ Assay, col=heat.colors(12))

library(ggsci)
ggplot(data=assays,aes(x=Assay,y=Mutation_Ratio,fill=Assay))+
  geom_bar(stat="identity")+
  geom_text(aes(label=Mutation_Ratio),vjust=0,angle = 90)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size=12, face="bold"),
        axis.text.y = element_text(size=12,face="bold"))+
  scale_y_continuous(trans='log2')
  #scale_y_sqrt()


#### ## 6 mutation count/type for any pcr primer

i=4
F1=assays[i,"F1"]
F2=assays[i,"F2"]
R1=assays[i,"R1"]
R2=assays[i,"R2"]
P1=assays[i,"P1"]
P2=assays[i,"P2"]


ggplot(data = nucmer[nucmer$rpos %in% c(F1:F2,P1:P2,R1:R2),],aes(x=rpos,color=M_type))+
  geom_point(stat="count")+
  theme_bw()+
  scale_y_continuous(trans='log10')+
  geom_vline(aes(xintercept=F1))+
  geom_vline(aes(xintercept=F2))+
  geom_vline(aes(xintercept=R1))+
  geom_vline(aes(xintercept=R2))+
  geom_vline(aes(xintercept=P1))+
  geom_vline(aes(xintercept=P2))+
  labs(x="SARS-CoV-2 Genomic position")
