
remove(list=ls())
options(stringsAsFactors=F)

library(ggplot2)

load("nucmer_RMD.Rdata")

assays<-read.csv("Assays.csv")

##################
i=3
F1=assays[i,"F1"]+3  # Not aac mutation
F2=assays[i,"F2"]
R1=assays[i,"R1"]
R2=assays[i,"R2"]
P1=assays[i,"P1"]
P2=assays[i,"P2"]

sub_nucmer<-subset(nucmer, nucmer$rpos %in% c(F1:F2,P1:P2,R1:R2))
write.csv(sub_nucmer,file="ChinCDC_n-NotAAC_SNP.csv")

sub_nucmer<-subset(nucmer, nucmer$rpos==28881:28883)

table(table(sub_nucmer$sample)>=3)

write.csv(sub_nucmer,file="ChinCDC_n-NotAAC_SNP.csv")

TMN<- table(table(unique(sub_nucmer$sample)))[[1]]
Total=11100  # Total Clear GISAID fasta data
Mutation_Rate<- round(TMN/Total*100,5)

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
       title=paste0("ChinaCDC-NotAAC-Total Mutant Samples:", TMN,"/Mutation_Rate:",Mutation_Rate,"%"))

ggsave(p,filename="ChinCDC_n-NotAAC.png",width = 12, height = 8, dpi=300)

###############
