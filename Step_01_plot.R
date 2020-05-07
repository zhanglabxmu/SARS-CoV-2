remove(list=ls())
options(stringsAsFactors = F)
library(tidyverse)

nucmer<-read.table("nucmer.snps",skip=4,as.is=TRUE)
nucmer <- nucmer[,c(1,2,3,4,14)]
colnames(nucmer)<-c("rpos","rvar","qvar","qpos","ID")
dim(nucmer)

nucmer$sample <-sapply(strsplit(as.character(nucmer$ID), "[|]"), function(x) x[2])
nucmer$time <-sapply(strsplit(as.character(nucmer$ID), "[|]"), function(x) x[3])
nucmer$country <-sapply(strsplit(as.character(nucmer$ID), "[/]"), function(x) x[2])
china <-read.csv("china.txt", header = F)
china
nucmer[nucmer$country %in% china$V1, ]$country <-"China"
nucmer$M_type <-str_c(nucmer$rvar,nucmer$qvar,sep ="->")
nucmer$PM_type <-str_c(nucmer$rpos,nucmer$M_type,sep =":")

save(nucmer, file="nucmer_RMD.Rdata")

##################

remove(list=ls())
load("nucmer_RMD.Rdata")

##1 global SNP profiling
library(ggplot2)
library(ggforce)
library(ggsci)

p<-ggplot(data=nucmer,aes(x=rpos, y=sample))+
  geom_point(size=0.001, alpha=2/3,aes(color=M_type))+
  theme_bw()+
  labs(x="SARS-CoV-2 Genomic position")
ggsave(p,filename = "Global_SNP.png",width = 12, height = 8, dpi=300)
  
# +ggforce::facet_zoom(xlim=c(2000,2020))



#####################

# mutation statics for Nucleiotide

##2 Most mutated sample
MMS<-head(sort(table(nucmer$sample),dec=TRUE),n=30)
par(las=3,mar=c(15,5,5,1))
barplot(MMS,ylab="nr of mutations",main="Top30 mutated samples",col="lightblue")
ggsave(filename = "Most_Mutation_Sample.png",width = 12, height = 8, dpi=300)

##3 Averagemutation per sample
MPS<- table(table(nucmer$sample))
par(las=1,mar=c(5,5,5,2))
barplot(MPS, ylab="Case Number",xlab="NR of Mutation",main="Summarized Mutation Per Sample",col="lightblue")
ggsave(filename = "Average_Mutation.png",width = 12, height = 8, dpi=300)

# unique(nucmer$sample)
# 
# table(table(unique(nucmer$sample)))

##4 Top mutation type
##4 Top mutation type

M_type_top10 <-head(sort(table(nucmer$M_type), decreasing = T),n=10)
par(las=2,mar=c(8,5,5,2))
barplot(M_type_top10,xlab="Mutation Type",ylab="number", main="Top10 Mutation Type",col=rainbow(length(M_type_top10)))



##4.2 global SNP Profilign for Top5 major SNP

FM_nucmer<-subset(nucmer, nucmer$M_type %in% c("C->T", "A->G","G->T", "G->A", "T->C", "G->C"))
library(ggsci)
ggplot(data = FM_nucmer,aes(x=rpos, y=sample))+
  geom_point(size=0.01, alpha=2/3,aes(color=M_type))+
  theme_bw()+
  scale_color_aaas()+
  labs(x="SARS-CoV-2 Genomic position")
ggsave(filename = "Global_SNP(Major5).png",width = 12, height = 8, dpi=300)


ggplot(data=FM_nucmer,aes(x=rpos))+
  geom_point(stat="count",aes(color=M_type, alpha=12/3))+
  theme_bw()+
  scale_y_continuous(trans='log10')+
  scale_color_aaas()+
  labs(x="SARS-CoV-2 Genomic position")
ggsave(filename = "Global_SNP(Major5)#2.png",width = 12, height = 8, dpi=300)


##4.3 Most mutated postion
MMP<-head(sort(table(nucmer$rpos),decreasing=TRUE),n=30)
par(las=2,mar=c(15,5,5,1))
barplot(MMP,ylab="nr of mutations",xlab="Postion",main="Top30 mutated sites",col=rainbow(length(MMP)))

##4.4 Top30 Mutation Sites
PM_Site_top30 <-head(sort(table(nucmer$PM_type), decreasing = T),n=30)
par(las=2,mar=c(8,5,5,2))
barplot(PM_Site_top30,ylab="number", main="Top30 Mutation site",col=rainbow(length(PM_Site_top30)))


  

###### Certain sample mutation
sub_nucmer1<-subset(nucmer, nucmer$PM_type=="28883:G->C" )
sub_nucmer2<-subset(nucmer, nucmer$sample %in% sub_nucmer1$sample )
ggplot(data=sub_nucmer2,aes(x=rpos, y=sample,colour=M_type))+
  geom_point(size=1,alpha=1/3)+
  theme_bw()+
  scale_color_aaas()+
  labs(x="SARS-CoV-2 Genomic position")

barplot(sort(table(sub_nucmer2$country), decreasing = T),
        col=heat.colors(length(table(sub_nucmer2$country))))


### mutation country distribution

country_sample<-nucmer[!duplicated(nucmer$sample),c("country", "sample")]

table(country_sample$country=="China")
par(las=2,mar=c(8,5,5,2))
barplot(sort(table(country_sample$country),decreasing =T),
        col=rainbow(1:length(country_sample$sample)))


## 5 mutation density profiling by position
MbP<-as.data.frame(table(nucmer$rpos))
library(ggplot2)

ggplot(data=nucmer[],aes(x=rpos))+
  geom_density()+
  theme_bw()
# country mutation type in Top10 country
Top10_country <- c("USA", "England", "Australia","Scotland","Iceland","Wales","Netherlands", "China","Belgium","Denmark")

ggplot(data=nucmer[nucmer$country %in% Top10_country,],aes(x=rpos))+
  geom_density()+
  theme_bw()+
  facet_grid(country~ .)

ggplot(data=nucmer[nucmer$country %in% Top10_country & nucmer$rpos %in% 28831:28931,],aes(x=rpos))+
  geom_density()+
  theme_bw()+
  facet_grid(country~ .)


