library(tidyverse)
options(stringsAsFactors = F)

china <-read.csv("china.txt", header = F)
china

AAC_country<- read.delim("AAC-country-infor.txt",header=F)
AAC_country[AAC_country$V1 %in% china$V1, "V1"]<-"China"
AAC_country<- as.data.frame(sort(table(AAC_country$V1),decreasing=T))
names(AAC_country)<-c("Country","AAC_number")
write.csv(AAC_country,file="AAC_country.csv")
barplot(data=AAC_country, AAC_number ~ Country,
        col=rainbow(length(AAC_country$Country)),
        xlab="")

all_country<-read.delim("country-infor.txt",header=F)
all_country[all_country$V1 %in% china$V1, "V1"]<-"China"
all_country<-as.data.frame(sort(table(all_country$V1),decreasing=T))
names(all_country)<-c("Country","Number")
write.csv(all_country,file="All_country.csv")
barplot(data=all_country, Number ~ Country, 
        col=heat.colors(length(all_country$Country)),
        xlab="")

AAC<-left_join(AAC_country, all_country)
AAC<- AAC[1:10,]
AAC$AAC_Ratio=AAC$AAC_number/AAC$Number
AAC<-AAC[order(AAC$AAC_Ratio,decreasing = T),]
par(las=2,mar=c(8,5,5,2))
barplot(data=AAC, AAC_Ratio ~ Country,
        col=rainbow(length(AAC$Country)),
        xlab="")
