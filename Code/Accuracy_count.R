list_file<-dir()
accuracy<-c()
specificity<-c()
sensitivity<-c()
for(i in 1:length(list_file)){
  temp<-read.table(list_file[i],sep='\t',stringsAsFactors = F,quote = "",header=T,fill = T)
  temp<-temp[!(grepl('human|homo',temp$contig_simplspecies,ignore.case = T)&!grepl('virus',temp$contig_simplspecies,ignore.case = T)
             |grepl('pan troglodytes|gorilla|Pongo abelii|pan paniscus',temp$contig_simplspecies,ignore.case = T)),]
  temp_Positive<-temp[temp$contig_length>=430,]
  Positive<-dim(temp_Positive)[1]
  TP<-dim(temp_Positive[grepl('Pseudo',temp_Positive$contig_simplspecies,ignore.case = T),])[1]
  temp_Negative<-temp[temp$contig_length<430,]
  Negative<-dim(temp_Negative)[1]
  TN<-dim(temp_Negative[!grepl('Pseudo',temp_Negative$contig_simplspecies,ignore.case = T),])[1]
  accuracy[i]<-(TP+TN)/(Positive+Negative)
  sensitivity[i]<-TP/(TP+Negative-TN)
  specificity[i]<-TN/(TN+Positive-TP)
  names(accuracy)[i]<-list_file[i]
}  
