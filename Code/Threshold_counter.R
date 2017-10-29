library(fitdistrplus)
list_file<-dir()
threshold<-c()
meant<-c()
sdt<-c()
for(i in 1:length(list_file)){
  temp<-read.table(list_file[i],sep='\t',stringsAsFactors = F,quote = "",header=T,fill = T)
  temp_contig_length<-temp$contig_length[!duplicated(temp$contig_length)]
  tempo<-fitdist(temp$contig_length,'lnorm')
  meant[i]<-tempo$estimate[1]
  sdt[i]<-tempo$estimate[2]
  thres<-temp_contig_length[order(temp_contig_length,decreasing = F)][length(temp_contig_length)*0.95]
  threshold[i]<-thres
}  
