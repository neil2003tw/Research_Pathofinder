library(ggplot2)
temp<-read.table('result_out.txt',sep='\t')
temp<-temp[!grepl('human|homo|h.sapiens',temp$V4,ignore.case = T),]
temp<-temp[!grepl('hepatitis',temp$V4,ignore.case = T),]
temp<-temp[-which(temp$V2==0),]
temp2<-data.frame(data=temp$V2,lognorm=rlnorm(758,5.5346392,0.3267225))
ggplot(temp2)+geom_density(aes(data,colour='Rawdata'))+geom_density(aes(lognorm,colour='Lognorm'))+xlab('Contig length')+scale_colour_manual(name='Distribution',values=c(Lognorm='blue',Rawdata='red'))
#ggplot(temp)+geom_line(aes(x=cov,temp$sensitivity,colour='Sensitivity'))+geom_line(aes(cov,temp$specificity,colour='Specificity'))+geom_line(aes(cov,temp$accuracy,colour='Accuracy'))+xlab('Coverage of pathogen reads')+scale_colour_manual(name=' ',values=c(Sensitivity='green',Specificity='blue',Accuracy='red'))+ylab('Value')