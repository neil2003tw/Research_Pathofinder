library(fitdistrplus)
library(ggplot2)
temp<-read.table('result_out.txt',sep='\t')
temp<-temp[!grepl('human|homo|h.sapiens',temp$V4,ignore.case = T),]
temp<-temp[!grepl('hepatitis',temp$V4,ignore.case = T),]
temp<-temp[-which(temp$V2==0),]
contig_length<-temp$V2
contig_length<-c(contig_length,seq(200,400,10),seq(300,400,5))
contig_length<-contig_length[order(contig_length,decreasing = T)]
pvalue<-(1:length(contig_length))/length(contig_length)
lookuptable<-data.frame(contig_length,pvalue)

fitlnorm<-fitdist(contig_length,'lnorm')
fitweibull<-fitdist(contig_length,'weibull')
fitgamma<-fitdist(contig_length,'gamma')
fitcauchy<-fitdist(contig_length,'cauchy')
fitnbinom<-fitdist(contig_length,'nbinom')
fitlogis<-fitdist(contig_length,'logis')
#fitpois<-fitdist(contig_length,'pois')

re_table<-data.frame()

for(i in 1:700){
  relnorm<-(1-plnorm(lookuptable[i,1],meanlog = fitlnorm$estimate[1],sdlog = fitlnorm$estimate[2])-lookuptable[i,2])/lookuptable[i,2]
  reweibull<-(1-pweibull(lookuptable[i,1],shape = fitweibull$estimate[1],scale = fitweibull$estimate[2])-lookuptable[i,2])/lookuptable[i,2]
  regamma<-(1-pgamma(lookuptable[i,1],shape = fitgamma$estimate[1],rate = fitgamma$estimate[2])-lookuptable[i,2])/lookuptable[i,2]
  recauchy<-(1-pcauchy(lookuptable[i,1],fitcauchy$estimate[1],fitcauchy$estimate[2])-lookuptable[i,2])/lookuptable[i,2]
  renbinom<-(1-pnbinom(lookuptable[i,1],fitnbinom$estimate[1],mu = fitnbinom$estimate[2])-lookuptable[i,2])/lookuptable[i,2]
  relogis<-(1-plogis(lookuptable[i,1],fitlogis$estimate[1],fitlogis$estimate[2])-lookuptable[i,2])/lookuptable[i,2]
  temp_table<-data.frame(reweibull=abs(reweibull)*100,relnorm=abs(relnorm)*100,regamma=abs(regamma)*100,
                         recauchy=abs(recauchy)*100,renbinom=abs(renbinom)*100,
                         relogis=abs(relogis)*100,pval=lookuptable[i,2])
  
  re_table<-rbind(re_table,temp_table)
}

ggplot(re_table,aes(x=pval))+geom_line(aes(y=relnorm,colour='lnorm'))+geom_line(aes(y=regamma,colour='gamma'))+
  geom_line(aes(y=reweibull,colour='weibull'))+geom_line(aes(y=renbinom,colour='nbinom'))+
  geom_line(aes(y=recauchy,colour='cauchy'))+geom_line(aes(y=relogis,colour='logistic'))+
  scale_colour_manual(name='Distribution',values=c(lnorm='blue',gamma='red',weibull='green',nbinom='purple',cauchy='yellow',logistic='orange'))+
  xlim(0.06,0.04)+ylim(0,100)+geom_vline(x=0.05)+ylab('Relative Error (%)')+xlab('p-value')


