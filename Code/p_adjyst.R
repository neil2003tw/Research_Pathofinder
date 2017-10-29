lookuptable<-data.frame(row=NA,col=NA,p=NA)

pb <- txtProgressBar(min = 0, max = 98 , style = 3)
for(i in 2:9){
  setTxtProgressBar(pb, i)
  for(j in 1:286){
    cl<-as.numeric(Length_table[j,i])
    if(cl!=0){
      p<-pweibull(cl,2.633996,301.2119)
      newline<-data.frame(row=j,col=i,p=(1-p))
      lookuptable<-rbind(lookuptable,newline)
    }
  }
}
lookuptable$p<-p.adjust(p = lookuptable$p,'BH')
lookuptable<-lookuptable[-1,]
P_table<-Length_table

for(i in 2:9){
  P_table[P_table[,i]==0,i]<-1
}

for(i in 1:nrow(lookuptable)){
  P_table[lookuptable$row[i],lookuptable$col[i]]<-lookuptable$p[i]
}

write.table(P_table,'P_table.txt',sep='\t')

for(i in 1:sum(lookuptable$p>0.05)){
  cl<-which(lookuptable$p>0.05)[i]
  print(Length_table[ lookuptable$row[cl],lookuptable$col[cl] ])
}

