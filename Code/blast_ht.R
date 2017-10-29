library(tm)
list_file<-dir()
list_file<-list_file[grepl('_blast.out',list_file)]

r_length_mine<-function(x){
  tcratio<-x[2]/x[4]
  chosen_contigs<-x[3][tcratio==max(tcratio)]
  return(chosen_contigs)
}

duplicate_lenremover<-function(x){
  len_vec<-x[2]
  return(len_vec[1])
}

for(i in 1:length(list_file)){
  temp<-read.table(list_file[i],sep='\t',stringsAsFactors = F,quote = "",header=F,fill = T)
  temp_contig_det<-temp[grepl('human|homo',temp$V3,ignore.case = T)&!
                      grepl('virus',temp$V3,ignore.case = T)|
                      grepl('pan troglodytes|gorilla|Pongo abelii|pan paniscus|Macaca mulatta|
                            Nomascus leucogeny|Papio anubis|Chlorocebus sabaeus|Saimiri boliviensis|
                            Rhinopithecus|H.sapiens|Nomascus leucogeny|Macaca fascicularis|
                            Pan paniscus',
                            temp$V3,ignore.case = T),1]
  temp_contig<-factor(unname(sapply(temp_contig_det,function(x) strsplit(x,split = ' ')[[1]][1])))
  list_homo_contigs<-levels(temp_contig)
  list_contig_name<-unname(sapply(temp$V1,function(x) strsplit(x,split = ' ')[[1]][1]))
  temp$contig_name<-list_contig_name
  temp2<-temp[!(list_contig_name %in% list_homo_contigs),]
  #weird_id<-which(grepl('human|homo',temp$contig_simplspecies,ignore.case = T)
  #                &!grepl('virus',temp$contig_simplspecies,ignore.case = T)
  #                |grepl('pan troglodytes|gorilla|Pongo abelii|pan paniscus|Macaca mulatta|
  #                       Nomascus leucogeny|Papio anubis|Chlorocebus sabaeus|Saimiri boliviensis|
  #                       Rhinopithecus roxellana',
  #                       temp$contig_simplspecies,ignore.case = T))
  #temp2<-temp[-weird_id,]
  temp2<-temp2[order(temp2$V2,decreasing = T),]
  temp2<-temp2[temp2$V2>430,]
  text_eval<-paste0('data',i,'<-temp2')
  eval(parse(text=text_eval))
  print(i)
}

Length_table<-data.frame(contig_id=NA)
Count_table<-data.frame(contig_id=NA)
Species_table<-data.frame()

for(i in 1:length(list_file)){
  text_eval<-paste0('dat<-data',i)
  eval(parse(text=text_eval))  
  temp<-as.list(by(dat,INDICES = dat$contig_name,function(x) r_length_mine(x)))
  temp<-sapply(temp,function(x) x[1])
  temp_len<-as.list(by(dat,INDICES = dat$contig_name,function(x) duplicate_lenremover(x)))
  temp_len<-sapply(temp_len,function(x) x[1])
  temp_dat<-data.frame(contig_names=names(unlist(temp)),contig_species=unname(unlist(temp)),stringsAsFactors = F)
  contig_id<-unname(sapply(temp_dat$contig_species,function(x) strsplit(x,' ')[[1]][2]))
  contig_species<-unname(sapply(temp_dat$contig_species,
                                function(x) paste(strsplit(x,' ')[[1]][-2],collapse = ' ')))
  sample_id<-strsplit(list_file[i],split = '_')[[1]][1]
  Species_table<-rbind(Species_table,data.frame(contig_id,contig_species))
  text_eval<-paste0('current_table<-data.frame(contig_id,',sample_id,'=1)')
  eval(parse(text=text_eval))
  raw_lentable<-data.frame(contig_id,unname(unlist(temp_len)))
  temp_len2<-as.list(by(raw_lentable,INDICES = contig_id,function(x) sum(x[2])))
  text_eval<-paste0('current_lentable<-data.frame(contig_id=names(temp_len2),',sample_id,'=unname(unlist(temp_len2)))')
  eval(parse(text=text_eval))
  Count_table<-merge(Count_table,current_table, by.y = 'contig_id',all = T)
  Length_table<-merge(Length_table,current_lentable, by.y = 'contig_id',all = T)
}
Count_table<-merge(Count_table,Species_table, by.y = 'contig_id',all=T)
Length_table<-merge(Length_table,Species_table, by.y = 'contig_id',all=T)
Count_table<-Count_table[!duplicated(Count_table),]
Length_table<-Length_table[!duplicated(Length_table),]
order_list<-order(rowSums(Count_table[,2:(length(list_file)+1)],na.rm = T),decreasing = T)
Count_table<-Count_table[order_list,]
Length_table<-Length_table[ match(Count_table$contig_id,Length_table$contig_id) ,]
row.names(Count_table)<-1:dim(Count_table)[1]
row.names(Length_table)<-1:dim(Length_table)[1]

#rm(list=ls()[ls()!='Count_table'])
write.csv(Count_table,'Result_table.csv',quote=F)
write.csv(Length_table,'Length_table.csv',quote=F)


