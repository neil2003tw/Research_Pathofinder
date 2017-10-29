library(tm)
source('/Users/NeilWu/Github/Research_Pathogen_war/Code/tool.R')
list_file<-dir()
list_file<-list_file[grepl('_out',list_file)]


for(i in 1:(length(list_file))){
  temp<-read.table(list_file[i],sep='\t',stringsAsFactors = F,quote = "",header=T,fill = T)
  temp_contig<-temp[(grepl('human|homo|H.sapiens',temp$contig_simplspecies,ignore.case = T)&!
                    grepl('virus',temp$contig_simplspecies,ignore.case = T))|
                    grepl('pan troglodytes|gorilla|Pongo abelii|pan paniscus|Macaca mulatta|
                           Nomascus leucogeny|Papio anubis|Chlorocebus sabaeus|Saimiri boliviensis|
                           Rhinopithecus|Nomascus leucogeny|Macaca fascicularis|
                           ',
                           temp$contig_simplspecies,ignore.case = T),1]
  temp_contig<-factor(temp_contig)
  list_homo_contigs<-levels(temp_contig)
  temp2<-temp[!(temp$contig_name %in% list_homo_contigs),]
  #weird_id<-which(grepl('human|homo',temp$contig_simplspecies,ignore.case = T)
  #                &!grepl('virus',temp$contig_simplspecies,ignore.case = T)
  #                |grepl('pan troglodytes|gorilla|Pongo abelii|pan paniscus|Macaca mulatta|
  #                       Nomascus leucogeny|Papio anubis|Chlorocebus sabaeus|Saimiri boliviensis|
  #                       Rhinopithecus roxellana',
  #                       temp$contig_simplspecies,ignore.case = T))
  #temp2<-temp[-weird_id,]
  temp2<-temp2[order(temp2$contig_length,decreasing = T),]
  temp2<-temp2[temp2$contig_length>430,]
  text_eval<-paste0('data',i,'<-temp2')
  eval(parse(text=text_eval))
  print(i)
}

Count_table<-data.frame(contig_id=NA)
Length_table<-data.frame(contig_id=NA)
Species_table<-data.frame()

for(i in 1:(length(list_file))){
  print(list_file[i])
  text_eval<-paste0('dat<-data',i)
  eval(parse(text=text_eval))  
  temp<-as.list(by(dat,INDICES = dat$contig_name,function(x) r_length_mine(x)))
  temp<-sapply(temp,function(x) x[1])
  temp_len<-as.list(by(dat,INDICES = dat$contig_name,function(x) duplicate_lenremover(x)))
  temp_len<-sapply(temp_len,function(x) x[1])
  temp_dat<-data.frame(contig_names=names(unlist(temp)),contig_species=unname(unlist(temp)),stringsAsFactors = F)
  contig_id<-unname(sapply(temp_dat$contig_species,function(x) strsplit(x,'>| ')[[1]][2]))
  contig_species<-unname(sapply(temp_dat$contig_species,
                                function(x) paste(strsplit(x,'>| ')[[1]][-c(1,2)],collapse = ' ')))
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
  Count_table<-Count_table[!duplicated(Count_table),]
  Length_table<-Length_table[!duplicated(Length_table),]
}
for(i in 1:ncol(Count_table)){
  Count_table[is.na(Count_table[,i]),i]<-0
  Length_table[is.na(Length_table[,i]),i]<-0
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


write.table(Count_table,'Result_table.txt',quote=F,sep='\t')
write.table(Length_table,'Length_table.txt',quote=F,sep='\t',row.names=F)
rm(list=ls()[ls()!='Count_table'])

#name_f<-unname(sapply(unname(sapply(data1$full_name,function(x) strsplit(x,' ')[[1]][1])),
#                      function(x) substr(x,2,nchar(x))[1]))
#name_g<-unname(sapply(data1$contig_simplspecies,function(x) strsplit(x,' ')[[1]][1]))
#sum(name_f!=name_g)