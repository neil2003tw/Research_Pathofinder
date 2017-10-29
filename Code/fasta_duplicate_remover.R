file_result<-file('../virus_add_sim.fa')
blast_result<-readLines(file_result)
close(file_result)
all_index2result<-1:length(blast_result)
header_index2result<-grep('>',blast_result)
header<-blast_result[header_index]
duplicate_index2header<-which(duplicated(header))
for(i in 1:length(duplicate_index2header)){
  duplicated_fasta_index<-(header_index2result[duplicate_index2header[i]]):(header_index2result[duplicate_index2header[i]+1]-1)
  print(duplicated_fasta_index)
  all_index2result[duplicated_fasta_index]<-NA
}
blast_result_clean<-blast_result[!is.na(all_index2result)]
