blastREF="/NFSShare/neil/Pathogen/Workstage/blast_db/nr/"
echo Reference Path = $REF
INPUTPREFIX=$1
echo Input Path = $INPUTPREFIX

SAMPLES=($(find $INPUTPREFIX -type d -name 'SRR*'))
ELEMENTS=${#SAMPLES[@]}

for (( i=0;i<$ELEMENTS;i++)); do
   echo ${SAMPLES[$i]}
   R1FILES=($(find ${SAMPLES[$i]} -type f -iname '*_1.fastq'))
   R2FILES=($(find ${SAMPLES[$i]} -type f -iname '*_2.fastq'))
   NUMPAIRS=${#R2FILES[@]}
   echo Read 1 is $R1FILES
   echo Read 2 is $R2FILES
   echo cd ${SAMPLES[$i]}
   cd ${SAMPLES[$i]}
   echo Start tophat ${SAMPLES[$i]}
   tophat2 -p 5 --no-coverage-search -G /data/iGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf /data/iGenome/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome $R1FILES $R2FILES

   echo Coverting bam2sam2fasta
   samtools view tophat_out/unmapped.bam > tophat_out/unmapped.sam
   cat tophat_out/unmapped.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > tophat_out/unmapped.fastq

   echo Trinity
   /data/NGSTools/trinityrnaseq_r20131110/Trinity.pl --JM 200G --seqType fq --single tophat_out/unmapped.fastq --CPU 5

   echo BLAST
   cd $blastREF
   blastn -query ${SAMPLES[$i]}/trinity_out_dir/Trinity.fasta -out ${SAMPLES[$i]}/trinity_out_dir/nr_result.out -db nt -num_descriptions 50 -num_alignments 30 -num_threads 5
done

cd $INPUTPREFIX
Rscript /NFSShare/neil/Pathogen/Codechunk/blast_result_categorize.R
