#################################################
########metagenomic reads control, annotation####
#################################################
#metagenomic reads control
nohup ~/software/fastp -i chen_m_SR_1.fq -o fastp/chen_m_SR_1_clean.fq -I chen_m_SR_2.fq -O fastp/chen_m_SR_2_clean.fq &
#singleM 
singlem pipe --singlem_packages ~/data/data_base/singlem_extra_packages/release1/4.40.2013_08_greengenes_97_otus.with_euks.spkg --forward chen_m_SR_1.fq --otu_table chen.otu_table.csv
#metagenomic reads assembling
nohup spades.py -t 144 -m 6000 --meta -1 ../fastp/chen_m_SR_1_clean.fq -2 ../fastp/chen_m_SR_2_clean.fq -o ./ &
#prodigal prediction
nohup prodigal -a protein_seq.faa -d nucleotide_seq.fasta -o gene_predicted.gff -s protein_scores.stat -f gff -i ../../chen_m_1kb.fasta &
#VPF alignment
nohup hmmsearch --tblout vpf.orfs.txt --cpu 5 --notextw ~/data/data_base/vpf/viral_protein_families.new.hmm ../../prodigal/chen/protein_seq.faa > vpf_hmm.txt &
PFAM alignment
nohup hmmsearch --tblout pfam.orfs.txt --cpu 5 --cut_nc --notextw ~/data/data_base/Pfam/Pfam-A.new.hmm ../../prodigal/chen/protein_seq.faa > pfam_hmm.txt &
#alignment to KEGG with diamond
nohup diamond blastx --db ~/data/data_base/kegg_without_virus/kegg_database_without_virus.dmnd --query ../../prodigal/chen/nucleotide_seq.fasta -o diamond.txt -p 20 --outfmt 6 &


#################################################
########virome reads control, annotation#########
#################################################
#porechop remove long-reads adpater
porechop -i chen_v_LR.fastq.gz -o chen_v_LR_clean.fastq.gz
#virome short-reads quality control
nohup ~/software/fastp -i chen_v_SR_1.fq -o chen_v_SR_1_clean.fq -I chen_v_SR_2.fq -O chen_v_SR_2_clean.fq &
#singlem
singlem pipe --singlem_packages ~/data/data_base/singlem_extra_packages/release1/4.40.2013_08_greengenes_97_otus.with_euks.spkg --sequence chen_v_LR.fastq 
singlem pipe --singlem_packages ~/data/data_base/singlem_extra_packages/release1/4.40.2013_08_greengenes_97_otus.with_euks.spkg --forward chen_v_SR_1.fq --otu_table chen.otu_table.csv
#hybrid assembling
nohup spades.py -t 80 --meta -1 ../SR_fastp/chen_v_SR_1_clean.fq -2 ../SR_fastp/chen_v_SR_2_clean.fq --nanopore ../LR_fastp/chen_v_LR_clean.fastq.gz -o . &
#prodigal prediction
nohup prodigal -a protein_seq.faa -d nucleotide_seq.fasta -o gene_predicted.gff -s protein_scores.stat -f gff -i ../../chen_v_1kb.fasta &
#VPF alignment
nohup hmmsearch --tblout vpf.orfs.txt --cpu 5 --notextw ~/data/data_base/vpf/viral_protein_families.new.hmm ../../prodigal/chen/protein_seq.faa > vpf_hmm.txt &
PFAM alignment
nohup hmmsearch --tblout pfam.orfs.txt --cpu 5 --cut_nc --notextw ~/data/data_base/Pfam/Pfam-A.new.hmm ../../prodigal/chen/protein_seq.faa > pfam_hmm.txt &
#alignment to KEGG with diamond
nohup diamond blastx --db ~/data/data_base/kegg_without_virus/kegg_database_without_virus.dmnd --query ../../prodigal/chen/nucleotide_seq.fasta -o diamond.txt -p 20 --outfmt 6 &

#################################################
############v-contig identification##############
#################################################
####virsorter
conda activate virsorter
wrapper_phage_contigs_sorter_iPlant.pl -f chen_v_1kb.fasta --db 1 --wdir chen --ncpu 40 --data-dir ~/data/data_base/virsorter-data
#keep category 1,2,4,5 and contig longer than 5 kb

####ORFpipe2 for phylogeny pipeline and taxonomy
bash ORFpipe2.sh -f input.fa -t 40 ####ORFpip2 is a unpublished pipeline use kranken and taxator for taxonomy annotation
#keep contig classified as virus and longer than 5 kb

####identify v-contigs from vHMM
grep -v '^#' ../../prodigal/chen/gene_predicted.gff > tmp.gff
grep '^>' ../../prodigal/chen/protein_seq.faa | sed 's/>//'  | cut -f1 -d ' ' > tmp.seq.fasta.name
paste tmp.gff tmp.seq.fasta.name > gene_predicted.gff.modify

####format adjust
cp ../diamond/diamond.txt ./
cp ../hmmsearch/*.orfs.txt ./
grep 'Viruse' ../../orf_tax/chen/scaffolds_1kb+.fasta.1kb+_taxa.tab > scaffolds_1kb+.fasta.1kb+_taxa.tab_Viruse
cp ../../orf_tax/chen/scaffolds_1kb+.fasta.1kb+_original.name .

cut -f 1,4,5,7,10 gene.predicted.gff.modify > gene.predicted.gff.modify2
grep -v "^#" pfam.orfs.txt | tr -s ' ' '\t' |  cut -f 1,4 > pfam.orfs.txt.modify
grep -v '^#' vpf.orfs.txt | tr -s ' ' '\t' | awk '$5<=1e-5' | cut -f 1,3 > vpf.orfs.txt.modify    

awk '$3>=70 && $4>=100' diamond.txt > diamond.txt.filter    
sort -u -k 1,1 diamond.txt.filter > diamond.txt.filter.uniq      

####R-script
getwd()
setwd("/home/disk4-8T/chengzhanwen/onefour_meta_10g/vir_id/chen")

rm(list = ls())

gff<-read.delim("gene.predicted.gff.modify2", header = F)
colnames(gff)<-c("contig","start","end","direction","gene.id")
dim(gff)

pfam<-read.delim("pfam.orfs.txt.modify", header = F)
colnames(pfam)<-c("gene.id","pfam") 

virus.pfam<-read.delim("/home/disk1/db/Pfam-A.hmm.acc.modify.virus", header = F)  
colnames(virus.pfam)<-c("pfam.name","pfam")

vfam<-read.delim("vpf.orfs.txt.modify",header = F)
colnames(vfam)<-c("gene.id","vfam")

ko<-read.delim("diamond.txt.filter.uniq", header = F)  
ko<-ko[,1:2]
colnames(ko)<-c("gene.id","ko.name")
dim(ko)

ko.lst<-read.delim("/home/disk1/db/KEGG.db.from.ZAN201612/genes/ko/ko_genes.list", header = F)
colnames(ko.lst)<-c("ko","ko.name")

ko2<-merge(ko,ko.lst,by = "ko.name")
dim(ko2)

ko2<-ko2[!duplicated(ko2$gene.id),]
dim(ko2)


dim(gff)
tmp<-merge(gff,ko2,by = "gene.id",all = T) 
tmp2<-merge(tmp,pfam,by = "gene.id", all = T)
dim(tmp2)
tmp3<-merge(tmp2,vfam,by = "gene.id", all = T)
dim(tmp3)

lookat1<-which(!is.na(tmp3$ko))
lookat2<-which(!is.na(tmp3$vfam))
lookat<-intersect(lookat1,lookat2)

tmp3$vfam[lookat]<-NA
# identify those genes have virus pfam annotaiton and change those pfam annotation to NA
lookat3<-which(tmp3$pfam %in% virus.pfam$pfam)
tmp3$pfam[lookat3]<-NA

# now genes with non-virus pfam annotations
lookat1<-which(!is.na(tmp3$pfam))

lookat2<-which(!is.na(tmp3$vfam)) # genes with vfam annotation
lookat<-intersect(lookat1,lookat2) # genes with both vfam and non-virus pfam annotation

tmp3$vfam[lookat]<-NA
df<-tmp3

df$contig.leng<-sapply(strsplit(as.character(df$contig),"_"),"[[",4)
df<-df[df$contig.leng>=5000,]
df.lst<-split(df,as.character(df$contig))

df.lst<-Filter(function(x) {!(all(is.na(x$ko)) & all(is.na(x$pfam)) & all(is.na(x$vfam)))},df.lst)
df.lst<-Filter(function(x) {any(!is.na(x$vfam))}, df.lst)

library(foreach)
library(doParallel)
cl<-makeCluster(20,outfile="")
registerDoParallel(cl)
pb <- txtProgressBar(min=1,max = length(df.lst), style = 3)

result<-list()
for(i in 1:length(df.lst)){
  x<-df.lst[[i]] 
  m<-length(unique(x$gene.id)) # total gene number  
  lookat<-which(!is.na(x$vfam))
  m.vfam<-length(unique(x$gene.id[lookat])) # nubmer of gene with vfam hit  
  lookat<-which(!is.na(x$ko))
  m.ko<-length(unique(x$gene.id[lookat])) # number of gene with ko hit  
  lookat<-which(!is.na(x$pfam))
  m.pfam<-length(unique(x$gene.id[lookat])) # number of gene with non-virus pfam hit  
  c1<-(m.vfam>=5 & m.ko/m<0.2 & m.pfam/m<0.4 )  
  c2<-m.vfam > m.pfam  
  c3<-m.vfam/m > 0.6  
  c<-all(c1,c2,c3)  
  result[[i]]<-c
  cat("finish",i,"out of ",length(df.lst),"\n")
}
stopCluster(cl)

lookat<-which(unlist(result))
vdf.lst<-df.lst[lookat]
length(vdf.lst)/length(df.lst) # 1 in 10,000 should be a normal value

vContig.name1<-names(vdf.lst)  ####only one contigs

tmp<-read.delim("scaffolds_1kb+.fasta.1kb+_original.name",header=F)
colnames(tmp)<-c("contig","ont")

taxa<-read.delim("scaffolds_1kb+.fasta.1kb+_taxa.tab_Viruse",header=F)
colnames(taxa)[1]<-"ont"

taxa<-merge(tmp,taxa,by="ont")

taxa$contig.leng<-sapply(strsplit(as.character(taxa$contig),"_"),"[[",4)
taxa$contig.leng<-as.numeric(taxa$contig.leng)
taxa<-taxa[taxa$contig.leng>=5000,] # only kept contig longer than 5kbp
vContig.name2<-taxa$contig

# combine the annotation method with phylogenetic method 
vContig.name<-unique(c(vContig.name1,as.character(vContig.name2)))
write.csv(vContig.name,file="vContig.name.csv")


