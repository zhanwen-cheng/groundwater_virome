####AMG calculated by vibrant summary####
getwd()
setwd("/media/chengzw_bk/revision/AMG/")

virus_kofam <- read.delim("./VIBRANT_annotations_tsne_modify.tsv")
virus_kofam <- read.delim("./vibrant_result.txt",header = F)
colnames(virus_kofam) <- c("protein","scaffold","KO","AMG",
                           "KO name","KO evalue","KO score","KO v-score",
                           "Pfam","Pfam name","Pfam evalue","Pfam score","Pfam v-score",
                           "VOG","VOG name","VOG evalue","VOG score","VOG v-score")
virus_kofam <- virus_kofam[virus_kofam$AMG=='AMG',]
virus_kofam <- virus_kofam[grepl('ych_',virus_kofam$protein),]
virus_kofam$protein <- gsub('ych_','',virus_kofam$protein)
virus_kofam$scaffold <- gsub('ych_','',virus_kofam$scaffold)
virus_kofam$protein <- gsub('_fragment_1','',virus_kofam$protein)
virus_kofam$scaffold <- gsub('_fragment_1','',virus_kofam$scaffold)
virus_kofam$protein <- gsub('_fragment_2','',virus_kofam$protein)
virus_kofam$scaffold <- gsub('_fragment_2','',virus_kofam$scaffold)
virus_kofam$protein <- gsub('_fragment_.*','',virus_kofam$protein)
virus_kofam$scaffold <- gsub('_fragment_.*','',virus_kofam$scaffold)
length(unique(virus_kofam$scaffold))

virus_kofam <- virus_kofam[,c(1:4)]

proviruse <- read.delim("./provirus.list",header = F)
virus_kofam <- virus_kofam[!virus_kofam$scaffold %in% proviruse$V1,]

kolist <- read.delim("./ko00001_20221203.txt",sep = '#')
kolist <- kolist[,c('level_b_des','level_c_des','level_d_num')]
kolist[kolist$level_b_des=='Protein families: genetic information processing',]$level_b_des='Others'
kolist[kolist$level_b_des=='Protein families: signaling and cellular processes',]$level_b_des='Others'
kolist[kolist$level_b_des=='Others',]$level_c_des='Others'
kolist <- kolist[!duplicated(kolist),]
virus_kofam <- merge(virus_kofam,kolist,by.x = 'KO',by.y = 'level_d_num')

virus_kofam1 <- virus_kofam[,c(1,2,5)]
virus_kofam1 <- virus_kofam1[!duplicated(virus_kofam1),]
virus_kofam2 <- virus_kofam[,c(1,2,6)]
virus_kofam2 <- virus_kofam2[!duplicated(virus_kofam2),]

amg_plot1 <- as.data.frame(table(virus_kofam1$level_b_des))
amg_plot1 <- amg_plot1[order(amg_plot1$Freq,decreasing = T),]
amg_plot1 <- amg_plot1[c(2:12,1),]
rownames(amg_plot1) <- c(1:12)
amg_plot1$Var1 <- factor(amg_plot1$Var1,levels = rev(amg_plot1$Var1))

library(ggplot2)
library(cowplot)
library(dplyr)

ggplot(amg_plot1,aes(x=Var1,y=Freq))+
  geom_bar(stat = 'identity')+
  theme(text = element_text(family = "arial",colour = 'black'))+
  coord_flip()+
  scale_y_continuous(expand = c(0,0),limits = c(0,150))+
  theme(panel.grid = element_blank(),panel.border = element_blank())+
  theme(text = element_text(family = "arial",colour = 'black'))
#1000*300

amg_plot2 <- as.data.frame(table(virus_kofam2$level_c_des))
amg_plot2 <- amg_plot2[order(amg_plot2$Freq,decreasing = T),]
amg_plot2 <- amg_plot2[c(2:21),]
rownames(amg_plot2) <- c(1:20)
amg_plot2$Var1 <- factor(amg_plot2$Var1,levels = rev(amg_plot2$Var1))

ggplot(amg_plot2,aes(x=Var1,y=Freq))+
  geom_bar(stat = 'identity')+
  theme(text = element_text(family = "arial",colour = 'black'))+
  coord_flip()+
  scale_y_continuous(expand = c(0,0),limits = c(0,120))+
  theme(panel.grid = element_blank(),panel.border = element_blank())+
  theme(text = element_text(family = "arial",colour = 'black'))
#1000*300

####only plot AMG in vibrant prediction####
virus_kofam <- read.delim("./VIBRANT_annotations_tsne_modify.tsv")
virus_kofam <- virus_kofam[virus_kofam$AMG=='AMG',]
virus_kofam <- virus_kofam[grepl('ych_',virus_kofam$protein),]
virus_kofam$protein <- gsub('ych_','',virus_kofam$protein)
virus_kofam$scaffold <- gsub('ych_','',virus_kofam$scaffold)
virus_kofam$protein <- gsub('_fragment_1','',virus_kofam$protein)
virus_kofam$scaffold <- gsub('_fragment_1','',virus_kofam$scaffold)
virus_kofam$protein <- gsub('_fragment_2','',virus_kofam$protein)
virus_kofam$scaffold <- gsub('_fragment_2','',virus_kofam$scaffold)
virus_kofam$protein <- gsub('_fragment_.*','',virus_kofam$protein)
virus_kofam$scaffold <- gsub('_fragment_.*','',virus_kofam$scaffold)
length(unique(virus_kofam$scaffold))

virus_kofam <- virus_kofam[,c(1:4)]

proviruse <- read.delim("./provirus.list",header = F)
virus_kofam <- virus_kofam[!virus_kofam$scaffold %in% proviruse$V1,]

virus_kofam$site <- 'chen'
virus_kofam[grepl('_chen_',virus_kofam$protein),'site'] <- 'chen'
virus_kofam[grepl('_market_',virus_kofam$protein),'site'] <- 'market'
virus_kofam[grepl('_yca23_',virus_kofam$protein),'site'] <- 'yca23'
virus_kofam[grepl('_yca30_',virus_kofam$protein),'site'] <- 'yca30'
virus_kofam[grepl('_yca30_',virus_kofam$protein),'site'] <- 'yca30'

library(ggplot2)
library(cowplot)
library(dplyr)

virus_kofam_count <- virus_kofam[,c('KO','site')] %>%
  group_by(KO,site) %>%
  dplyr::summarize(gene_num = n())

kolist <- read.delim("./ko00001_20221203.txt",sep = '#')
kolist <- kolist[,c('level_d_num','level_d_des')]
kolist <- kolist[!duplicated(kolist),]
virus_kofam_count <- merge(virus_kofam_count,kolist,by.x = 'KO',by.y = 'level_d_num',all.x = T)
virus_kofam_count[virus_kofam_count$KO=='K21029',]$level_d_des='moeB'
virus_kofam_count[virus_kofam_count$KO=='K21140',]$level_d_des='mec'
virus_kofam_count$gene_name <- sapply(strsplit(as.character(virus_kofam_count$level_d_des),";"),"[[",1)
virus_kofam_count <- virus_kofam_count[order(virus_kofam_count$gene_num,decreasing = T),]
rownames(virus_kofam_count) <- c(1:nrow(virus_kofam_count))
virus_kofam_count$KO_re <- factor(virus_kofam_count$KO,levels = unique(virus_kofam_count$KO))
virus_kofam_count$gene_name_re <- factor(virus_kofam_count$gene_name,levels = unique(virus_kofam_count$gene_name))
top_20 <- unique(virus_kofam_count$KO)[1:30]
virus_kofam_count <- virus_kofam_count[virus_kofam_count$KO %in% top_20,]
virus_kofam_count$gene_num_re <- '>20'
virus_kofam_count[virus_kofam_count$gene_num <=20,]$gene_num_re <- '<20'
virus_kofam_count[virus_kofam_count$gene_num <= 10,]$gene_num_re <- '<10'
virus_kofam_count[virus_kofam_count$gene_num <= 5,]$gene_num_re <- '<5'

p_v_1 <- ggplot(virus_kofam_count, aes(y = site, x = KO_re,size = gene_num_re,color=site))+
  geom_point(shape = 19)+
  #scale_size(name = "gene number")+
  ylab("ko") +
  xlab("host class")+
  scale_size_manual(values = c('<5' =1,'<10' =2,'<20'=3,'>20'=4))+
  #scale_size_manual(values = c('1' =1,'2' =2,'4'=4))+
  scale_color_manual(values = c('chen'='#F8766D',
    'market'='#7CAE00',
    'yca23'='#00BFC4',
    'yca30'='#C77CFF'))+
  theme(panel.border = element_rect(fill = NA, color = 'black',size = 1, linetype = "solid"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = 'gray'),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 90,size =12),
        axis.text.y = element_text(size = 10,),
        axis.title = element_text(size = 14))+
  theme(legend.position = 'none')+
  theme(text = element_text(family = 'arial',color = 'black'))
p_v_1
p_v_2 <- p_v_1+ scale_x_discrete(position = 'top',breaks = c(virus_kofam_count$KO_re),
                                     labels = c(virus_kofam_count$gene_name_re))+
  theme(text = element_text(family = 'serif',color = 'black'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0,size =12),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14))+
  theme(text = element_text(family = 'arial',color = 'black'))
p_v_2
p_v <- ggdraw(insert_xaxis_grob(p_v_1,get_x_axis(p_v_2,position='top')))
p_v
#1500*600
#800*300

######################################################################################
#############################viral gene arrange plot##################################
######################################################################################
getwd()
setwd("/media/chengzw_bk/revision/AMG")

library(Biostrings)

vibrant <- read.delim("/media/chengzw_bk/revision/AMG/vibrant_result.txt",header = F)
vibrant <- vibrant[,c(1,2,3,4,5,9,10,14,15)]
# virus_kw <- c("capsid", "phage", "terminase", "base plate", "baseplate", 
#               "prohead", "virion", "virus", "viral", "tape measure", "tapemeasure neck", 
#               "tail", "head", "bacteriophage", "prophage", "portal", "DNA packaging", 
#               "T4", "p22", "holin")
viral_vib <- read.delim("/media/chengzw_bk/revision/AMG/vibrant_result_viralgene.txt",header = F)

vibrant[vibrant$V1 %in% viral_vib$V1,'V16'] <- 'virus_gene'
AMG_virus <- unique(vibrant[vibrant$V4=='AMG',]$V2)
vibrant <- vibrant[vibrant$V2 %in% AMG_virus,]
rownames(vibrant) <- c(1:nrow(vibrant))

inter_v <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(inter_v) <- c('V1','V2','V3','V4','V5','V9','V10','V14','V15','V16')

k=5 #k means how far should AMG be included in viral genes
for (i in rownames(vibrant[vibrant$V4=='AMG',])) {
  i=as.numeric(i)
  inter_v2=vibrant[c((i-k):i),]
  inter_v3=vibrant[c(i:(i+k)),]
  inter_v4=vibrant[c((i-k):(i+k)),]
  print(i)
  if ((!isEmpty(table(inter_v2$V16))) &
      (!isEmpty(table(inter_v3$V16))) &
      (length(unique(inter_v4$V2))==1)
      ) {
    inter_v = rbind(inter_v,inter_v4)
  }
  inter_v=inter_v[!duplicated(inter_v),]
}
inter_v <- inter_v[,c(1,2,3,4,10,5,6,7,8,9)]

AMG_virus2 <- unique(inter_v[inter_v$V4=='AMG',]$V2)
inter_v <- inter_v[inter_v$V2 %in% AMG_virus2,]
inter_v$V1 <- gsub("_fragment_1","",inter_v$V1)
inter_v$V2 <- gsub("_fragment_1","",inter_v$V2)
inter_v$V1 <- gsub("_fragment_2","",inter_v$V1)
inter_v$V2 <- gsub("_fragment_2","",inter_v$V2)
inter_v$V2 <- gsub("_length_.*","",inter_v$V2)

gff <- read.delim("/media/chengzw_bk/vibrant/VIBRANT_tsne/gene2.gff",header = F,sep = ' ')
colnames(gff) <- c('protein','start','end','direction')
inter_v <- merge(inter_v,gff,by.x='V1',by.y='protein',all.x=T)
inter_v$type <- 'Others'
inter_v[!is.na(inter_v$V16),'type'] <- 'virus_gene'
inter_v[inter_v$V4=='AMG','type'] <- 'AMG'
inter_v <- inter_v[order(inter_v$V2,inter_v$start),]
write.csv(inter_v,file='AMG_5viral_gene.csv')


library(gggenes)
library(ggplot2)
inter_v[inter_v$direction==-1,'direction'] = 0
inter_v$V2 <- gsub("ych_","",inter_v$V2)

ggplot(inter_v[c(200:210),], aes(xmin = start, xmax = end, y = 'V2', 
                              fill = type, forward=direction)) +
  #geom_gene_arrow()+
  geom_gene_arrow(arrowhead_width = grid::unit(6, "mm"),
                  arrowhead_height = grid::unit(6, "mm"),
                  arrow_body_height = grid::unit(4, "mm"))+
  facet_wrap(~ V2,scales = 'free', ncol = 1)+
  theme(axis.text.y = element_text(size = 0.01),
        legend.position = 'none')+
  #theme_genes()
  scale_fill_manual(values = c("AMG"="#3887BC",
                               "Others"="#e5e5e5",
                               "virus_gene"="#E52829"))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        text = element_text(family = "arial"))
#800*600