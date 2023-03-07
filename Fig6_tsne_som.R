####tsne+som analysis####
getwd()
setwd("/media/chengzw_bk/revision/tsne_som/")


library(mmgenome)
library(Rtsne)
library(dplyr)
library(kohonen)
library(vegan)
library(ggplot2)
library(ggrepel)

#################################
####claculate 5-mer frequence####
#################################
set.seed(100)
assembly<-readDNAStringSet("/media/chengzw_bk/revision/vcontact2_mix/merge.fasta", format = 'fasta')
out<- cbind.data.frame(names(assembly),
                       width(assembly),
                       letterFrequency(assembly, letters = c("GC"), as.prob = T)*100)
colnames(out)<-c("scaffold","length","gc")
out$scaffold<-as.character(out$scaffold)

out$color[grep("GOV2.0_",out$scaffold)]<-"GOV2.0"
out$color[grep("datong_",out$scaffold)]<-"GW"
out$color[grep("glen_",out$scaffold)]<-"GW"
out$color[grep("us_",out$scaffold)]<-"GW"
out$color[grep("wwtp_",out$scaffold)]<-"WWTP"
out$color[grep("ych_",out$scaffold)]<-"GW"

kmer_forC <- oligonucleotideFrequency(assembly, 5, as.prob = T)
kmer_revC <- oligonucleotideFrequency(reverseComplement(assembly),5, as.prob = T)
kmer <- (kmer_forC + kmer_revC)/2 * 100
rownames(kmer) <- out$scaffold

######################
####tsne calculate####
######################
print("Be patient calculating BH-SNE")
tsne_out <- Rtsne(as.matrix(kmer),check_duplicates = F)
tsne_plot <- data.frame(x_tsne = tsne_out$Y[, 1], y_tsne = tsne_out$Y[,2])
res <- cbind.data.frame(scaffold = as.character(names(assembly)),tsne_plot)
out_all_plot <- merge(x = out, y = res, by = "scaffold")

cbbPalette <- c("FE8308","3F8EAA","DDD399","E52823")

tsne_label <- data.frame(a1=c(40,0,-30),
                         a2=c(0,-20,30),
                         a3=c('WWTP','GW','GOV2.0'))
#tsne plot
ggplot()+
  geom_point(data = out_all_plot,
             aes(x = x_tsne, y = y_tsne,color=color),
             size=0.8)+
  #theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_color_manual(values = c("GOV2.0"="#3F8EAA",
                                "GW"="#E52829",
                                "WWTP"="#DDD399"))+
  geom_label_repel(data = tsne_label,aes(x=a1,y=a2,color=a3),label=tsne_label$a3,family='arial')+
  theme(text = element_text(family = 'arial',color = 'black'))+
  theme(axis.title = element_text(size = 20))+
  theme(axis.text = element_text(size=15))+
  theme(legend.position = 'none')
####550*550

##########################
####som training model####
##########################
grid.size <- ceiling((dim(kmer)[1]) ^ (1/2.5))
#grid.size-=36
#som.grid <- somgrid(xdim = 10, ydim = 10, topo = 'hexagonal', toroidal = T)
som.grid <- somgrid(xdim = 36, ydim = 36, topo = 'hexagonal', toroidal = T)
set.seed(100)
som.model <- som(data.matrix(kmer), grid = som.grid)
#som.model.dim_2.5 <- som.model
plot(som.model,type="changes")####鍒涘缓浜唖om.grid缁村害澶珮浜嗭紝璇漷raining progress鏃朵笉鏀舵暃, 33*33涓嶆敹鏁?
#600*400 
mydata <- as.matrix(as.data.frame(som.model$codes))
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:50) wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
par(mar=c(6,5,5,3))
plot(1:50, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares", main="Within cluster sum of squares (WCSS)")
#600*400
som.events <- som.model$codes[[1]]
set.seed(100)
#som_cluster_k <- kmeans(som.events, centers = 30, iter.max = 100, nstart = 10)$cluster
som_cluster_k <- kmeans(som.events, centers = 20, iter.max = 100, nstart = 10)$cluster
# Colour palette definition
cluster_palette <- function(x, alpha = 0.6) {
  n = length(unique(x)) * 2
  rainbow(n, start=2/6, end=6/6, alpha=alpha)[seq(n,0,-2)]
}
cluster_palette_init = cluster_palette(som_cluster_k)
bgcol = cluster_palette_init[som_cluster_k]
#show the same plot with the codes instead of just colours
# plot(som.model, type="codes", 
#      bgcol = bgcol, 
#      main = " ", codeRendering="stars")
# add.cluster.boundaries(som.model, som_cluster_k)

contig2node <- cbind(rownames(som.model$data[[1]]),
                     som.model$unit.classif)
colnames(contig2node) <- c('contig','node')
contig2node <- as.data.frame(contig2node)
contig2node$node <- sub("^", "V", contig2node$node )

node_cluster <- as.data.frame(som_cluster_k)
colnames(node_cluster) <- 'group'
node_cluster$group <- sub("^","group_",node_cluster$group)
node_cluster$node <- rownames(node_cluster)
rownames(node_cluster) <- NULL
contig2node2group <- merge(contig2node,node_cluster, by='node')

####som plot group and node####
contig2node2group[grep("GOV2.0_",contig2node2group$contig),'site']<-"GOV2.0"
contig2node2group[grep("datong_",contig2node2group$contig),'site']<-"GW"
contig2node2group[grep("glen_",contig2node2group$contig),'site']<-"GW"
contig2node2group[grep("us_",contig2node2group$contig),'site']<-"GW"
contig2node2group[grep("wwtp_",contig2node2group$contig),'site']<-"WWTP"
#contig2node2group[grep("ych_",contig2node2group$contig),'site']<-"YCH"
contig2node2group[grep("ych_",contig2node2group$contig),'site']<-"GW"
test <- contig2node2group %>%
  group_by(group,site) %>%
  dplyr::summarise(count = n())

for (i in unique(test$group)) {
  test[test$group==i,'total_contig']= sum(test[test$group == i,'count'])
}
test$log_contig <- log(test$total_contig, base = 2)

test$group<-factor(test$group, levels = c('group_1','group_2','group_3','group_4','group_5',
                                          'group_6','group_7','group_8','group_9','group_10',
                                          'group_11','group_12','group_13','group_14','group_15',
                                          'group_16','group_17','group_18','group_19','group_20'))
test$site<-factor(test$site)

test$contig_num_mod <- cut(test$total_contig,breaks = c(0,50,100,200,400,600),
                           labels = c(2,3,4,5,6))
test$contig_num_mod <- as.numeric(test$contig_num_mod)

#som plot
ggplot(data = test,aes(x = contig_num_mod/2, y=count, group=site, fill=site, width = contig_num_mod))+
  geom_bar(width = 1,stat = "identity",position = "fill")+
  #scale_fill_manual(values = c("#89D0F5", "#CC0033","#00CC00","#FFFF00"))+
  #scale_fill_manual(values = c("#89D0F5","#CC0033","#00CC00")) +
  scale_fill_manual(values = c("#3F8EAA","#E52829","#DDD399")) +
  coord_polar("y",start = 0)+
  #facet_grid(.~ group)+
  facet_wrap(~ factor(test$group), ncol = 10)+
  #theme_void()+
  #ggtitle("5 mer 100 node and 30 group seed_100")+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x = '',y = '')+
  theme(text = element_text(family = 'arial',color = 'black'))+
  theme(legend.position = 'left')
#800*200


som_tsne <- merge(contig2node2group, out_all_plot,by.x = 'contig', by.y = 'scaffold')
hc.norm.cent = som_tsne %>% group_by(group) %>% select(x_tsne, 
                                                       y_tsne) %>% summarize_all(mean)
hc.norm.cent$label <- NA

hc.norm.cent[hc.norm.cent$group=='group_1','label'] = 'group1'
#hc.norm.cent[hc.norm.cent$group=='group_4','label'] = 'group4'
#hc.norm.cent[hc.norm.cent$group=='group_8','label'] = 'group8'
#hc.norm.cent[hc.norm.cent$group=='group_15','label'] = 'group15'
hc.norm.cent[hc.norm.cent$group=='group_19','label'] = 'group19'

hc.norm.cent[hc.norm.cent$group=='group_9','label'] = 'group9'
hc.norm.cent[hc.norm.cent$group=='group_18','label'] = 'group18'

hc.norm.cent[hc.norm.cent$group=='group_10','label'] = 'group10'
hc.norm.cent[hc.norm.cent$group=='group_17','label'] = 'group17'

#tsne + som group
ggplot(som_tsne, aes(x = x_tsne, y = y_tsne , colour = group))+
  geom_point(size=0.8) + 
  scale_colour_manual(values = c("group_1" = "#E52829","group_2" = "#B2B1B9",
                                 "group_3" = "#B2B1B9","group_4" = "#B2B1B9",
                                 "group_5" = "#B2B1B9", "group_6" = "#B2B1B9",
                                 "group_7" = "#B2B1B9","group_8" = "#B2B1B9",
                                 "group_9" = "#A6CEE3","group_10" = "#DDD399",
                                 "group_11" = "#B2B1B9","group_12" = "#B2B1B9",
                                 "group_13" = "#B2B1B9","group_14" = "#B2B1B9",
                                 "group_15" = "#B2B1B9","group_16" = "#B2B1B9",
                                 "group_17" = "#D09C56","group_18" = "#3F8EAA",
                                 "group_19" = "#EA247E","group_20" = "#B2B1B9"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_label_repel(aes(label = label), data = hc.norm.cent,family='arial') +
  guides(colour = FALSE)+
  theme(text = element_text(family = 'arial',color = 'black'))+
  theme(axis.title = element_text(size = 20))+
  theme(axis.text = element_text(size=15))+
  theme(legend.position = 'none')
#550*550
write.csv(contig2node2group,file = './contig2node2group.txt', col.names = T, row.names = F)

####plot heatmap####
contig2node2group <- read.delim("./contig2node2group.txt", sep = ',')
contig2node2group$node <- gsub("V",'node_',contig2node2group$node)
kofam <- read.delim("/media/chengzw_bk/revision/tsne_som/kofam_modify2_new.txt",sep = ';',quote = '', header = F)
#28415
colnames(kofam) <- c('protein','KO','KO_des')
kofam <- kofam[!grepl('YP_',kofam$protein),]
kofam <- kofam[!grepl('NP_',kofam$protein),]
#65153 protein-->KO relationship
g2g <- read.delim("/media/chengzw_bk/vcontact2/viral_genomes_g2g.csv", sep = ',')
colnames(g2g) <- c('protein','contig','keywords')
kofam <- merge(kofam,g2g,by = 'protein')
kofam$keywords <-NULL

som_kofam <- merge(kofam, contig2node2group, by = 'contig')
group_kofam <- group_by(som_kofam[c('group','KO')],group,KO)
group_kofam <- dplyr::summarise(group_kofam,count = n())
group_kofam_origin <- group_kofam
for (j in unique(group_kofam$KO)){
  group_kofam[group_kofam$KO == j,'KOhitttedgroup'] = dim(group_kofam[group_kofam$KO == j,])[1]
}
kohitgroup <- group_kofam[,c('KO','KOhitttedgroup')]
kohitgroup <- kohitgroup[!duplicated(kohitgroup),]
kohitgroup <- kohitgroup[order(kohitgroup$KOhitttedgroup,decreasing = T),]
#kohitgroup <- kohitgroup[c(1:1000),]#鍙栧墠1000涓熀鍥犲惂
kohitgroup <- kohitgroup[c(1:100),]
group_kofam <- group_kofam[(group_kofam$KO %in% kohitgroup$KO),]

group_kofam_g <- dcast(group_kofam,KO~group, value.var = 'count')
rownames(group_kofam_g) <- group_kofam_g$KO
group_kofam_g <- group_kofam_g[kohitgroup$KO,]
group_kofam_g$KO <- NULL
group_kofam_g[is.na(group_kofam_g)] = 0
group_kofam_g[group_kofam_g >= 1] =1
group_kofam_g <- group_kofam_g[,c('group_1','group_19','group_4',
                                  'group_9','group_18',
                                  'group_10','group_17',
                                  'group_2','group_3',
                                  'group_5','group_6',
                                  'group_7','group_8','group_11',
                                  'group_12','group_13','group_14','group_15',
                                  'group_16','group_20')]


####plot heatmap####
library('pheatmap')
pheatmap(group_kofam_g, cluster_cols = F,cluster_rows = FALSE,scale = 'none',
         show_rownames = F,show_colnames = T,legend = F,
         #color = colorRampPalette(c('white','red'))(100),
         color = c('#A6CEE3','#EA4A34'),
         cellwidth = 10, cellheight = 2,
         #annotation_row = ko_an,annotation_col = node_an,
         gaps_col = c(1:19),
         fontfamily='arial',angle_col = 45
)
#666*400



















