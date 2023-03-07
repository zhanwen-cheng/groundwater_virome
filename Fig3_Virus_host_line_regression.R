getwd()
setwd("/media/chengzw_bk/revision/host_predict")


library(ggplot2)
library(hrbrthemes)

require(devtools)
install_version("ggpmisc", version = "0.4.5")
library('ggpmisc')
library('ggplot2')

merge <- read.delim("./class_cov2vc_cov.txt")
merge['id_method'] <- 'merge'
hic <- read.delim("./class_cov2vc_cov_hic.txt")
hic['id_method'] <- 'hic based'
# seq_ho <- read.delim("./order_cov2vc_cov_seqhomology.txt")
# seq_ho['id_method'] <- 'sequence homology'
crispr <- read.delim("./class_cov2vc_cov_crisp.txt")
crispr['id_method'] <- 'crispr'
# #crispr['id_method'] <- 'blastn'
# crispr['id_method'] <- 'sequence homology'
blastn <- read.delim("./class_cov2vc_cov_blastn.txt")
blastn['id_method'] <- 'blastn'


line_plot <- rbind(merge,hic)
line_plot <- rbind(line_plot,blastn)
line_plot <- rbind(line_plot,crispr)

ggplot(line_plot,aes(x = vc_cov, y = class_cov, colour = class))+
  geom_point(size = 2)+
  scale_shape_manual(values = c(15,16,17,18))+
  stat_smooth(method='lm',formula = y~x,colour='red')+
  facet_wrap(~ id_method)


my.formula <- vc_cov ~ class_cov

####设置颜色
library(RColorBrewer)
#colourCount = length(unique(line_plot$order))
#getPalette = colorRampPalette(brewer.pal(12, "Paired"))
# cb_palette <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
#                 "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
#                 "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
#                 "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
#                 "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
#                 "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#处理后有73种差异还比较明显的颜色，基本够用
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

phylum2order <- read.delim("./phylum2order.txt", header = F,
                           col.names = c('phylum','order'))
class2order <- read.delim("./class2order.txt", header = F,
                           col.names = c('class','order'))

# line_plot <- merge(line_plot,class2order, by = 'class')
# line_plot <- merge(line_plot,phylum2order, by = 'order')

p <- ggplot(data = line_plot, aes(x = sqrt(vc_cov), y = sqrt(class_cov))) +
  geom_point(aes(color = class),size = 3) +
  #stat_smooth(method=lm, level=0.99)+
  geom_smooth(method = "lm", se=TRUE, color="red")+
  facet_wrap(~ id_method,scales = "free",ncol = 2)+
  scale_shape_manual(values = c(15,16,17,18))+
  stat_poly_eq(
    aes(label =paste(..eq.label.., ..adj.rr.label.., sep = "~~")),
    formula = y ~ x,  parse = TRUE,
    family="arial",
    size = 3.6,
    color="black",
    label.x = 0.1,  #0-1之间的比例确定位???
    label.y = 0.9)+
  theme(panel.background=element_rect(fill='transparent', color='black'),
        axis.text=element_text(color = 'black'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(legend.position = "bottom")+
  theme(legend.key.height = unit(10, "pt"),legend.key.width = unit(10, "pt"))+
  guides(color=guide_legend(nrow=4))+
  theme(legend.title = element_blank())+
  theme(text = element_text(family = 'arial',size=15),legend.text = element_text(face="italic"))+
#  theme(legend.position="bottom") +
#  guides(fill=guide_legend(nrow=4))+
  scale_color_manual(values = col_vector)
p
#1100*800

library(Hmisc)
merge_p = rcorr(as.matrix(line_plot[line_plot$id_method == 'merge',2:3]),type = 'spearman')
merge_p
#p = 0
merge_hic = rcorr(as.matrix(line_plot[line_plot$id_method == 'hic based',2:3]),type = 'spearman')
merge_hic
#p =  2e-04
merge_crispr = rcorr(as.matrix(line_plot[line_plot$id_method == 'crispr',2:3]),type = 'spearman')
merge_crispr
#NO P VALUE
merge_blastn = rcorr(as.matrix(line_plot[line_plot$id_method == 'blastn',2:3]),type = 'spearman')
merge_blastn
#p = 0.0171


ggplot(data = merge, aes(x = sqrt(vc_cov), y = sqrt(class_cov))) +
  geom_point(aes(color = class),size = 3) +
  #stat_smooth(method=lm, level=0.99)+
  geom_smooth(method = "lm", se=TRUE, color="red")+
  #facet_wrap(~ id_method,scales = "free",ncol = 2)+
  scale_shape_manual(values = c(15,16,17,18))+
  stat_poly_eq(
    aes(label =paste(..eq.label.., ..adj.rr.label.., sep = "~~")),
    formula = y ~ x,  parse = TRUE,
    family="arial",
    size = 6,
    color="black",
    label.x = 0.1,  #0-1之间的比例确定位???
    label.y = 0.9)+
  theme(panel.background=element_rect(fill='transparent', color='black'),
        axis.text=element_text(color = 'black'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(legend.position = "right")+
  theme(legend.key.height = unit(10, "pt"),legend.key.width = unit(10, "pt"))+
  guides(color=guide_legend(ncol = 1))+
  theme(legend.title = element_blank())+
  theme(text = element_text(family = 'arial',size=20),legend.text = element_text(face="italic"))+
  #  theme(legend.position="bottom") +
  #  guides(fill=guide_legend(nrow=4))+
  scale_color_manual(values = col_vector)
#900*600


cdhit <- read.delim("/media/chengzw_bk/revision/5kb/cdhit_0.95_0.8/v+m_5kb_cdhit.txt",header = F)
####
df <- read.delim("./v_2_h_merge_tax_final.txt")
df <- df[df$viruse %in% cdhit$V1,]
df <- df[(df$kingdom != 'Viruses'),]
df <- df[(df$VC != ''),]
df1 <- df[(df$class != ''),]
df1 <- df1[,c('class','VC')]
df1 <- df1[!duplicated(df1),]
vc_fre <- as.data.frame(table(df1$VC))
colnames(vc_fre) <- c('VC','VC_fre')
df1 <- merge(df1,vc_fre)
length(unique(df1$VC))
#95个VC
length(unique(df1[(df1$VC_fre >= 2),'VC']))
#44个VC???2个class


df2 <- df[(df$phylum != ''),]
df2 <- df2[,c('phylum','VC')]
df2 <- df2[!duplicated(df2),]
vc_fre <- as.data.frame(table(df2$VC))
colnames(vc_fre) <- c('VC','VC_fre')
df2 <- merge(df2,vc_fre)
length(unique(df2$VC))
#96个VC
length(unique(df2[(df2$VC_fre >= 2),'VC']))
#19个VC???2个class

VC_250 <- df[df$VC=='VC_250_0',]
VC_250 <- VC_250[VC_250$class!='',]
VC_250 <- VC_250[,c('VC','phylum','class')]
VC_250 <- VC_250[!duplicated(VC_250),]
length(unique(VC_250$phylum))
length(unique(VC_250$class))

####crispr summary####
cdhit <- read.delim("/media/chengzw_bk/revision/5kb/cdhit_0.95_0.8/v+m_5kb_cdhit.txt",header = F)
####
df <- read.delim("./v_2_h_merge_tax_final.txt")
df <- df[df$viruse %in% cdhit$V1,]
df_crispr <- df[df$crisp==1,]
df_hic <- df[df$hic==1,]
df_blastn <- df[df$blastn==1,]
length(unique(df_blastn$viruse))
length(unique(df_hic$viruse))
View(df[((df$hic==1) & (df$blastn==1)),])
dim(df[((df$hic==1) & (df$blastn==1)),])
dim(df[((df$hic==1) & (df$crisp==1)),])

df_VC <- df[df$VC!='',]
View(df_VC[df_VC$crisp==1,])
length(unique(df_VC[((df_VC$crisp==1) &(df_VC$class!='')),]$VC))
length(unique(df_VC[((df_VC$hic==1) &(df_VC$class!='')),]$VC))
length(unique(df_VC[((df_VC$blastn==1) &(df_VC$class!='')),]$VC))
length(unique(df_VC[(df_VC$class!=''),]$VC))
length(unique(df_VC[(df_VC$class!=''),]$class))
length(unique(df_VC[(df_VC$class!=''),]$class))
test1 <- df_VC[,c('class','VC')]
test1 <- test1[test1$class!='',]
test1 <- test1[!duplicated(test1),]
View(table(test1$VC))
test2 <- df_VC[,c('phylum','VC')]
test2 <- test2[test2$phylum!='',]
test2 <- test2[!duplicated(test2),]
View(table(test2$VC))



