#summary the taxonomed v-contig and plot sankey in this script
getwd()
setwd("/media/chengzw_bk/revision/vcontact2/")

gbg <- read.delim("/media/chengzw_bk/revision/vcontact2/genome_by_genome_overview_exclude5kb_NODE.csv",
                  sep = ",",header = F)
gbg <- gbg[,c('V2','V7','V9')]
VC_tax <- read.delim("/media/chengzw_bk/revision/vcontact2/viral_cluster_overview_27taxedVC.csv",
                     sep = ",")
VC_tax <- VC_tax[,c(1,3,4,5,6,7,8)]

gbg <- merge(gbg,VC_tax,by.x='V9',by.y='VC',all.x=T)
dim(gbg[!is.na(gbg$Realm),])

library(networkD3)
library(dplyr)

links <- read.csv("/media/chengzw_bk/revision/5kb/vcontact2/sankey.csv")
links <- data.frame(links)
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source),
         as.character(links$target)) %>% unique()
)
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "contig", NodeID = "name",
                   sinksRight=FALSE, fontFamily = 'arial',
                   fontSize = 15)
p

library(htmlwidgets)
saveWidget(p, file="./sankey.html")
library(webshot)
webshot("./sankeyNetwork.html", "sankeyNetwork.pdf")
