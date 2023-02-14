rm(list = ls())
options(stringsAsFactors = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
# options(BioC_mirror="https://anaconda.org/bioconda/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

getwd()
setwd('C:/Users/LWH/Downloads/bioinfo/COVID&IS')
source("C:/Users/LWH/Downloads/bioinfo/AD/AD_functions.R")
source("PPI.R")
library(dplyr)
library(plyr)
library(tidyr)
library(R.utils)
library(biomaRt)
library(RcisTarget)
IRG<-read.table("immune/GeneList.txt",header=T,sep="\t")

gene<-intersect(gene_id,IRG$ID)
geneLists <- list(geneListName=fd_COV[fd_COV$entrezgene_id%in%
                                        gene,]$hgnc_symbol)

data(motifAnnotations_hgnc)
motifRankings <- importRankings("../AD_single/01_getdata/hg19-tss-centered-10kb-7species.mc9nr.feather")
motifEnrichmentTable_wGenes <- cisTarget(geneSets = geneLists, motifRankings,
                                         motifAnnot=motifAnnotations_hgnc)
motifs_AUC <- calcAUC(geneLists, motifRankings)
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, 
                                           motifAnnot=motifAnnotations_hgnc)
motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable, 
                                                   geneSets=geneLists,
                                                   rankings=motifRankings, 
                                                   nCores=3,
                                                   method="aprox")
motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable)
resultsSubset <- motifEnrichmentTable_wGenes_wLogo[1:5,]
library(DT)
datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=5))

anotatedTfs <- lapply(split(motifEnrichmentTable_wGenes$TF_highConf,
                            motifEnrichmentTable$geneSet),
                      function(x) {
                        genes <- gsub(" \\(.*\\). ", "; ", x, fixed=FALSE)
                        genesSplit <- unique(unlist(strsplit(genes, "; ")))
                        return(genesSplit)
                      })

get_tf<-function(k, motifEnrichmentTable){
  c = Reduce(intersect, list(k[order(k$NES,decreasing = T),]$motif[1:10],
                             k[order(k$AUC,decreasing = T),]$motif[1:10],
                             k[order(k$rankAtMax,decreasing = T),]$motif[1:15]))
  # c=c(c,"taipale_tf_pairs__GCM1_NHLH1_NCAGCTGNNNNNNNNTRCGGG_CAP_repr",
  #     "taipale_cyt_meth__SPIB_RAWWGRGGAAGTN_FL_meth")
  k = k[k$motif %in% c,]
  signifMotifNames <- motifEnrichmentTable %>% .[.$TF_highConf != "",] %>% 
    .$motif %>% .[.%in% c]
  incidenceMatrix <- getSignificantGenes(geneLists$geneListName, 
                                         motifRankings,
                                         signifRankingNames=signifMotifNames,
                                         plotCurve=TRUE, maxRank=5000, 
                                         genesFormat="incidMatrix",
                                         method="aprox", nCores = 3)$incidMatrix
  tf_list = k$TF_highConf %>% gsub("\\(.*\\)\\.", "", ., perl=T) %>% gsub(" ","",.)  %>% 
    str_split(., ";")
  tf_Matrix=matrix(ncol=ncol(incidenceMatrix))
  for (i in 1:length(tf_list)) {
    for (j in 1:length(tf_list[[i]])) {
      print(i)
      tf_Matrix=rbind(tf_Matrix,incidenceMatrix[i,])
    }
  }
  kf = cbind(matrix(tf_list %>% unlist()), tf_Matrix[-1,]) %>% as.data.frame()
  kf = aggregate.data.frame(kf[,2:ncol(kf)] %>% apply(2,as.numeric),
                       list(ID=kf$V1),sum)
  rownames(kf)<-kf$ID
  TF<-kf[,-1] %>% mutate_all(as.numeric) %>% as.matrix()
  return(TF)
}

k = motifEnrichmentTable_wGenes %>% .[.$TF_highConf != "",]
TF <- get_tf(k, motifEnrichmentTable)
deg_IS<-differ(ad_IS,pd_IS$group)
deg_IS[rownames(TF),]
deg_COV2[rownames(TF),]
TF<-TF[-1,]

library(reshape2)
library(STRINGdb)
library(magrittr)
edges=melt(TF) %>% .[.$value!=0,]
colnames(edges) <- c("from","to","value")

string_db <- STRINGdb$new( version="11.5", species=9606,
                           score_threshold=400,
                           input_directory="../AD/01_getdata")

get_network = function(edges,string_db,fd=fd_COV) {
  motifs <- unique(as.character(edges[,1]))
  genes <- unique(as.character(edges[,2]))
  ENTREZID = fd[genes,]$entrezgene_id
  genes_df <- data.frame(ENTREZID = ENTREZID, SYMBOL = genes)
  data_mapped <- genes_df %>% 
    string_db$map(my_data_frame_id_col_names = "ENTREZID",
                  removeUnmappedRows = TRUE)
  genes_input <- get_input(data_mapped) %>% mutate(value=round(combined_score)) %>% 
    dplyr::select(c(1,2,4))
  colnames(genes_input) <- c("from","to","value")
  tmp <- rbind(edges,genes_input) %>% as.data.frame() %>% 
    graph.data.frame(directed = T)
  tmp1 <- rbind(edges,genes_input) %>% as.data.frame() %>% 
    graph.data.frame(directed = F)
  group <- data.frame(row.names = c(motifs,genes),
                      group=c(rep("TF", length(motifs)),rep("DEG", length(genes))),
                      image=c(rep("04_PPI/gene-mutation.png", length(motifs)),
                              rep("04_PPI/genes1.png", length(genes))))
  igraph::V(tmp)$deg <- igraph::degree(tmp) # 每个节点连接的节点数
  # igraph::V(tmp)$size <- igraph::degree(tmp)/5 #
  igraph::E(tmp)$width <- igraph::E(tmp)$value
  igraph::V(tmp)$grp <- igraph::V(tmp) %>% names() %>% group[.,] %>%.$group
  igraph::V(tmp)$image <- igraph::V(tmp) %>% names() %>% group[.,] %>% .$image
  igraph::E(tmp1)$width <- igraph::E(tmp)$value
  igraph::V(tmp1)$grp <- igraph::V(tmp) %>% names() %>% group[.,] %>%.$group
  igraph::V(tmp1)$image <- igraph::V(tmp) %>% names() %>% group[.,] %>% .$image
  return(list(one=tmp,two=tmp1))
}

net <- get_network(edges,string_db,fd_COV) %>% .$one
net1 <- get_network(edges,string_db,fd_COV)%>% .$two
edges1<-get.data.frame(net1)
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
values<-data.frame(row.names=V(net1)%>%names(),value=degree(net1))
nodes1 <- data.frame(id=c(motifs, genes),   
                    label=c(motifs, genes),    
                    title=c(motifs, genes), # tooltip 
                    shape=c(rep("diamond", length(motifs)),
                            rep("elypse", length(genes))),
                    color=c(rep("purple", length(motifs)), 
                            rep("skyblue", length(genes))),
                    value=values[c(motifs, genes),])


grouplist<-lapply(pd_COV2$Group, function(x)ifelse(x=="Control","Control","Case")) %>% 
  unlist()
data_IS<-ad_IS[motifs,]
hub_COV_box = immunebox(ad_COV2[motifs,], grouplist)
hub_IS_box = immunebox(data_IS, pd_IS$group)
source("flat_violin.R")
library(ggpubr)

hub_COV_box1<-hub_COV_box
hub_COV_box1$value<-sign(hub_COV_box1$value)*(abs(hub_COV_box1$value))/10000
distributions1 <- 
  ggplot(data = hub_COV_box1, 
         aes(x = Attribute, y = value, fill = Group)) +
  geom_flat_violin(position=position_nudge(x=0.2,y=0),alpha = 0.8,trim=T) +
  geom_point(aes(y = value, color = Group), 
             position = position_jitter(width = 0.15), size = 1, alpha = 0.1) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  stat_compare_means(aes(group=Attribute),label="p.signif",label.x=1.5,
                     size=5,hjust=1)+
  labs(y = "Expression Level", x = NULL) +
  guides(fill = 'none', color = 'none') +
  scale_y_continuous(limits = c(0, 3)) + #
  scale_fill_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C",
                               "#00FF66FF","#00FFFFFF")) +
  scale_colour_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C",
                                 "#00FF66FF","#00FFFFFF")) +
  scale_x_discrete(breaks=c("Control","Case"),labels=c("Control","COVID"))+
  facet_wrap(vars(Group),nrow=2,scales="free_y")+coord_flip()+
  theme_niwot()+theme(strip.text=element_text(size=12,face="bold.italic"))

hub_IS_box1<-hub_IS_box
k<-aggregate(hub_IS_box$value,by=list(Group=hub_IS_box$Group),min)
for (i in 1:nrow(k)) {
  hub_IS_box1[hub_IS_box1$Group==k$Group[i],]$value<-
    (hub_IS_box1[hub_IS_box1$Group==k$Group[i],]$value-k$x[i])
}

distributions2 <- 
  ggplot(data = hub_IS_box1, 
         aes(x = Attribute, y = value, fill = Group)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.8) +
  geom_point(aes(y = value, color = Group), 
             position = position_jitter(width = 0.15), size = 1, alpha = 0.1) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  stat_compare_means(aes(group=Attribute),label="p.signif",label.x=1.5,
                     size=5,hjust=1)+
  labs(y = "Expression Level", x = NULL) +
  guides(fill = 'none', color = 'none') +
  scale_y_continuous(limits = c(0, 4)) + #
  scale_fill_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C",
                               "#00FF66FF","#00FFFFFF")) +
  scale_colour_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C",
                                 "#00FF66FF","#00FFFFFF")) +
  scale_x_discrete(breaks=c("Control","Case"),labels=c("Control","IS"))+
  facet_wrap(vars(Group),nrow=2,scales="free_y")+coord_flip()+
  theme_niwot()+theme(strip.text=element_text(size=12,face="bold.italic"))

