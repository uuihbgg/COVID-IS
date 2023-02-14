rm(list = ls())
options(stringsAsFactors = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
# options(BioC_mirror="https://anaconda.org/bioconda/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

getwd()
setwd('C:/Users/LWH/Downloads/bioinfo/COVID&IS')

library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(fgsea)
library(customLayout)
library(msigdbr)
library(GOplot)
library(stringr)
R.utils::setOption( "clusterProfiler.download.method",'auto' )

load(file = "Rdatas/path.Rdata")
load(file="Rdatas/COV_exp.Rdata")
res = read.csv("C:/Users/LWH/Downloads/bioinfo/验证数据集/volcano_res.csv",
               row.names = 1)
universe = fd_COV[fd_COV$hgnc_symbol%in%rownames(ad_COV2),]$entrezgene_id %>%
  unique()
  
GOplot=dotplot(ego,x="GeneRatio",showCategory =10,split="ONTOLOGY")+ 
  theme(legend.text = element_text(size=10),
        legend.title=element_text(size=10,face="bold"),
        axis.text.y = element_text(size=10,lineheight=0.8,face = "bold"),
        axis.title.x=element_text(size=10,face="bold"),
        axis.text.x = element_text(size=10), 
        strip.text=element_text(size=10,face="bold"),
        strip.background = element_rect(color="black", fill="#F5F5F5"),
        legend.key.size = unit(3,"mm"),
        legend.margin = unit(0.5,"mm"),
        plot.margin = margin(t = 0,  
                             r = 1,  
                             b = 0,  
                             l = 0,  
                             unit = "mm"))+ 
  scale_size_continuous(range=c(1,6))+
  facet_grid(ONTOLOGY~., scale='free_y', space = "free")+ 
  scale_color_continuous (name="FDR", low = "red", high ="blue", 
                          labels = function(col) format(col, scientific = TRUE)) + 
  scale_y_discrete(labels=function(x) x %>% str_to_sentence() %>% 
                     str_wrap(width=40),
                   expand = c(0,0.6))  
 
 KEGGplot <- barplot(kk, color = "pvalue",showCategory=10)+
  theme(legend.text=element_text(size=10),legend.title=element_text(size=10,face = "bold"),
        axis.text.y=element_text(size=10,face = "bold"), 
        axis.text.x=element_text(size=10),
        axis.title.x=element_text(size=10,face = "bold"))+ 
  scale_fill_continuous (name="FDR", low = "red", high ="blue", 
            labels = function(col) format(col, scientific = TRUE))+
  # theme(legend.position = "none")+
  theme(plot.margin = margin(t=0, r=0, b=0, l=1, unit = "mm"),
        legend.key.size = unit(3, "mm"),
        legend.margin= unit(0.5,"mm"))

Reactome = barplot(x, color = "pvalue",showCategory=10)+
  theme(legend.text=element_text(size=10),
        legend.title=element_text(size=10,face = "bold"),
        legend.key.size= unit(3, "mm") ,
        legend.margin= unit(0.5,"mm"),
        plot.margin = margin(t=0, r=0, b=0, l=1, unit="mm"),
        axis.text.y=element_text(size=10,face = "bold"), 
        axis.text.x=element_text(size=10),
        axis.title.x=element_text(size=10,face = "bold"))+ 
  scale_fill_continuous (name="FDR", low = "red", high ="blue", 
                labels = function(col) format(col, scientific = TRUE))

(GOplot | (KEGGplot/Reactome))+
  plot_layout(widths = c(3,2))+
  # plot_annotation(tag_levels = 'A')&
  theme(plot.tag = element_text(size=9,face="bold"))



fd<-fd_COV[fd_COV$entrezgene_id%in%gene_id&fd_COV$hgnc_symbol!="HLA-DQA2",]
rownames(fd)<-fd$entrezgene_id
genelist_DO=data.frame(ID=gene_id,
                logFC=deg_COV2[fd[as.character(gene_id),]$hgnc_symbol,]$logFC)
gene<-genelist_DO$ID
DO <- DOSE::enrichDO(
  gene = gene,
  pvalueCutoff = 0.05,
  universe = as.character(universe)
  # pAdjustMethod = "none"
)
# DO <- setReadable(DO, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
# t(DO@result[1:10,])
DO_result = DO@result
DO_result<-DO_result[DO_result$pvalue<0.05,]
DO@result = DO_result
do=data.frame(Category = "DO",ID = DO_result$ID,
              Term = DO_result$Description %>% str_to_title(), 
              Genes = gsub("/", ", ", DO_result$geneID), adj_pval = DO_result$pvalue)
circ <- circle_dat(do, genelist_DO)
circ$genes = fd[as.character(circ$genes),]$hgnc_symbol
rownames(genelist_DO) = genelist_DO$ID 
genelist_DO$ID = fd[rownames(genelist_DO),]$hgnc_symbol
termNum = 9	#
termNum=ifelse(nrow(do)<termNum,nrow(do),termNum)
geneNum = nrow(genelist_DO)
chord <- chord_dat(circ, genelist_DO[1:geneNum,], do$Term[1:termNum])
# colnames(chord)[9]="Parasitic Protozoa Infection"
DOplot= GOChord(chord, lfc.max = 5,lfc.min = -2,process.label = 10,
                space = 0.02,
                gene.order = 'logFC', 
                gene.space = 0.25,      
                gene.size = 4.5,
                lfc.col=c('firebrick3', 'white','royalblue3'))+
  theme(legend.text=element_text(size=12,face="bold",family="serif"),
        legend.position = "bottom",legend.direction = "vertical",
        legend.text.align = 0, legend.box.just = "top", 
        legend.title=element_text(size=13,face="bold",family="serif"),
        plot.title=element_text(face = "bold")) + 
  scale_fill_continuous(breaks = seq(-2,5,by=2), limits=c(-2,5), 
                        low = "blue", high ="red")

DOplot$guides$size$title<-"DO Terms"
DOplot$guides$size$ncol<-3
# file.create("fig")
save_as_pdf({print(DOplot)},
            file.name = "DO.pdf",
            width = 9,
            height = 10.5)

