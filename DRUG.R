rm(list = ls())
options(stringsAsFactors = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
# options(BioC_mirror="https://anaconda.org/bioconda/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

getwd()
setwd('C:/Users/LWH/Downloads/bioinfo/COVID&IS')
source("../AD/AD_functions.R")
source("differ.R")
library(tidyverse)
library(clusterProfiler) 
library(org.Hs.eg.db) 
library(STRINGdb)
library(igraph)
library(ggraph)
library(biomaRt)
library(ggpubr)
library(stringr)
library(forcats)
library(ggvenn)
library(ComplexUpset)
library(ggplot2)
library(ggVennDiagram)
load(file = "Rdatas/gene_id.Rdata")
load(file = "Drug.Rdata")
load(file = "Rdatas/AD_COVGSE198449.Rdata")
load(file = "C:/Users/LWH/Downloads/bioinfo/AD/00_drugbank/drug_link.Rdata")

uniprot = drug_targets_polypep_ex_ident %>% .[.$resource == "UniProtKB",]
rownames(uniprot) = uniprot$parent_key
drug_targets = drug_targets %>% .[.$id %in% uniprot$parent_key,]
ensembl = biomaRt::useDataset("hsapiens_gene_ensembl",
                              mart=biomaRt::useMart("ensembl"))
gene_targets = biomaRt::getBM(attributes=c('uniprotswissprot','hgnc_symbol'),
                              filters = "uniprotswissprot",
                              values = uniprot$identifier,
                              mart = ensembl)
gene_targets = gene_targets[!is.na(gene_targets$hgnc_symbol),]
DTGs = gene_targets$hgnc_symbol
x = list('DEGs_COVID' = DEGs_COVID,
         'DEGs_IS' = DEGs_IS,
         "DTGs" = DTGs)

venn.plot <- venn.diagram(
  x = x,
  filename = "3triple_Venn.tiff",
  col = "transparent",
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  label.col = c("darkred", "white", "darkblue", "white",
                "white", "white", "darkgreen"),
  cex = 2.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("darkred", "darkblue", "darkgreen"),
  cat.cex = 2.5,
  cat.fontfamily = "serif",
  cat.dist = c(0.06, 0.06, 0.03),
  cat.pos = 0
)

abc_mat = data.frame(row.names = union(union(DEGs_COVID, DEGs_IS), DTGs))
abc_data= data.frame(abc_mat, DEGs_COVID = rownames(abc_mat) %in% DEGs_COVID,
                     DEGs_IS = rownames(abc_mat) %in% DEGs_IS,
                     DTGs =  rownames(abc_mat) %in% DTGs)
data<-arrange_venn(abc_data,verbose=T,outwards_adjust=1.2,
                   starting_grid_size = 116,
                   max_iterations=200)
abc_venn = (
  ggplot(data)
  + coord_fixed()
  + theme_void()
  + scale_color_venn_mix(abc_data)
  + geom_curve(
    data=data['ORM1', ],
    aes(xend=x+0.01, yend=y+0.01), x=1.5, y=2.5, curvature=.2,ncp=10,size =1.05,
    arrow = arrow(length = unit(0.025, "npc"))
  )
  + annotate(
    geom='text', x=0.8, y=2.6, size=8,
    label= "common DEGs with Drug Targets"
  )
)

(
  ggvenn =  abc_venn
  + geom_venn_region(data=abc_data, alpha=0.05)
  + geom_point(aes(x=x, y=y, color=region), size=1)
  + geom_venn_circle(abc_data)
  + geom_venn_label_set(abc_data, aes(label=region), outwards_adjust = 2, size =8)
  + geom_venn_label_region(
    abc_data, aes(label=size), 
    outwards_adjust=1.25, size =7,
    position=position_nudge(y=0.2)
  )
  + scale_fill_venn_mix(abc_data, guide='none')
  + guides(color="none")
  #+ theme(legend.title = element_text(size =16), legend.text = element_text(size =14))
)


rownames(proteins)<-proteins$identifier
drugs<-drug_targets[!duplicated(drug_targets$id,),];rownames(drugs)<-drugs$id
proteins$name<-drugs[proteins$parent_key,]$name
proteins1<-proteins[!is.na(proteins$name),]
proteins1<-proteins1[-13,]
rownames(proteins1)<-proteins1$identifier
DF$Protein<-proteins1[DF$Identifier,]$name
library(openxlsx)
write.xlsx(DF,"immune/DF.xlsx")

drug2 = drug1[drug1$Drug_score>80|drug1$Gene_number>1,]
drug2<-drug2[!drug2$SName%in%c("Zinc","Zinc sulfate, unspecified form",
                               "Copper"),]
drug2$Drug_score <- drug2$Drug_score %>% sqrt()
drug2 = drug2[order(drug2$Gene_number, drug2$Drug_score, decreasing = T),]
p <- ggplot(data = drug2, mapping = aes(
  x = fct_reorder(Name, Gene_number),
  y = Gene_number,
  fill = Drug_score))+
  labs(y = "Numbers of Target Genes", x = "Drug Name") + geom_col() + coord_flip()  +
  theme(axis.text = element_text(size = 16, face = "bold"), 
        axis.title = element_text(size = 18, face = "bold"), 
        axis.title.y = element_text(margin = margin(r = 15)), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14, face = "bold")) +
  scale_fill_distiller(palette = "Spectral")












