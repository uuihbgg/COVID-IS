rm(list = ls())
options(stringsAsFactors = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
# options(BioC_mirror="https://anaconda.org/bioconda/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

getwd()
setwd('C:/Users/LWH/Downloads/bioinfo/COVID&IS')

library(ggVennDiagram)
library(VennDiagram)
library(ggvenn)
library(ComplexUpset)
library(sva)
library(limma)
library(dplyr)
library(plyr)
library(tidyr)
library(org.Hs.eg.db)
library(customLayout)
library(EnhancedVolcano)
source("../AD/AD_functions.R")
library(multiMiR)

get_deggene<-function(ad,pd,a="Control",b="Asymptomatic"){
  pd<-pd[pd$Group%in%c(a,b),]
  group_list<-pd$Group %>% lapply(function(x)ifelse(x==a,"Control","Case")) %>% 
    unlist()
  deg<-differ(group_list,ad[,rownames(pd)])
  gene<-differ(group_list,ad[,rownames(pd)])%>%
    .[abs(.$logFC)>0.5&.$adj.P.Val<0.01,] %>% rownames()
  return(list(one=gene,two=deg))
}
get_b=function(x,n){
  list_b<-list()
  for (i in 1:n) {
    list_b[[i]]<- combn(x,i)
  }
  return(list_b)
}
differ = function(group_list, exp) {
  design <- model.matrix(~0+factor(group_list))#制作分组矩阵
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exp)
  
  cont.matrix=makeContrasts(contrasts=c('Case-Control'), levels = design)#制作差异比较矩阵
  
  fit <- lmFit(exp, design)#构建线性拟合模型
  
  fit2=contrasts.fit(fit, cont.matrix)#根据线性模型和比较矩阵进行差值运算
  
  fit2=eBayes(fit2)#贝叶斯检验
  
  tempOutput = topTable(fit2, coef='Case-Control', n=Inf)#生成结果报告
  
  deg = na.omit(tempOutput)
  return(deg)
}

load(file="Rdatas/COV_exp.Rdata")
load(file = "Rdatas/IS_exp.Rdata")

deg_IS<-differ(pd_IS$group,ad_IS)

V1 = EnhancedVolcano(deg_IS,
                     lab = rownames(deg_IS),
                     x = 'logFC',
                     y = 'adj.P.Val',
                     title = "Volcano plot of IS",
                     subtitle = NULL,
                     caption = NULL,
                     #selectLab = selectLab_italics,
                     xlab = bquote(~Log[2]~ 'fold change'),
                     ylab = bquote(~Log[10]~ italic('FDR')),
                     xlim = c(-1.5,2),
                     ylim = c(-.1,22),
                     pCutoff = 10e-3,
                     FCcutoff = 0.5,
                     pointSize = 3.0,
                     captionLabSize = 20,
                     labSize = 5.0,
                     labCol = 'black',
                     labFace = 'bold',
                     boxedLabels = TRUE,
                     parseLabels = TRUE,
                     col = c('black', 'pink', 'purple', 'red3'),
                     colAlpha = 0.7,
                     legendLabels = c("NS", expression(Log[2] ~ FC), "FDR", 
                                      expression(FDR ~ and ~ log[2] ~ FC)),
                     legendPosition = 'bottom',
                     legendLabSize = 18,
                     legendIconSize = 5.0,
                     drawConnectors = TRUE,
                     widthConnectors = 1.0,
                     colConnectors = 'black') + coord_flip()


V2 = EnhancedVolcano(deg_COV2,
                     lab = rownames(deg_COV2),
                     selectLab=genename,
                     x = 'logFC',
                     y = 'adj.P.Val',
                     title = "Volcano plot of COVID",
                     subtitle = NULL,
                     caption = NULL,
                     #selectLab = selectLab_italics,
                     xlab = bquote(~Log[2]~ 'fold change'),
                     ylab = bquote(~Log[10]~ italic('FDR')),
                     xlim = c(-2,2),
                     ylim = c(-.1,30),
                     pCutoff = 10e-3,
                     FCcutoff = 0.5,
                     pointSize = 3.0,
                     captionLabSize = 20,
                     labSize = 5.0,
                     labCol = 'black',
                     labFace = 'bold',
                     boxedLabels = TRUE,
                     parseLabels = TRUE,
                     col = c('black', 'pink', 'purple', 'red3'),
                     colAlpha = 0.7,
                     legendLabels = c("NS", expression(Log[2] ~ FC), "FDR", 
                                      expression(FDR ~ and ~ log[2] ~ FC)),
                     legendPosition = 'bottom',
                     legendLabSize = 18,
                     legendIconSize = 5.0,
                     drawConnectors = TRUE,
                     widthConnectors = 1.0,
                     colConnectors = 'black') + coord_flip()

gene_COV<-deg_COV2 %>% .[.$adj.P.Val<.01&abs(.$logFC)>0.5,] %>% rownames()
gene_IS<-deg_IS[abs(deg_IS$logFC)>0.5&deg_IS$adj.P.Val<0.01,]%>%rownames() 
genename = intersect(gene_IS, gene_COV)
x = list('DEGs_COVID' = DEGs_COVID,
         'DEGs_IS' = DEGs_IS)

venn <- Venn(x)
data <- process_data(venn)
ggvenn = ggplot() +
  geom_sf(aes(fill=id %>% as.character()), data = venn_region(data)) +
  geom_sf(size = 2, lty = "dashed", color = "grey", data = venn_setedge(data), show.legend = F) +
  geom_sf_text(aes(label = c("DEGs_COVID","DEGs_IS")), size =6, data = venn_setlabel(data)) +
  geom_sf_label(aes(label=count), fontface = "bold", 
                data = venn_region(data), label.size = 0.5, size = 5)+
  theme_void() +scale_fill_manual(bquote(bold("Group")), 
                                  labels= c("DEGs_COVID","common DEGs","DEGs_IS"),
                                  values = c(rgb(160,33,240,max = 255), 
                                             rgb(220,78,78,max = 255), 
                                             rgb(255,197,207,max = 255)))+
  theme(legend.text = element_text(size = 9, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"))
V1/V2/ggvenn


            
            
