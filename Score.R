rm(list=ls())

library(Seurat)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(stringr)
library(msigdbr)
library(RColorBrewer)

Treg <- readRDS("D:/R/0902/Treg.intergrated.rds")
setwd("D:/R/0902/Score")

df <- read.csv('Core.csv',sep = ',',check.names = FALSE)
features <- list(papermarker)

sce_score<- AddModuleScore(Treg,
                           features = features,
                           ctrl = 100,
                           name = "features")
head(sce_score@meta.data)
colnames(sce_score@meta.data)[12] <- 'Score'

mydata<- FetchData(sce_score,vars = c("umap_1","umap_2","Score"))
colours = c("#2166AC", "#4393C3","#D1E5F0", "#F7F7F7", "#F4A582", "#D6604D", "#B2182B", "#67001F")

##UMAP
ggplot(mydata,aes(x = umap_1,y =umap_2,colour = Score))+
  geom_point(size = 2)+
  scale_color_gradientn(values = seq(0,1,0.2),colours = colours )+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
        legend.position = "none")& NoAxes()
ggsave("Score.pdf",width = 4 , height = 4, dpi=600)


##Vlnplot
library(ggpubr)
vln.dat=FetchData(sce_score,c("Score","orig.ident","seurat_clusters"))
vln.dat$seurat_clusters <- factor(vln.dat$seurat_clusters, 
                                  levels = c("0", "1", "2"),
                                  labels = c("C0", "C1", "C2")) 

ggviolin(vln.dat, "seurat_clusters", "Score",
         fill = "seurat_clusters", 
         color = "seurat_clusters", 
         trim = T, 
         palette = c("#71ACD8","#F4DA90", "#E9A66E", "#60C2A6","#DC7656","#DCE49F"), 
         font.y = 1,  
         font.tickslab = c(20,"plain","black"), 
         add = "boxplot", 
         add.params = list(   
           fill = "white", 
           color = "black",  
           width = 0.2,  
           linetype = 1)) +
  theme(legend.position = "none")+
labs(x = NULL, 
     y = "Score",
     title = NULL)
ggsave("vln_score4.pdf",width = 4, height = 2, dpi=600)

