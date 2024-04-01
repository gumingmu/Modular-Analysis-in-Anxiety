# title: cellType.R
# author: gumingmu
# date: 2024-03-31
# descrption: This script is for cellType analysis
rm(list = ls())
gc()
options(stringsAsFactors = F)

library(tidyverse)
library(data.table)

load("D:/source/color/color_single.Rdata")
source("D:/source/R function/analysis_function.R")
source("D:/source/R function/plot_function.R")
source("D:/work/SLE+circRNA/bioanalysis/Rdata/.R/function/function.R")

dir <- "D:/work/2.honghao_wgcna/3.de/1.Region8/"
outdir <- "D:/work/2.honghao_wgcna/3.de/1.Region8/"

library(readxl)
library(aplot)

outdir1 <- paste0(outdir, "4.articlePlot/")
dir.create(outdir1)

sex.c <- c("Male", "Female")
region.c <- group.anno$Region %>% unique()
status.c <- c("Susceptible", "Resceptible")
change.c <- c("Up", "Down")

anno_cell.df <- read_xlsx("D:/work/2.honghao_wgcna/articleData/geneList/BrainCellType.xlsx") %>% 
  as.list()
allGene.c <- fread("D:/work/2.honghao_wgcna/3.de/1.Region8/deResceptibleMerge.txt")$V1 %>% unlist()
cell.c <- names(anno_cell.df)

fish.df <- lapply(cell.c, function(x){
  gene.cell <- anno_cell.df[[x]] %>% na.omit()
  
  lapply(sex.c, function(y){
    lapply(region.c, function(b){
      lapply(status.c, function(z){
        if(z == "Resceptible"){
          label.1 <- "Resilient"
        }else{
          label.1 <- "Susceptible"
        }
        lapply(change.c, function(a){
          gene.sign <- fread(paste0("D:/work/2.honghao_wgcna/3.de/1.Region8/",b,"/", z, "/", y, "/", "/de", a, "Gene.txt"), header = F) %>% unlist()
          
          gene1 <- intersect(allGene.c, gene.cell)
          gene.inter <- intersect(gene.sign, gene1)
          
          num1.c <- length(gene1)
          num2.c <- length(gene.sign)
          num3.c <- length(gene.inter)
          num4.c <- length(allGene.c)
          
          p.c <- fisher.test(matrix(c(num3.c , num1.c - num3.c,num2.c - num3.c, num4.c - num1.c - num2.c + num3.c), nrow = 2), alternative = "greater")$p.value
          
          data.frame(cell = x, sex = y, region = b, status = label.1, change = a, p = p.c)
        }) %>% Reduce(rbind, .)
      }) %>% Reduce(rbind, .)
    }) %>% Reduce(rbind, .)
  }) %>% Reduce(rbind, .)
}) %>% Reduce(rbind, .) %>% 
  mutate(fdr = p.adjust(p, method = "BH"))

fwrite(fish.df, paste0(outdir1, "fisher.txt"), sep = "\t")

plot.df <- fish.df %>% 
  mutate(sex = factor(sex, levels = sex.c)
         , region = factor(region, levels = region.c)
         , status = factor(status, levels = c("Susceptible", "Resilient"))
         , change = factor(change, levels = change.c)) %>% 
  arrange(sex, region, status, change, cell) %>% 
  mutate(rank = rep(1:(nrow(fish.df)/5), each = 5)) %>% 
  mutate(signp = ifelse(p < 0.05, "*", "")
         , signfdr = ifelse(fdr < 0.05, "*", ""))

p1 <- ggplot(plot.df, aes(x = rank, y = cell))+
  geom_tile(aes(fill = -log10(p), color = -log10(p)))+
  scale_color_gradient(low = "white", high = color_mutli[3])+
  scale_fill_gradient(low = "white", high = color_mutli[3])+
  scale_x_continuous(expand = c(0, 0))+
  theme_void()+
  theme(text = element_text(size = 6), plot.margin = unit(rep(0, 4), "cm"), panel.background = element_rect(color = "grey60"), axis.text.y.left = element_text(size = 6, color = "black", vjust = 1, hjust = 1),legend.key.size =unit(rep(0.3, 2), "cm"))

p1.1 <- ggplot(plot.df, aes(x = rank, y = cell))+
  geom_tile(aes(fill = -log10(p), color = -log10(p)))+
  geom_text(aes(label = signp))+
  scale_x_continuous(expand = c(0, 0))+
  scale_color_gradient(low = "white", high = color_mutli[3])+
  scale_fill_gradient(low = "white", high = color_mutli[3])+
  theme_void()+
  theme(text = element_text(size = 6), plot.margin = unit(rep(0, 4), "cm"), panel.background = element_rect(color = "grey60"), axis.text.y.left = element_text(size = 6, color = "black", vjust = 1, hjust = 1),legend.key.size =unit(rep(0.3, 2), "cm"))

p1.2 <- ggplot(plot.df, aes(x = rank, y = cell))+
  geom_tile(aes(fill = -log10(fdr), color = -log10(fdr)))+
  geom_text(aes(label = signfdr))+
  scale_x_continuous(expand = c(0, 0))+
  scale_color_gradient(low = "white", high = color_mutli[3])+
  scale_fill_gradient(low = "white", high = color_mutli[3])+
  theme_void()+
  theme(text = element_text(size = 6), plot.margin = unit(rep(0, 4), "cm"), panel.background = element_rect(color = "grey60"), axis.text.y.left = element_text(size = 6, color = "black", vjust = 1, hjust = 1),legend.key.size =unit(rep(0.3, 2), "cm"))

p2 <- ggplot(plot.df, aes(x = rank, y = "Sex"))+
  geom_tile(aes(fill = sex, color = sex))+
  scale_y_discrete(expand = c(0, 0))+
  scale_color_manual(values = c("#4A708B", "#EE6AA7"))+
  scale_fill_manual(values = c("#4A708B", "#EE6AA7"))+
  theme_void()+
  theme(text = element_text(size = 6), plot.margin = unit(rep(0, 4), "cm"), axis.text.y.left = element_text(size = 6, color = "black"),legend.key.size =unit(rep(0.3, 2), "cm"))
p3 <- ggplot(plot.df, aes(x = rank, y = "Region"))+
  geom_tile(aes(fill = region, color = region))+
  scale_y_discrete(expand = c(0, 0))+
  scale_color_manual(values = color_mutli[6:12])+
  scale_fill_manual(values = color_mutli[6:12])+
  guides(color = F)+
  theme_void()+
  theme(text = element_text(size = 6), plot.margin = unit(rep(0, 4), "cm"), axis.text.y.left = element_text(size = 6, color = "black"),legend.key.size =unit(rep(0.3, 2), "cm"))
p4 <- ggplot(plot.df, aes(x = rank, y = "Phenotype"))+
  geom_tile(aes(fill = status, color = status))+
  scale_y_discrete(expand = c(0, 0))+
  scale_color_manual(values = alpha(c("#FF8C00", "#6495ED"), 0.1))+
  scale_fill_manual(name = "Phenotype",values = alpha(c("#FF8C00", "#6495ED"), 0.1))+
  guides(color = F)+
  theme_void()+
  theme(text = element_text(size = 6), plot.margin = unit(rep(0, 4), "cm"), axis.text.y.left = element_text(size = 6, color = "black"),legend.key.size =unit(rep(0.3, 2), "cm"))
p5 <- ggplot(plot.df, aes(x = rank, y = "Direction"))+
  geom_tile(aes(fill = change, color = change))+
  scale_y_discrete(expand = c(0, 0))+
  scale_color_manual(values = alpha(color_mutli[3:2], 0.1))+
  scale_fill_manual(name = "Phenotype",values = alpha(color_mutli[3:2], 0.1))+
  guides(color = F)+
  theme_void()+
  theme(text = element_text(size = 6), plot.margin = unit(rep(0, 4), "cm"), axis.text.y.left = element_text(size = 6, color = "black"),legend.key.size =unit(rep(0.3, 2), "cm"))

p <- insert_bottom(p1.2, p2, height = 0.1)
p <- insert_bottom(p, p3, height = 0.1)
p <- insert_bottom(p, p4, height = 0.1)
p <- insert_bottom(p, p5, height = 0.1)
p

ggsave(paste0(outdir1, "cellTypeSignFdr.pdf"), width = 9*2.5, height = 4*2.5, units = "cm")
