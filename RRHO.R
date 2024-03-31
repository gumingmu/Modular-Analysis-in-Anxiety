# title: RRHO.R
# author: gumingmu
# date: 2024-03-31
# descrption: This script is for RRHO analysis
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



library(RRHO)
library(ggpubr)

outdir1 <- paste0(outdir, "3.RRHO/")
dir.create(outdir1)

region.c <- group.anno$Region %>% unique()
region.c <- region.c

## Susceptible---------
status.c <- "Susceptible"

outdir2 <- paste0(outdir1, status.c, "/")
dir.create(outdir2)
## Male--------------
sex.c <- "Male"

outdir3 <- paste0(outdir2, sex.c, "/")
dir.create(outdir3)

de_all.df <- fread(paste0(outdir, "deMergeData/",sex.c, "/de", status.c, "Merge.txt"))

compare.l <- combn(region.c, 2) %>% as.data.frame()

p.l <- lapply(1:ncol(compare.l), function(x){
  c1 <- compare.l[1, x]
  c2 <- compare.l[2, x]
  
  genel1 <- de_all.df %>% 
    select("V1"|paste0("logFC_", c1)|paste0("PValue_", c1)) %>% 
    rename_all(.funs = funs(c("V1", "logFC", "PValue"))) %>% 
    mutate(value = -log10(PValue)*logFC/abs(logFC)) %>%
    mutate(value = replace_na(value , 0)) %>% 
    select(V1,value)
  genel2 <- de_all.df %>% 
    select("V1"|paste0("logFC_", c2)|paste0("PValue_", c2)) %>% 
    rename_all(.funs = funs(c("V1", "logFC", "PValue"))) %>% 
    mutate(value = -log10(PValue)*logFC/abs(logFC)) %>% 
    mutate(value = replace_na(value , 0)) %>% 
    select(V1,value)
  
  rrho.r <- RRHO(genel1, genel2, BY = T, alternative = "enrichment", log10.ind = T)
  
  plot.rrhoHeat(rrho.r, label1 = c1, label2 = c2, scale = 0:1000, col = colorRampPalette(c(rep(c("#153e6a"), 1),"#32aecc","yellow","yellow","#f45b54","#f45b54",rep(c("#89111c"), 2)))(50))
})

ggarrange(plotlist = p.l, ncol = 7, nrow = 3, common.legend = T)
ggsave(paste0(outdir3, "allRegionPair.pdf"), width = 7*3, height =3.5*3, units = "cm")


## Female--------------
sex.c <- "Female"

outdir3 <- paste0(outdir2, sex.c, "/")
dir.create(outdir3)

de_all.df <- fread(paste0(outdir, "deMergeData/",sex.c, "/de", status.c, "Merge.txt")) 
compare.l <- combn(region.c, 2) %>% as.data.frame()

p.l <- lapply(1:ncol(compare.l), function(x){
  c1 <- compare.l[1, x]
  c2 <- compare.l[2, x]
  
  genel1 <- de_all.df %>% 
    select("V1"|paste0("logFC_", c1)|paste0("PValue_", c1)) %>% 
    rename_all(.funs = funs(c("V1", "logFC", "PValue"))) %>% 
    mutate(value = -log10(PValue)*logFC/abs(logFC)) %>%
    mutate(value = replace_na(value , 0)) %>% 
    select(V1,value)
  genel2 <- de_all.df %>% 
    select("V1"|paste0("logFC_", c2)|paste0("PValue_", c2)) %>% 
    rename_all(.funs = funs(c("V1", "logFC", "PValue"))) %>% 
    mutate(value = -log10(PValue)*logFC/abs(logFC)) %>% 
    mutate(value = replace_na(value , 0)) %>% 
    select(V1,value)
  
  rrho.r <- RRHO(genel1, genel2, BY = T, alternative = "enrichment", log10.ind = T)
  
  plot.rrhoHeat(rrho.r, label1 = c1, label2 = c2, scale = 0:1000, col = colorRampPalette(c(rep(c("#153e6a"), 1),"#32aecc","yellow","yellow","#f45b54","#f45b54",rep(c("#89111c"), 2)))(50))
})

ggarrange(plotlist = p.l, ncol = 7, nrow = 3, common.legend = T)
ggsave(paste0(outdir3, "allRegionPair.pdf"), width = 7*3, height =3.5*3, units = "cm")



## Resceptible---------
status.c <- "Resceptible"

outdir2 <- paste0(outdir1, status.c, "/")
dir.create(outdir2)
## Male--------------
sex.c <- "Male"

outdir3 <- paste0(outdir2, sex.c, "/")
dir.create(outdir3)

de_all.df <- fread(paste0(outdir, "deMergeData/",sex.c, "/de", status.c, "Merge.txt"))

compare.l <- combn(region.c, 2) %>% as.data.frame()

p.l <- lapply(1:ncol(compare.l), function(x){
  c1 <- compare.l[1, x]
  c2 <- compare.l[2, x]
  
  genel1 <- de_all.df %>% 
    select("V1"|paste0("logFC_", c1)|paste0("PValue_", c1)) %>% 
    rename_all(.funs = funs(c("V1", "logFC", "PValue"))) %>% 
    mutate(value = -log10(PValue)*logFC/abs(logFC)) %>%
    mutate(value = replace_na(value , 0)) %>% 
    select(V1,value)
  genel2 <- de_all.df %>% 
    select("V1"|paste0("logFC_", c2)|paste0("PValue_", c2)) %>% 
    rename_all(.funs = funs(c("V1", "logFC", "PValue"))) %>% 
    mutate(value = -log10(PValue)*logFC/abs(logFC)) %>% 
    mutate(value = replace_na(value , 0)) %>% 
    select(V1,value)
  
  rrho.r <- RRHO(genel1, genel2, BY = T, alternative = "enrichment", log10.ind = T)
  
  plot.rrhoHeat(rrho.r, label1 = c1, label2 = c2, scale = 0:1000, col = colorRampPalette(c(rep(c("#153e6a"), 1),"#32aecc","yellow","yellow","#f45b54","#f45b54",rep(c("#89111c"), 2)))(50))
})

ggarrange(plotlist = p.l, ncol = 7, nrow = 3, common.legend = T)
ggsave(paste0(outdir3, "allRegionPair.pdf"), width = 7*3, height =3.5*3, units = "cm")



## Female--------------
sex.c <- "Female"

outdir3 <- paste0(outdir2, sex.c, "/")
dir.create(outdir3)

de_all.df <- fread(paste0(outdir, "deMergeData/",sex.c, "/de", status.c, "Merge.txt"))
compare.l <- combn(region.c, 2) %>% as.data.frame()

p.l <- lapply(1:ncol(compare.l), function(x){
  c1 <- compare.l[1, x]
  c2 <- compare.l[2, x]
  
  genel1 <- de_all.df %>% 
    select("V1"|paste0("logFC_", c1)|paste0("PValue_", c1)) %>% 
    rename_all(.funs = funs(c("V1", "logFC", "PValue"))) %>% 
    mutate(value = -log10(PValue)*logFC/abs(logFC)) %>%
    mutate(value = replace_na(value , 0)) %>% 
    select(V1,value)
  genel2 <- de_all.df %>% 
    select("V1"|paste0("logFC_", c2)|paste0("PValue_", c2)) %>% 
    rename_all(.funs = funs(c("V1", "logFC", "PValue"))) %>% 
    mutate(value = -log10(PValue)*logFC/abs(logFC)) %>% 
    mutate(value = replace_na(value , 0)) %>% 
    select(V1,value)
  
  rrho.r <- RRHO(genel1, genel2, BY = T, alternative = "enrichment", log10.ind = T)
  
  plot.rrhoHeat(rrho.r, label1 = c1, label2 = c2, scale = 0:1000, col = colorRampPalette(c(rep(c("#153e6a"), 1),"#32aecc","yellow","yellow","#f45b54","#f45b54",rep(c("#89111c"), 2)))(50))
})

ggarrange(plotlist = p.l, ncol = 7, nrow = 3, common.legend = T)
ggsave(paste0(outdir3, "allRegionPair.pdf"), width = 7*3, height =3.5*3, units = "cm")








