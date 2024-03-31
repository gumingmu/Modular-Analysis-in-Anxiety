# title: aracne.R
# author: gumingmu
# date: 2024-03-31
# descrption: This script is for GSEA analysis
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


# aracne
outdir1 <- paste0(outdir, "1009sup/aracne/")
dir.create(outdir1)

library(bnlearn)

moduleStatus.df <- fread(paste0(outdir, "1009sup/10module.csv"))
module_gene1.df <- module_gene.df %>% 
  mutate(colour = ifelse(colour == "blue", "cyan", colour)) %>% 
  mutate(colour = ifelse(colour == "white", "blue", colour))

th.01scale <- function(x){
  (x-min(x))/(max(x) - min(x))
}

lapply(1:nrow(moduleStatus.df), function(x){
  module <- moduleStatus.df[x, 1] %>% unlist()
  sex <- moduleStatus.df[x, 2] %>% unlist()
  region <- moduleStatus.df[x, 4] %>% unlist()
  status <- moduleStatus.df[x, 3] %>% unlist()
  status1 <- ifelse(status == "Resilient", "Resceptible", status)
  gene <- module_gene1.df$symbol[module_gene1.df$colour == module] %>% unlist()
  
  
  index <- group.anno$sample[group.anno$sex == sex&group.anno$Region == region&group.anno$status %in% c("Control", status1)]
  index <- group.anno$sample
  data <- exp.d %>% 
    filter(symbol %in% gene) %>% 
    column_to_rownames("symbol") %>% 
    select(all_of(index)) %>% 
    t() %>% as.data.frame()
  
  aracne.r <- aracne(data, mi = "mi-g")
  aracne.r1 <- aracne.r$arcs %>% as.data.frame() %>% select(1) %>% table() %>% as.data.frame() %>% arrange(-Freq)
  
  de.df <- fread(paste0("D:/work/2.honghao_wgcna/3.de/1.Region8/", region, "/", status1, "/", sex, "/", "de.txt")) %>% 
    select(1, 2,5) %>% 
    rename("symbol" = "V1")
  sign.c <- fread(paste0("D:/work/2.honghao_wgcna/3.de/1.Region8/", region, "/", status1, "/", sex, "/", "deAllGene.txt"), header = F) %>% unlist()
  
  
  ## n-hop
  n <- 6
  for(i in 2:n){
    if(i == 2){
      temp.df <- aracne.r$arcs %>% as.data.frame() %>% 
        left_join(aracne.r$arcs %>% as.data.frame() %>% rename_all(.funs = list(~c("to", "a")))) %>% select(1, 3) %>% rename_all(.funs = list(~c("from", "to"))) %>% distinct()
      
      result1.df <- rbind(aracne.r$arcs %>% as.data.frame(), temp.df) %>% distinct()
    }else{
      temp.df <- temp.df %>% 
        left_join(aracne.r$arcs %>% as.data.frame() %>% rename_all(.funs = list(~c("to", "a")))) %>% select(1, 3) %>% rename_all(.funs = list(~c("from", "to"))) %>% distinct()
      result1.df <-  rbind(result1.df, temp.df) %>% distinct()
    }
  }
  
  result2.df <- result1.df %>% select(1) %>% table() %>% as.data.frame()
  
  result.df <- aracne.r1 %>% 
    rename_all(.funs = list(~c("symbol", "outDegree"))) %>% 
    left_join(result2.df %>% rename_all(.funs = list(~c("symbol", "6hopDegree")))) %>% 
    mutate(outhub = ifelse(symbol %in% sign.c&outDegree > mean(outDegree) + 2*sd(outDegree), "Y", "N")
           , hub = ifelse(symbol %in% sign.c&`6hopDegree` > mean(`6hopDegree`) + 1*sd(`6hopDegree`), "Y", "N")
           , chub = ifelse(hub == "Y"&outhub == "Y", "Y", "N")) %>% 
    arrange(desc(chub)) %>% 
    left_join(de.df)
  
  fwrite(result.df, paste0(outdir1, module, "-", region, "-", sex, "-", status, ".txt"), sep = "\t")
  print(table(result.df$chub))
})
