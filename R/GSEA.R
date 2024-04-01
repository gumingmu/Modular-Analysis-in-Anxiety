# title: GSEA.R
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

# function-------------------
analysis.goGsea <- function(logfc, genetype = "SYMBOL", max = 30, outdir = NULL, species = "human", p = "p.adjust", goterm = "BP", color = c("#153e6a","gray70","#89111c"), labelWrap = 40, point.size = 1, text.size = 6, pcut = 0.05, plot.scale = 1){
  require(clusterProfiler)
  require(org.Hs.eg.db)
  require(org.Mm.eg.db)
  
  if(species == "human"){
    lib = org.Hs.eg.db
  }else if(species == "mouse"){
    lib = org.Mm.eg.db
  }

  if(genetype != "ENTREZID"){
    temp.gene <- bitr(names(logfc), fromType = genetype,
                      toType = c("ENTREZID"),
                      OrgDb = lib, drop = F) %>% 
      distinct(UQ(rlang::sym(genetype)), .keep_all = T) %>% 
      select(2) %>% unlist()
    names(logfc) <- temp.gene
    logfc <- logfc[!is.na(names(logfc))]
  }

  logfc <- sort(logfc, decreasing = T)
  
  go.result <- lapply(goterm, function(x){
    go_each.result <- gseGO(geneList = logfc,
                            OrgDb = lib,
                            ont = x,
                            minGSSize = 20,
                            maxGSSize = 500,
                            pvalueCutoff = 1,
                            pAdjustMethod = "BH",
                            verbose = FALSE) %>% setReadable(., OrgDb = lib, keyType = "ENTREZID")
    
    go.df <- go_each.result@result
    fwrite(go.df, paste0(outdir, "/gseaGo", x, ".txt"), sep = "\t")
    

    if(p == "p.adjust"){
      xlabel <- "-log10(p.adjust)"
      colnames(go.df)[7] <- "p.type"
    }else if(p == "p"){
      xlabel <- "-log10(Pvalue)"
      colnames(go.df)[6] <- "p.type"
    }
    

    if(min(go.df$p.type) >= pcut){
      shape.type <- c(4)
    }else{
      shape.type <- c(18, 4)
    }
    
    input.p <- go.df %>% 
      dplyr::mutate(sign = ifelse(p.type < pcut, "Y", "N") %>% factor(., levels = c("Y", "N"))) %>% 
      arrange(p.type) %>% 
      dplyr::filter(row_number() %in% 1:max) %>% 
      dplyr::mutate(Description = factor(Description, levels = rev(Description)))
    ggplot(input.p, aes(x = -log10(p.type), y = Description))+
      geom_point(aes(shape = sign, size = `setSize`, color = NES))+
      scale_shape_manual(values = shape.type)+
      scale_y_discrete(labels = function(x)str_wrap(x, width = labelWrap))+
      scale_color_gradientn(colours =  color, breaks = c(-3,-2, -1, 0,1, 2, 3))+
      scale_size_continuous(range = c(2,3)*point.size)+
      labs(y = "", x = xlabel, title = paste0("GSEA_GO", x))+
      guides(shape = F)+
      theme_bw()+
      theme(text = element_text(size = text.size), axis.text.y = element_text(color = "black"), legend.key.size = unit(0.2, "cm"), legend.text = element_text(size = text.size))
    
    if(max >= 20){
      plot.s <- (max)/20*plot.scale
    }else{
      plot.s <- (max+45)/60*plot.scale
    }
    ggsave(paste0(outdir, "/gseaGo", x, max, ".pdf"), width = 6+(10*plot.s-plot.s*2)*0.5, height = 9*plot.s, units = "cm")
    
    return(go_each.result)
  })
  
}





# analysis
outdir1 <- paste0(outdir, "4.articlePlot/")
dir.create(outdir1)

name <- unique(group.anno$Region)

#Susceptible
file <- paste0(outdir, name, "/Susceptible/de.txt")

lapply(1:length(name), function(x){
  de.df <- fread(file[x]) %>%
    mutate(rrho = -log10(PValue)*logFC/abs(logFC)) %>% 
    dplyr::select(1, rrho)
  num <- de.df %>% dplyr::select(2) %>% unlist()
  logfc <- num
  names(logfc) <- de.df$V1
  outdir1 <- paste0(outdir, "rrhoGsea/", name[x], "/Susceptible/")
  dir.create(outdir1, recursive = T)
  
  analysis.goGsea(logfc = logfc, species = "mouse", max = 30, outdir = outdir1)
})


#Susceptible
file <- paste0(outdir, name, "/Resceptible/de.txt")

lapply(1:length(name), function(x){
  de.df <- fread(file[x]) %>%
    mutate(rrho = -log10(PValue)*logFC/abs(logFC)) %>% 
    dplyr::select(1, rrho)
  num <- de.df %>% dplyr::select(2) %>% unlist()
  logfc <- num
  names(logfc) <- de.df$V1
  outdir1 <- paste0(outdir, "rrhoGsea/", name[x], "/Resceptible/")
  dir.create(outdir1, recursive = T)
  
  analysis.goGsea(logfc = logfc, species = "mouse", max = 30, outdir = outdir1)
})







## Male-------------------
sex.c <- "Male"
#Susceptible
file <- paste0(outdir, name, "/Susceptible/", sex.c, "/de.txt")

backgene.c <- exp.d$symbol
lapply(1:length(name), function(x){
  de.df <- fread(file[x]) %>%
    mutate(rrho = -log10(PValue)*logFC/abs(logFC)) %>% 
    dplyr::select(1, rrho)
  num <- de.df %>% dplyr::select(2) %>% unlist()
  logfc <- num
  names(logfc) <- de.df$V1
  outdir1 <- paste0(outdir, "rrhoGsea/", name[x], "/Susceptible/", sex.c, "/")
  dir.create(outdir1, recursive = T)
  
  analysis.goGsea(logfc = logfc, species = "mouse", max = 30, outdir = outdir1)
})


#Resceptible
file <- paste0(outdir, name, "/Resceptible/", sex.c, "/de.txt")

backgene.c <- exp.d$symbol
lapply(1:length(name), function(x){
  de.df <- fread(file[x]) %>%
    mutate(rrho = -log10(PValue)*logFC/abs(logFC)) %>% 
    dplyr::select(1, rrho)
  num <- de.df %>% dplyr::select(2) %>% unlist()
  logfc <- num
  names(logfc) <- de.df$V1
  outdir1 <- paste0(outdir, "rrhoGsea/", name[x], "/Resceptible/", sex.c, "/")
  dir.create(outdir1, recursive = T)
  
  analysis.goGsea(logfc = logfc, species = "mouse", max = 30, outdir = outdir1)
})
