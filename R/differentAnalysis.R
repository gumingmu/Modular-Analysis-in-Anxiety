# title: differentAnalysis.R
# author: gumingmu
# date: 2024-03-31
# descrption: This script is for different analysis, Enrichment analysis


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
exp_protein.f <- "D:/work/2.honghao_wgcna/metadata/matrix/table_PCG_all_two_batch_20230419/normalized_CPM_edgeR.txt"
count_protein.f <- "D:/work/2.honghao_wgcna/metadata/matrix/table_PCG_all_two_batch_20230419/cutoff200_PCG_symbol.txt"
group.f <- "D:/work/2.honghao_wgcna/metadata/anno/group_anno1.txt"
rename <- dplyr::rename
select <- dplyr::select
mutate <- dplyr::mutate
filter <- dplyr::filter

dir.create(dir)



# loading ------------------------
exp.d <- fread("D:/work/2.honghao_wgcna/metadata/matrix/table_PCG_all_two_batch_20230419/normalized_CPM_edgeR_symbol.txt")
group.anno <- fread("D:/work/2.honghao_wgcna/metadata/anno/group_anno1.txt") %>% 
  mutate(status = factor(status, levels = c("Control", "Susceptible", "Resceptible"))
         , Region = factor(Region, levels = c("PFC", "DG", "CA1", "CA3", "LS", "LH", "Hb")))

name <- unique(group.anno$Region)
file <- paste0("D:/work/2.honghao_wgcna/3.de/1.Region8/", name, "/Susceptible/deAllGene.txt")
gene.union <- lapply(file, function(x){
  fread(x, header = F) %>% unlist()
}) %>% unlist() %>% unique()

file_Male <- paste0("D:/work/2.honghao_wgcna/3.de/1.Region8/", name, "/Susceptible/", "Male","/deAllGene.txt")
gene_Male.union <- lapply(file_Male, function(x){
  fread(x, header = F) %>% unlist()
}) %>% unlist() %>% unique()

file_Female <- paste0("D:/work/2.honghao_wgcna/3.de/1.Region8/",name, "/Susceptible/", "Female","/deAllGene.txt")
gene_Female.union <- lapply(file_Female, function(x){
  fread(x, header = F) %>% unlist()
}) %>% unlist() %>% unique()

gene_sex.union <- unique(c(gene_Male.union, gene_Female.union))
#Resceptible
file_Male_res <- paste0("D:/work/2.honghao_wgcna/3.de/1.Region8/", name, "/Resceptible/", "Male","/deAllGene.txt")
gene_Male_res.union <- lapply(file_Male_res, function(x){
  fread(x, header = F) %>% unlist()
}) %>% unlist() %>% unique()

file_Female_res <- paste0("D:/work/2.honghao_wgcna/3.de/1.Region8/",name, "/Resceptible/", "Female","/deAllGene.txt")
gene_Female_res.union <- lapply(file_Female_res, function(x){
  fread(x, header = F) %>% unlist()
}) %>% unlist() %>% unique()

gene_sex_res.union <- unique(c(gene_Male_res.union, gene_Female_res.union))


#after get merge de data
gene_change.df <- fread(paste0(outdir, "pipedata/deMerge.txt")) %>% 
  rename("symbol" = "V1")






# function--------------------------------
#1. edgeR-------------
this.edgeR <- function(data, group, outdir, p = "p", cutOff = 0.05, lfc = 0){
  temp <- data %>% 
    column_to_rownames("symbol")
  
  exp.dge <- DGEList(counts = temp, group = group)
  exp.dge <- calcNormFactors(exp.dge)
  design <- model.matrix(~exp.dge$samples$group)
  temp.dge <- estimateDisp(exp.dge, design, robust = T)
  temp.fit <- glmFit(temp.dge)
  temp.lrt <- glmLRT(temp.fit)
  temp.de <- topTags(temp.lrt,n=nrow(temp.lrt$table))$table
  
  if(p == "p"){
    temp.de <- temp.de %>% 
      mutate(sign = ifelse(PValue < cutOff&abs(logFC) > lfc, "Y", "N"), change = ifelse(logFC > 0, ifelse(logFC > 1, "2x Up", "Up"), ifelse(logFC < -1, "2x Down", "Down"))
             , change1 = ifelse(PValue < cutOff&logFC > 0, "Up", ifelse(PValue < cutOff&logFC < 0, "Down", "NS")))
  }else if(p == "fdr"){
    temp.de <- temp.de %>% 
      mutate(sign = ifelse(FDR < cutOff&abs(logFC) > lfc, "Y", "N"), change = ifelse(logFC > 0, ifelse(logFC > 1, "2x Up", "Up"), ifelse(logFC < -1, "2x Down", "Down"))
             , change1 = ifelse(PValue < cutOff&logFC > 0, "Up", ifelse(PValue < cutOff&logFC < 0, "Down", "NS")))
  }
  
  temp.summ <- table(temp.de[, c(6, 8)])
  
  write.table(temp.de, paste0(outdir, "de.txt"), sep = "\t", quote = F)
  write.table(temp.summ, paste0(outdir, "de_summary.txt"), sep = "\t", quote = F)
  
  up_gene <- rownames(temp.de)[temp.de$sign == "Y"&temp.de$logFC > 0]
  down_gene <- rownames(temp.de)[temp.de$sign == "Y"&temp.de$logFC < 0]
  all_gene <- rownames(temp.de)[temp.de$sign == "Y"]
  
  fwrite(list(up_gene),  paste0(outdir, "deUpGene.txt"))
  fwrite(list(down_gene),  paste0(outdir, "deDownGene.txt"))
  fwrite(list(all_gene),  paste0(outdir, "deAllGene.txt"))
}

#2. merge de data-------------
this.deSum <- function(file, name, outfile){
  de_sum.l <- lapply(1:length(file), function(x){
    sum.df <- fread(file[x])
    up <- sum.df[2, 4] %>% unlist()
    down <- sum.df[2, 2] %>% unlist()
    
    data.frame(num = c(up, down),Region = rep(name[x], 2), type = c("Up", "Down"))
  })
  de_sum.df <- Reduce(rbind, de_sum.l)
  fwrite(de_sum.df, outfile, sep = "\t")
  
  return(de_sum.df)
}









# analysis--------------------------
outdir1 <- paste0(outdir, "edgeR/region/")
dir.create(outdir1, recursive = T)
library(edgeR)

#loading 
exp.d <- fread(count_protein.f)
group.anno <- fread(group.f)


# de
library(edgeR)

#Stress  Male
region.c <- unique(group.anno$Region)
sex.c <- "Male"
lapply(region.c, function(x){
  index <- which(group.anno$status != "Resceptible"&group.anno$Region == x&group.anno$sex == sex.c)+1
  temp.df <- exp.d[, c(1, index)]
  
  outdir1 <- paste0(outdir, x, "/Susceptible/", sex.c, "/")
  dir.create(outdir1)
  this.edgeR(temp.df, group = group.anno$status[index-1], outdir = outdir1, lfc = log2(1.3))
})
#Stress  Female
region.c <- unique(group.anno$Region)
sex.c <- "Female"
lapply(region.c, function(x){
  index <- which(group.anno$status != "Resceptible"&group.anno$Region == x&group.anno$sex == sex.c)+1
  temp.df <- exp.d[, c(1, index)]
  
  outdir1 <- paste0(outdir, x, "/Susceptible/", sex.c, "/")
  dir.create(outdir1)
  this.edgeR(temp.df, group = group.anno$status[index-1], outdir = outdir1, lfc = log2(1.3))
})

#Resceptible  Male
region.c <- unique(group.anno$Region)
sex.c <- "Male"
lapply(region.c, function(x){
  index <- which(group.anno$status != "Susceptible"&group.anno$Region == x&group.anno$sex == sex.c)+1
  temp.df <- exp.d[, c(1, index)]
  
  outdir1 <- paste0(outdir, x, "/Resceptible/", sex.c, "/")
  dir.create(outdir1)
  this.edgeR(temp.df, group = group.anno$status[index-1], outdir = outdir1, lfc = log2(1.3))
})
#Resceptible  Female
region.c <- unique(group.anno$Region)
sex.c <- "Female"
lapply(region.c, function(x){
  index <- which(group.anno$status != "Susceptible"&group.anno$Region == x&group.anno$sex == sex.c)+1
  temp.df <- exp.d[, c(1, index)]
  
  outdir1 <- paste0(outdir, x, "/Resceptible/", sex.c, "/")
  dir.create(outdir1)
  this.edgeR(temp.df, group = group.anno$status[index-1], outdir = outdir1, lfc = log2(1.3))
})



# plot--------------------------------------------------------
#Susceptible
name <- unique(group.anno$Region)
file <- paste0(outdir, name, "/Susceptible/de_summary.txt")
sus_sum.df <- this.deSum(file, name, outfile = paste0(outdir, "susceptible_sum.txt"))
#Resceptible
name <- unique(group.anno$Region)
file <- paste0(outdir, name, "/Resceptible/de_summary.txt")
res_sum.df <- this.deSum(file, name, outfile = paste0(outdir, "resceptible_sum.txt"))

all.df <- rbind(sus_sum.df, res_sum.df) %>% 
  mutate(status = rep(c("Susceptible", "Resceptible"), each = nrow(sus_sum.df))) %>% 
  mutate(Region = factor(Region, levels = name)) %>% 
  mutate(status = factor(status, levels = c("Susceptible", "Resceptible")))
fwrite(all.df, paste0(outdir, "allDe_sum.txt"), sep = "\t")

#plot
this.dePlot(all.df, outfile = paste0(outdir, "deSum.pdf"))
