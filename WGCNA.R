# title: WGCNA.R
# author: gumingmu
# date: 2024-03-31
# descrption: This script is for WGCNA analysis

rm(list = ls())
dir <- "./wgcna0.2/"
outdir <- dir
dir.create(outdir, recursive = T)

cut = 0.2
pick = F
tom = F

library(tidyverse)
library(data.table)
library(WGCNA)
allowWGCNAThreads()
options(stringsAsFactors = FALSE);

exp.d <- fread("./metadata/normalized_CPM_edgeR_symbol.txt")
group.anno <- fread("./metadata/group_anno1.txt")
exp.df <- exp.d %>%
  column_to_rownames("symbol")

gene_wgc = t(exp.df)
#check
gsg <- goodSamplesGenes(gene_wgc)
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(gene_wgc)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(gene_wgc)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  gene_wgc = gene_wgc[gsg$goodSamples, gsg$goodGenes]
}

if(pick == T){
  powers <- c(c(1:10), seq(12, 30, 2))
  sft <- pickSoftThreshold(gene_wgc, powerVector = powers, networkType = "signed", verbose = 5)
  sft
  pdf(paste0(outdir, "soft.pdf"), width = 8, height = 6)
  par(mfrow = c(1,2))
  cex1 <- 0.8
  
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
       , xlab = "Soft Threshold (power)"
       , ylab = "Scale Free Topology Model Fit,signed R^2"
       , type = "n"
       , main = "Scale independence"
       , ylim = c(-1,1)
  )
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
       ,labels = powers, cex = cex1, col = "red"
  )
  abline(h = 0.8, col = "red")
  abline(h = 0.85, col = "orange")
  abline(h = 0.9, col = "blue")
  data.frame(power = powers, R2 = -sign(sft$fitIndices[,3])*sft$fitIndices[,2], connect = sft$fitIndices[,5])
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
  
  stop("finish pickpower")
}


if(tom == F){
  power <- 12
  TOM = TOMsimilarityFromExpr(gene_wgc,
                              networkType = 'signed',
                              TOMType = 'signed',
                              power = power)
  colnames(TOM) =rownames(TOM) = colnames(gene_wgc)
  save.image(paste0(outdir, "./wgcna.Rdata"))
}else{
  load(paste0(outdir, "./wgcna.Rdata"))
}


dissTOM=1-TOM
geneTree = hclust(as.dist(dissTOM),method='average')
minModule

dynamicMods = cutreeDynamic(dendro = geneTree,
                            distM = dissTOM,
                            method='hybrid',
                            deepSplit = 2,
                            pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize) #
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods) 
table(dynamicColors)


#module eigengenes
MEList = moduleEigengenes (gene_wgc, colors = dynamicColors) 
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
pdf(paste0(outdir, "moduleCor.pdf"))
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap",
                      marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE,
                      xLabelsAngle = 90)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h=cut, col="red")
dev.off()


# merge module
merge_modules = mergeCloseModules(gene_wgc, dynamicColors, cutHeight = cut, verbose = 3)
mergedColors = merge_modules$colors
mergedMEs = merge_modules$newMEs
pdf(paste0(outdir, "module_cluster.pdf"))
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()



# extract geneset
gene.names <- gene_wgc %>% colnames()
cl=c()
module_gene=c()
module_colors=unique(mergedColors)
for (color in module_colors){
  module=gene.names[which(mergedColors==color)]
  group1=rep(color,times=length(module))
  cl=c(cl,group1)
  module_gene=c(module_gene,module)
}
my_module=data.frame("colour"=cl,"gene"=module_gene)
write.csv(my_module,paste0(outdir, "Modules_gene.csv"))



# MEs
MEs0 = moduleEigengenes(gene_wgc, mergedColors)$eigengenes 
MEs = orderMEs(MEs0)
write.table(MEs, paste0(outdir, "MEs.txt"), quote = F, sep = "\t")

nSamples = nrow(gene_wgc)
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(gene_wgc, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

MM.df=cbind(geneModuleMembership, MMPvalue)


fwrite(MM.df, paste0(outdir, "MM.txt"), sep = "\t", row.names = T)


# KIM
adjacency = adjacency(gene_wgc, power = power);
KIM=intramodularConnectivity(adjacency, mergedColors)
fwrite(KIM, paste0(outdir, "KIM.txt"), sep = "\t", row.names = T)

#算KME
datKME=signedKME(gene_wgc, MEs, outputColumnName="MM.")
fwrite(datKME, paste0(outdir, "KME.txt"), sep = "\t", row.names = T)

## merge module between blue and cyan
mergedColors1 = gsub( "blue", "cyan", mergedColors) %>% gsub("white", "blue", .)
MEs01 = moduleEigengenes(gene_wgc, mergedColors1)$eigengenes #计算合并后的模块的第一主成分ME。
MEs1 = orderMEs(MEs01);
write.table(MEs1, paste0(outdir, "MEsMerge.txt"), quote = F, sep = "\t")

datKME1=signedKME(gene_wgc, MEs1, outputColumnName="MM.")
fwrite(datKME1, paste0(outdir, "KMEMerge.txt"), sep = "\t", row.names = T)





# calculate module change by wilcox test ----
outdir1 <- paste0(outdir, "geneChange/1.sex/mergeCyanBlue/")
dir.create(outdir1)

# function
th.change2 <- function(outdir, module, MEs, group.anno, name = "change", method = "wilcox"){
  data1 <- MEs %>% 
    dplyr::rename_all(.funs = funs(str_sub(., 3))) %>% 
    rownames_to_column("sample") %>% 
    left_join(group.anno)
  
  library(ggpubr)
  
  lapply(module, function(x){
    outdir1 <- paste0(outdir, x, "/")
    dir.create(outdir1)
    
    ggplot(data1, aes(x = status, y = get(x)))+
      geom_boxplot(aes(color = status), outlier.size = 0.5, fill = "white", width = 0.6)+
      scale_y_continuous(expand = c(0.1,0, 0.2, 0))+
      scale_color_manual(values = color_mutli[c(2, 3, 6)])+
      stat_compare_means(comparisons = list(c("Control", "Susceptible"), c("Control", "Resceptible")), size = 6*5/14, method = method, method.args = list(exact = F))+
      theme_minimal()+
      labs(x = "", y = "The relation expression", title = x)+
      facet_wrap(~Region)+
      theme_bw()+
      theme(text = element_text(size = 7), legend.position = "none")
    
    ggsave(paste0(outdir1, name, ".pdf"), width = 4*3, height = 5*3, units = "cm")
  })
  
}





# analysis--------
module_gene1.df <- module_gene.df %>% 
  mutate(colour = gsub("blue", "cyan", colour))

# MEs
exp.df <- exp.d %>%
  column_to_rownames("symbol")

#Male
sex.c <- "Male"
outdir2 <- paste0(outdir1, sex.c, "/")
dir.create(outdir2)
module.c <- module_gene1.df$colour %>% unique()

gene_wgc = t(exp.df)
colour.c <- module_gene1.df %>%
  filter(symbol %in% get(paste0("gene_", sex.c, ".union"))) %>%
  select(colour) %>% unlist
names(colour.c) <- module_gene1.df %>%
  filter(symbol %in% get(paste0("gene_", sex.c, ".union"))) %>%
  select(symbol) %>% unlist
colour.c <- colour.c[colnames(gene_wgc)]
table(names(colour.c) == colnames(gene_wgc))
MEs.sub <- moduleEigengenes(gene_wgc, colour.c)$eigengenes
fwrite(MEs.sub, paste0(outdir2, "MEs_sub.txt"), sep = "\t", row.names = T)

MEs <- fread(paste0(outdir1, "MEs_sub.txt")) %>% 
  column_to_rownames("V1")
sample.c <- group.anno$sample[group.anno$sex == sex.c]

th.change2(outdir = outdir2, module = module.c, MEs = MEs %>% filter(rownames(MEs) %in% sample.c), group.anno = group.anno %>% filter(sample %in% sample.c), method = "wilcox.test", name = "wilcox")

#Female
sex.c <- "Female"
outdir2 <- paste0(outdir1, sex.c, "/")
dir.create(outdir2)
module.c <- module_gene1.df$colour %>% unique()

gene_wgc = t(exp.df)
colour.c <- module_gene1.df %>%
  filter(symbol %in% get(paste0("gene_", sex.c, ".union"))) %>%
  select(colour) %>% unlist
names(colour.c) <- module_gene1.df %>%
  filter(symbol %in% get(paste0("gene_", sex.c, ".union"))) %>%
  select(symbol) %>% unlist
colour.c <- colour.c[colnames(gene_wgc)]
table(names(colour.c) == colnames(gene_wgc))
MEs.sub <- moduleEigengenes(gene_wgc, colour.c)$eigengenes
fwrite(MEs.sub, paste0(outdir2, "MEs_sub.txt"), sep = "\t", row.names = T)

MEs <- fread(paste0(outdir1, "MEs_sub.txt")) %>% 
  column_to_rownames("V1")
sample.c <- group.anno$sample[group.anno$sex == sex.c]

th.change2(outdir = outdir2, module = module.c, MEs = MEs %>% filter(rownames(MEs) %in% sample.c), group.anno = group.anno %>% filter(sample %in% sample.c), method = "wilcox.test", name = "wilcox")





# function------------------------
plot.wgcnaMuduleSum <- function(file, data, outdir, outfile = "moduleSumBar.pdf", color = "orange", title = NULL){
  #判断文件或者数据输入
  if(!is.null(file)){
    data <- fread(file, header = T)
  }else if(!is.null(data)){
    data <- data
  }
  
  #处理数据
  input.df <- data %>% 
    dplyr::group_by(colour) %>% 
    dplyr::summarise(count = n()) %>% 
    dplyr::arrange(-count) %>% 
    dplyr::mutate(colour = factor(colour, levels = colour %>% rev))
  fwrite(input.df, paste0(outdir, "moduleSum.txt"), sep = "\t")
  
  #画图
  plot <- ggplot(input.df, aes(y = colour, x = count))+
    geom_bar(stat = "identity", color= color, width = 0.6, fill = color)+
    geom_tile(aes(x = -100, y = colour, fill = colour), width = 100, height = 0.9)+
    scale_fill_manual(values = levels(input.df$colour))+
    scale_x_continuous(expand = c(0, 0, 0.1, 0))+
    geom_text(aes(label = count), nudge_x = 350, size = 5*5/14)+
    labs(y = "", x = "The numbers of gene in modules", title = title)+
    theme_minimal()+
    theme(legend.position = "none", text = element_text(size = 7),panel.grid = element_line(size = 0.1), plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
  
  ggsave(paste0(outdir, outfile), width = 6, height = 6, units = "cm")
  
  return(plot)
}


# plot---------------------
file.c <- paste0(outdir, "/Modules_gene.csv")
outdir1 <- paste0(outdir, "moduleSum/")
dir.create(outdir1)

plot.wgcnaMuduleSum(file = file.c, outdir = outdir1) 