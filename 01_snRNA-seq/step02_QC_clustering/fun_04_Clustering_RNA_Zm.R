##########################################
####        Seurat Clustering         ####
##########################################

##=============================== Load package ===============================
#module load R/4.3.1-foss-2022a

#----to install different version of seurate:
#remotes::install_version("SeuratObject", "4.1.4")
#remotes::install_version("Seurat", "4.4.0", upgrade = FALSE)

library(Seurat)
library(dplyr)
library(ggplot2)
library(monocle3)
#library(scCustomize)
#library(qs)   #quickly writing and reading any R object to and from disk.
#library(psych)  ## for stats describeby()
library(DoubletFinder)
library(harmony)
library(purrr) #for reduce funciton
library(vioplot)
#library(scCustomize)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(patchwork)
library(future)

options(future.globals.maxSize = 8 * 1024^3)  # 8 GB


args <- commandArgs(T)
if(length(args) != 1){stop("Rscript snRNA-seq_clustering_xz.R <prefix>")}

#path <- as.character(args[1])
#tissue <- as.character(args[2])
prefix <- as.character(args[1])

#path <- "/scratch/xz24199/07USB_infected_root/snRNA-seq/step02_QC_cells"
#tissue <- "Gm_Infected_root"
#prefix <- "Gm_Infected_root_integrated"

#get obj data
obj <- paste0(prefix,".raw.seurat.rds")
#objs <- list.files(path=path, pattern = "*.raw.seurat.rds", full.names=T)
#objs <- objs[grepl(tissue, objs)]

#get meta data
meta <- paste0(prefix,".updated_metadata_v2.txt")
#metas <- list.files(path=path, pattern = "*.updated_metadata_v3.txt", full.names=T)
#metas <- metas[grepl(tissue, metas)]

#---update obj based on meta and merge
#objs <- lapply(objs, function(x){
  #x <- "/scratch/xz24199/07USB_infected_root/snRNA-seq/step02_QC_cells/Gm_Infected_root-rep1.raw.seurat.rds"
  #name <- sub(".raw.seurat.rds","",basename(x))
  #meta <- metas[grepl(name, metas)]
  obj <- readRDS(obj)
  meta <- read.table(meta)
  #update counts and meta
  counts <- GetAssayData(obj, slot="counts", assay="RNA")
  counts <- counts[!grepl("^Ze",rownames(counts)),] #filter nematode gene
  overlap <- intersect(rownames(meta), colnames(counts))
  counts <- counts[,overlap]
  meta <- meta[overlap,]
  # create new Seurat object
  obj <- CreateSeuratObject(counts=counts, project = prefix, meta.data=meta)
#  return(obj.n)
#})

#saveRDS(objs, file=paste0(tissue,".seuratObj.raw.rds"))

message("------normalize and clustering---")

#objs <- readRDS("Gm_Infected_root.seuratObjs.list.RDS")
# normalize and filter doublets.
#objs <- lapply(objs, function(x){
  obj <- obj %>% 
    SCTransform(vars.to.regress = 'pOrg') %>% 
    RunPCA(verbose = F, npcs = 25) %>% 
    RunUMAP(reduction = "pca", dims = 1:25)
  
  nCell = length(obj$cellID)  
  nExpRate = round(nCell/1000) * 0.008
  nExp <- round(nCell * nExpRate)         # expected number of doublets  (0.8% double rate, 8 doublets in 1k cells called)
  message(Project(obj),' Cell number: ', nCell, ". Expected doublets: ", nExp, "(", nExpRate, ")")
  obj <- doubletFinder(obj, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:25, sct=T)   # find doublet
  colnames(obj@meta.data)[ncol(obj@meta.data)] <- "DoubletFinder"   # change the colname of the last column (DoubletFinder result)
#  return(obj)
#})

#merge object.
#obj.m <- purrr::reduce(objs, function(x, y) {
#  merge(x = x, y = y, project = prefix)})
saveRDS(obj, file=paste0(prefix,".seuratObj.raw.rds"))

#---Plot QC doublet calling---
meta.s <- obj@meta.data
write.table(meta.s, file=paste0(prefix,".rna.DoubletFinder.metadata.txt"), sep = "\t", quote=F)

pdf(file=paste0(prefix, '_Doublet_QC.pdf'), width = 12, height = 3)
p1 <- ggplot(meta.s, aes(x = DoubletFinder, y = log10nUMI, fill = library)) +
  geom_violin(position = position_dodge(0.7), trim = FALSE) + 
  geom_boxplot(position = position_dodge(0.7), color = "white", 
               width = 0.05, show.legend = FALSE) + theme_classic()+xlab("") + ylab("Total nUMI (log10)") +
  theme(legend.position='top')
p2 <- ggplot(meta.s, aes(x = DoubletFinder, y = log10nGene, fill = library)) + 
  geom_violin(position = position_dodge(0.7), trim = FALSE) + 
  geom_boxplot(position = position_dodge(0.7), color = "white", 
               width = 0.05, show.legend = FALSE) + theme_classic()+xlab("") + ylab("Total nGene (log10)") +
  theme(legend.position='top')
p3 <- ggplot(meta.s, aes(x = DoubletFinder, y = pOrg, fill = library)) + 
  geom_violin(position = position_dodge(0.7), trim = FALSE) + 
  geom_boxplot(position = position_dodge(0.7), color = "white", 
               width = 0.05, show.legend = FALSE) + theme_classic()+xlab("") + ylab("Proportion of organelle") +
  theme(legend.position='top')
p4 <- ggplot(meta.s, aes(x = DoubletFinder, y = spliced_rate, fill = library)) + 
  geom_violin(position = position_dodge(0.7), trim = FALSE) + 
  geom_boxplot(position = position_dodge(0.7), color = "white", 
               width = 0.05, show.legend = FALSE) + theme_classic()+xlab("") + ylab("Proportion of spliced reads") +
  theme(legend.position='top')
wrap_plots(p1, p2, p3, p4, nrow=1)
dev.off()



#---Keep singlet and Plot QC---
message("---Keep singlet and Plot QC---")
obj.m <- subset(x = obj, subset = DoubletFinder == "Singlet")

meta.m <- obj.m@meta.data
write.table(meta.m, file=paste0(prefix,".rna.Singlet.metadata.txt"), sep = "\t", quote=F)

# number <- length(unique(meta.m$library))
# #cols <- colorRampPalette(brewer.pal(12, "Paired"))(number)
# cols <- brewer.pal(12, "Paired")[1:number]
# pdf(file=paste0(prefix,"_QC_violin_Rmdoublet.pdf"), height = 3, width = 6)
# layout(matrix(c(1:3), nrow=1))
# vioplot(meta.m$log10nUMI ~ meta.m$library,col = cols, las=3, ylab="nUMI (log10)", xlab="")
# vioplot(meta.m$log10nGene ~ meta.m$library,col = cols, las=3, ylab="nGene (log10)", xlab="")
# vioplot(meta.m$pOrg ~ meta.m$library,col = cols, las=3, ylab="Proportion of organelle", xlab="")
# #vioplot(meta.m$pNmt ~ meta.m$library,col = cols, las=3, ylab="Proportion of nematode", xlab="")
# dev.off()

# #---Integration with harmony
# message("---Integration with harmony---")
# #select features that are repeatedly variable across datasets for integration
# objs <- lapply(objs, function(x){
  # x <- subset(x = x, subset = DoubletFinder == "Singlet")
  # x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  # return(x)
# })
# features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 3000)

# #run integration
# VariableFeatures(obj.m) <- features
# obj.m <- obj.m %>% 
  # RunPCA(npcs = 25, verbose = FALSE) %>%
  # RunHarmony(group.by.vars = "library", assay.use = "SCT")

#---find clusters
#obj.m <- readRDS("Gm_Infected_root.MergedSeuratObjs.tmp.RDS")
obj.m <- FindVariableFeatures(obj.m, selection.method = "vst", nfeatures = 3000)
message("---Finding clusters---")
obj.m <- obj.m %>%
  RunUMAP(assay = "SCT",reduction = "pca", dims = 1:25) %>%
  FindNeighbors(assay = "SCT",reduction = "pca", dims = 1:25) %>%
  FindClusters(resolution = 0.5)
#obj.m$tech <- "dmulti_rna"
Idents(obj.m) <- obj.m$seurat_clusters
obj.m$celltype <- paste0("RNA_",Idents(obj.m))
saveRDS(obj.m, file=paste0(prefix,".SeuratObj.cluster.rds"))

#---Exporting rds, meta, harmony---
message("---Exporting rds, meta, harmony---")
#export gene x cell matrix
data <- as.matrix(GetAssayData(object = obj.m, slot = "counts")) #export count matrix data(gene x cells)
saveRDS(data, file = paste0(prefix,".rna.sparse.rds"))

#export meta data with umap
meta <- obj.m@meta.data #export count matrix data(cells x attributes)
umap <- as.data.frame(Embeddings(object = obj.m[["umap"]]))
colnames(umap) <- c("umap1","umap2") 
#merge meta
id <- intersect(row.names(meta), row.names(umap))
meta <- meta[id,]
umap <- umap[id,]
meta <- cbind(meta,umap)
write.table(meta, file=paste0(prefix,".rna.Singlet.Cluster.metadata.txt"), sep = "\t", quote=F)

#---export pca---
pca <- as.data.frame(Embeddings(object = obj.m, reduction = "pca"))
write.table(pca, file=paste0(prefix,".rna.PCA.txt"), sep = "\t", quote=F)

# #---export harmony---
# pca <- as.data.frame(Embeddings(object = obj.m, reduction = "harmony"))
# write.table(pca, file=paste0(prefix,".rna.Harmony.txt"), sep = "\t", quote=F)

message("---Plot clusters---")
p <- DimPlot(obj.m, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend() + ggtitle(prefix)
ggsave(paste0(prefix,"_RNA_UMAP.pdf"), plot = p, width = 4, height = 4)

#PlotUMAP2 function modified from Socrate(https://github.com/plantformatics/Socrates/blob/main/R/visualizations.R)
plotUMAP2 <- function(meta,
                     column="seurat_clusters",
                     cex=0.3,
                     opaque=1,
                     xlab="umap1",
                     ylab="umap2",
                     main=""){
  
  # # set b as meta data
  # if(is.null(obj[[cluster_slotName]])){
  #   #stop(" - ERROR: final.meta slot from callClusters is missing from object ...")
  #   cluster_slotName <- "meta"
  # }
  b <- meta
  b <- b[complete.cases(b$umap1),]
  b <- b[complete.cases(b$umap2),]
  
  # test if column is present
  if(!column %in% colnames(b)){
    stop(" - ERROR: column header, ", column, ", is missing from ",cluster_slotName, " ...")
  }
  
  # cols
  if(is.factor(b[,c(column)])){
    b <- b[sample(nrow(b)),]
    cols <- colorRampPalette(brewer.pal(12,"Paired")[1:10])(length(unique(b[,column])))
    colv <- cols[factor(b[,column])]
  }else if(is.character(b[,column])){
    b[,column] <- factor(b[,column])
    b <- b[sample(nrow(b)),]
    cols <- colorRampPalette(brewer.pal(12,"Paired")[1:10])(length(unique(b[,column])))
    colv <- cols[factor(b[,column])]
  }else if(is.numeric(b[,column])){
    b <- b[order(b[,column], decreasing=F),]
    cols <- viridis(100)
    colv <- cols[cut(b[,column], breaks=101)]
  }
  
  # plot
  plot(b$umap1, b$umap2, pch=16, cex=cex, col=alpha(colv,opaque),
       xlab=xlab,
       ylab=ylab,
       main=main,
       xlim=c(min(b$umap1), max(b$umap1)+(abs(max(b$umap1))*0.5)))
  
  if(is.factor(b[,column])){
    legend("right", legend=sort(unique(b[,column])),
           fill=cols[sort(unique(b[,column]))])
  }
}

pdf(paste0(prefix,".UMAP.harmony.QC.pdf"), width=12, height=3)
layout(matrix(c(1:4), nrow=1, ncol = 4))
plotUMAP2(meta, column="log10nUMI", cex=0.2, main = "Total nUMI(log10)")
plotUMAP2(meta, column="log10nGene", cex=0.2, main = "Total nGene(log10)")
plotUMAP2(meta, column="pOrg", cex=0.2, main = "Fraction reads in organelle")
#plotUMAP2(meta, column="pNmt", cex=0.2, main = "Fraction reads in nematode")
dev.off()

# message("---QC clusters---")
# #-------Plot distribution----
# library(dplyr)
# meta <- as.data.frame(meta)
# counts <- meta %>% group_by(celltype, library) %>% summarise(n = n())
# #counts <- meta %>% count("celltype","library")
# #?count
# #counts
# counts.w <- data.frame(c(unique(counts$library)))
# colnames(counts.w) <- "library"
# #head(counts.w)
# #convert to count matrix
# for (i in unique(counts$celltype)) {
  # df <- counts[counts$celltype == i, c("library","n")]
  # colnames(df) <- c("library",i)
  # counts.w <- full_join(counts.w, df, by = "library")
# }

# counts.w[is.na(counts.w)] <- 0
# rownames(counts.w) <- counts.w$library
# counts.w <- counts.w[,-1]
# #head(counts.w)

# counts.w <- data.matrix(counts.w) #convert to numberic

# counts.r <- apply(counts.w, 1, function(x){x*100/sum(x,na.rm=T)})
# write.table(counts.r, file = paste0(prefix,"_nuclei_proprotion_clusters.txt"),sep = "\t", quote = F)

# counts.t <- t(counts.r)
# #mycol <- c("#C6DBEF", "#22B2AE", "#F4B9A2", "#CB181D")
# mycol <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta$library)))
# #mycol

# pdf(paste0(prefix,"_nuclei_proprotion_clusters.pdf"), width=8, height=4)
# barplot(counts.t, col=mycol, xlab = "cluster_name", border = "white", ylab = "Proportion of nuclei(%)", cex.names = 0.8, las = 2, beside=T)
# legend("topright",unique(meta$library), fill=mycol, cex=0.5)
# dev.off()

#Cluster QC distributions
#meta <- read.table("Gm_Infected_root_integrated.rna.Singlet.metadata.txt")
meta$seurat_clusters <- as.character(meta$seurat_clusters)
pdf(file=paste0(prefix, '_Clustering_QC.pdf'), width = 9, height = 12)
p1 <- ggplot(meta, aes(x = seurat_clusters, y = log10nUMI, fill = library)) +
  geom_violin(position = position_dodge(0.7), trim = FALSE) + 
  geom_boxplot(position = position_dodge(0.7), color = "white", 
               width = 0.05, show.legend = FALSE) + theme_classic()+xlab("") + ylab("Total nUMI (log10)") +
  theme(legend.position='top')
p2 <- ggplot(meta, aes(x = seurat_clusters, y = log10nGene, fill = library)) + 
  geom_violin(position = position_dodge(0.7), trim = FALSE) + 
  geom_boxplot(position = position_dodge(0.7), color = "white", 
               width = 0.05, show.legend = FALSE) + theme_classic()+xlab("") + ylab("Total nGene (log10)") +
  theme(legend.position='top')
p3 <- ggplot(meta, aes(x = seurat_clusters, y = pOrg, fill = library)) + 
  geom_violin(position = position_dodge(0.7), trim = FALSE) + 
  geom_boxplot(position = position_dodge(0.7), color = "white", 
               width = 0.05, show.legend = FALSE) + theme_classic()+xlab("") + ylab("Proportion of organelle") +
  theme(legend.position='top')
p4 <- ggplot(meta, aes(x = seurat_clusters, y = spliced_rate, fill = library)) + 
  geom_violin(position = position_dodge(0.7), trim = FALSE) + 
  geom_boxplot(position = position_dodge(0.7), color = "white", 
               width = 0.05, show.legend = FALSE) + theme_classic()+xlab("") + ylab("Proportion of spliced reads") +
  theme(legend.position='top')
wrap_plots(p1, p2, p3, p4, ncol=1)
dev.off()



