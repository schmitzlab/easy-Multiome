## plot marker heatmap ##

#module load R/4.1.0-foss-2019b

# load libraries
library(pheatmap)
library(RColorBrewer)
library(edgeR)
library(Matrix)



# load arguments
args <- commandArgs(T)
if(length(args)!=5){stop("Rscript PlotClusterZscore.R <metadata> <gene.sparse.rds> <markers> <cluster_id> <name>")}

#a <- read.table("Gm_Cotyledon_stage_seeds_50000_feature.corrected.metadata.txt")
#markers <- read.table("/scratch/xz24199/02Will82/scatac-seq/markers/test_final_markers_C.txt", header=T)
#rds <- readRDS("Gm_Cotyledon_stage_seeds_50000_feature_imputed_data.rds")
#prefix <- "Gm_Cotyledon_stage_seeds_50000_feature"

# load data
a <- read.table(as.character(args[1]), comment.char="")
markers <- read.table(as.character(args[3]), header=T)
rds <- readRDS(as.character(args[2]))
cluster_id <- as.character(args[4])
prefix <- as.character(args[5])

b <- rds$impute.activity

#a <- subset(a, type == type)

#b <- rds$impute.activity

# aggregate clusters
clusts <- sort(unique(a[,cluster_id]))
out <- lapply(clusts, function(x){
    df <- a[a[,cluster_id] == x,]
    Matrix::rowSums(b[,colnames(b) %in% rownames(df)])
})
out <- do.call(cbind, out)
#colnames(out) <- paste0(cluster_id,"_", clusts)
colnames(out) <- clusts

# normalize
out <- cpm(out, log=F)
gene.aves <- rowMeans(out)
out <- out[gene.aves > 0,]
vars <- apply(out, 1, var)
out <- out[vars > 0,]

# estimate log2FC
ids <- colnames(out)
fc <- lapply(ids, function(x){
    log2((out[,x])/rowMeans(out[,!colnames(out) %in% x]))
})
fc <- do.call(cbind, fc)
colnames(fc) <- ids
fc[is.infinite(fc) & fc < 0] <- -1*max(fc[is.finite(fc)])
fc[is.infinite(fc) & fc > 0] <- max(fc[is.finite(fc)])
fc[is.na(fc)] <- 0

# get zscores
zscore <- as.matrix(t(scale(t(out))))

# save tables
#write.table(out, file=paste0(prefix,"_clusters.CPM.txt"), quote=F, row.names=T, col.names=T, sep="\t")
#write.table(fc, file=paste0(prefix,"_clusters.log2FC.txt"), quote=F, row.names=T, col.names=T, sep="\t")
#write.table(zscore, file=paste0(prefix,"_clusters.zscore.txt"), quote=F, row.names=T, col.names=T, sep="\t")

# cluster columns
fc.clust <- hclust(dist(t(fc)))$order
z.clust <- hclust(dist(t(zscore)))$order
fc <- fc[,fc.clust]
zscore <- zscore[,z.clust]

# subset to markers
zscore.m <- zscore[rownames(zscore) %in% as.character(markers$geneID),]
fc.m <- fc[rownames(fc) %in% as.character(markers$geneID),]

# cluster rows
f.row <- apply(fc.m, 1, which.max)
z.row <- apply(zscore.m, 1, which.max)
fc.m <- fc.m[order(f.row, decreasing=F),]
zscore.m <- zscore.m[order(z.row, decreasing=F),]

# geneIDs
gname <- paste(markers$name, markers$type, markers$geneID, sep = ":")
names(gname) <- markers$geneID
rownames(fc.m) <- gname[rownames(fc.m)]
rownames(zscore.m) <- gname[rownames(zscore.m)]

# plot
write.table(fc.m, file=paste0(prefix,"_Log2FC_marker_genes_imputed.heatmap.txt"), quote=F, sep="\t")
pdf(paste0(prefix,"_Log2FC_marker_genes_imputed.heatmap.pdf"), width=8, height=24)
pheatmap(fc.m, col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         cluster_col=F, cluster_row=F)
dev.off()

write.table(zscore.m, file=paste0(prefix,"_Zscore_marker_genes_imputed.heatmap.txt"), quote=F, sep="\t")
pdf(paste0(prefix,"_Zscore_marker_genes_imputed.heatmap.pdf"), width=8, height=24)
val <- max(abs(zscore.m))
pheatmap(zscore.m, col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         cluster_col=F, cluster_row=F, breaks=seq(from= -val, to=val, length.out=101),
         fontsize_row=8)
dev.off()




