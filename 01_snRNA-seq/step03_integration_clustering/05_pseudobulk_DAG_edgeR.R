#-----------------------------------------
## modified from the script for pseudo-bulk DEG analysis from Alex.
#  the region could be gene or ACRs.
#-----------------------------------------

#module load R/4.1.0-foss-2019b

# load libraries
library(SingleCellExperiment)
library(Matrix.utils)
library(edgeR)
library(parallel)
library(dplyr)

args <- commandArgs(T)
if(length(args) != 5){stop("Rscript pseudobulk_DAR_analysis.R <meta> <gene_sparse.rds> <markers> <column_id> <prefix>")}

metafile <- as.character(args[1])
gene.sparse <- as.character(args[2])
markerfile <- as.character(args[3])
column_id <- as.character(args[4])
prefix <- as.character(args[5])

# gene.sparse <- "/scratch/xz24199/02Will82/scatac-seq/bed/Gm_atlas_Cotyledon_stage_seeds_metav3_gene_sparse.rds"
# metafile <- "/scratch/xz24199/02Will82/scatac-seq/step03_clustering/v1-2/Gm_Cotyledon_stage_seeds_metav3.corrected.metadata.txt"
# markerfile <- "/scratch/xz24199/02Will82/scatac-seq/markers/test_final_markers_C.txt"

#prefix <- sub("_metav3_gene_sparse.rds","",basename(metafile))

# functions
#edgeRscDEG <- function(sce, ids=c("LouvainClusters", "library"), threads=1, use.ref.cells=F){
edgeRscDEG <- function(sce, ids, threads=1, use.ref.cells=F){
    
    # set up vars
    #ids=c("LouvainClusters", "library")
    groups <- colData(sce)[, ids]
    groups[,1] <- as.factor(groups[,1])
    groups[,2] <- as.factor(groups[,2])
    clusters <- levels(groups[,1])
    threads <- ifelse(length(clusters) > threads, threads, length(clusters))
    
    # iterate
    outs <- mclapply(clusters, function(z){
        
        # set-up reference cells
        if(use.ref.cells){
            dff <- as.data.frame(groups)
            rownames(dff) <- rownames(colData(sce))
            not.cluster <- dff[dff[,1] != z,]
            cells.per.cluster <- nrow(dff[dff[,1] == z,])
            ran.cells <- sample(rownames(not.cluster), cells.per.cluster)
            #message(" - selected ", length(ran.cells), " control cells ...")
        }
        
        # verbose
        message(" - identifying DEG from cluster ", z)
        
        # rename groups
        df <- as.data.frame(groups)
        rownames(df) <- rownames(colData(sce))
        df[,1] <- as.character(df[,1])
        df$cluster_id <- as.factor(ifelse(df[,1] == z, 1, 0))
        df[,1] <- NULL
        if(use.ref.cells){
            df <- df[ifelse(df$cluster_id==1 | rownames(df) %in% ran.cells, T, F),]
            message(" - number of cluster cells = ", nrow(subset(df, df$cluster_id==1)))
            message(" - number of control cells = ", nrow(subset(df, df$cluster_id==0)))
        }
        #message(" - total nummber of cells in test = ", nrow(df))
        n.genes <- Matrix::colSums(counts(sce)[,rownames(df)] > 0)
        sample.data <- aggregate(n.genes~df$cluster_id+df[,ids[2]], FUN=mean)
        sample.data <- sample.data[order(sample.data[,1], decreasing=F),]
        
        # aggregate by cluster/library
        pb <- aggregate.Matrix(t(counts(sce)[,rownames(df)]), 
                               groupings = df, fun = "sum")
        pb <- t(pb)
        colnames(sample.data) <- c("cluster_id",ids[2], "n.genes")
        #print(head(pb))
        #print(head(sample.data))
        
        # create edgeR object
        group <- factor(sample.data$cluster_id, levels=c(0,1))
        #batch <- factor(sample.data[,ids[2]])  #no use, might could put in the model.matrix function by add model.matrix(~ group+batch)
        #ave.genes <- scale(sample.data$n.genes) # no use, might could put in the DGEList by add norm.factors = ave.genes
        dge <- DGEList(counts = pb, 
                       norm.factors = rep(1, length(pb[1,])), 
                       group = group)

        # design experiment and estimate norm factors/dispersion
        design <- model.matrix(~ group)
        dge <- calcNormFactors(dge, method = "TMM", logratioTrim=0.1)
        dge <- estimateDisp(dge, design = design)
        
        # estimate Differential expression
        fit <- glmQLFit(dge, design)
        res <- glmQLFTest(fit, coef=ncol(fit$coefficients))
        res$table$FDR <- p.adjust(res$table[,4], method="fdr")
        res$table$cluster_id <- z
        res$table$geneID <- rownames(res$table)
        rownames(res$table) <- seq(1:nrow(res$table))
        message("   ~ returning results for cluster ",z)
        return(res$table)
        
    }, mc.cores=threads)
    outs <- do.call(rbind, outs)
	outs <- as.data.frame(outs)
    return(outs)
    
}

# load data
message(" - loading data ...")
counts <- readRDS(gene.sparse)
metadata <- read.table(metafile)
shared <- intersect(rownames(metadata), colnames(counts))
metadata <- metadata[shared,]
counts <- counts[,shared]
counts <- counts[Matrix::rowSums(counts) > 0,]
sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = metadata)
sce <- sce[rowSums(counts(sce) > 0) >= 0, ]
colData(sce)$batch <- colData(sce)$library

#----
markers <- read.delim(markerfile, header = T)
#row.names(markers) <- markers$geneID
#markers$label <- paste(markers$cell_type,markers$Symbol, sep = ':') #format marker label
#markers <- read.table(as.character(args[3]), header=T)
#markers <- markers[!duplicated(markers$geneID),]
#type <- markers$type
#name <- markers$name
#names(type) <- markers$geneID
#names(name) <- markers$geneID
#functions <- read.delim("/scratch/apm25309/reference_genomes/Zmays/v5/B73v5AGP_tair10_functions.txt", quote="", header=F)
#-------

# filter low frequency, in-accessible genes
sce <- sce[rowMeans(counts(sce) > 0) > 0, ]

# iterate over all clusters
message(" - running DAR analysis by cluster with replicates ...")
deg <- edgeRscDEG(sce, threads=16, ids=c(column_id, "library"), use.ref.cells=T)

# add gene functional annotation
#cur <- as.character(functions$V2)
#comp <- as.character(functions$V3)
#names(cur) <- as.character(functions$V1)
#names(comp) <- as.character(functions$V1)

# filter for significant genes (positive)
deg <- deg[order(deg$cluster_id, deg$PValue, decreasing=F),]
#deg$desc <- cur[as.character(deg$geneID)]
#deg$comp_desc <- comp[as.character(deg$geneID)]
write.table(deg, file=paste0(prefix,"_DAG_pseudobulk.txt"), quote=F, row.names=T, sep="\t", col.names=T)
sig <- subset(deg, deg$logFC > 0)
write.table(sig, file=paste0(prefix,"_DAG_pseudobulk.log2fc0.txt"), quote=F, row.names=T, col.names=T, sep="\t")

# filter markers
sig.markers <- deg[as.character(deg$geneID) %in% as.character(markers$geneID),]
sig.markers <- left_join(sig.markers, markers, by = "geneID")
#sig.markers$type <- type[sig.markers$geneID]
#sig.markers$name <- name[sig.markers$geneID]
write.table(sig.markers, file=paste0(prefix,"_DAG_pseudobulk.marker_genes.txt"), quote=F, row.names=T, col.names=T, sep="\t")
