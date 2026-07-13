# run GSEA analysis #

# arguments
args <- commandArgs(T)
if(length(args) != 4){stop("Rscript gsea_analysis.R <DESeq2.output> <GO_ids> <prefix> <threads>")}
input <- as.character(args[1])
goids <- as.character(args[2])
prefix <- as.character(args[3])
threads <- as.numeric(args[4])

# load libraries
library(dplyr)
library(gage)
library(fgsea)
library(reshape2)
library(parallel)

# functions
convertTerms <- function(x){
    out <- lapply(seq(1:nrow(x)), function(z){
        gids <- strsplit(as.character(x$V2[z]), "\\,")
        gids <- do.call(c, gids)
        return(gids)
    })
    names(out) <- as.character(x$V1)
    return(out)
}
formatInput <- function(x, type="F"){
    x$cluster_id <- as.character(x$cluster_id)
    outs <- lapply(unique(x$cluster_id), function(z){
        clust <- subset(x, x$cluster_id==z)
        if(type=="F"){
            metric <- clust$F
        }else if(type=="logFC"){
            metric <- clust$logFC
        }else{
		metric <- clust[,type]
	}
        names(metric) <- clust$geneID
        metric <- metric[order(metric, decreasing=T)]
        return(metric)
    })
    names(outs) <- unique(x$cluster_id)
    return(outs)
}
GSEA <- function(gene_list, myGO, pval) {
    set.seed(54321)
    
    if ( any( duplicated(names(gene_list)) )  ) {
        warning("Duplicates in gene names")
        gene_list = gene_list[!duplicated(names(gene_list))]
    }
    if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
        warning("Gene list not sorted")
        gene_list = sort(gene_list, decreasing = TRUE)
    }
    print(head(gene_list))
    print(head(myGO))
    message(" - running FGSEA ... ")
    fgRes <- as.data.frame(fgsea::fgsea(pathways = myGO, 
                          stats = gene_list,
                          minSize=10,
                          maxSize=600,
                          nperm=10000))
    
    ## Filter FGSEA by using gage results. Must be significant and in same direction to keep 
    fgRes <- fgRes[order(fgRes$NES, decreasing=T),]
    fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
    print(head(fgRes))   
    return(fgRes)
}
run_GSEA_parallel <- function(clusters, GO, threshold=0.05, threads=1){
    
    # check
    print(str(clusters))
    print(str(GO))

    # run GSEA in parallel
    outs <- lapply(names(clusters), function(x){
        message(" - running GSEA for ", x, " ...")
        input <- clusters[[x]]
        results <- GSEA(input, GO, threshold)
        results$cluster <- x
	print(head(results))
        return(results)
    })#, mc.cores=threads)
    names(outs) <- names(clusters)
    
    # merge
    all <- do.call(rbind, outs)
    return(all)
}

# read data
message(" - reading edgeR input ...")
a <- read.delim(input, quote="")
message(" - reformating edgeR data frame ...")
a <- formatInput(a, type="logFC")
message(" - read geneID GO terms ...")
b <- read.table(goids)
message(" - reformating geneID GO terms ...")
b <- convertTerms(b)

# run in parallel
res <- run_GSEA_parallel(a, b, threads=threads)
res$leadingEdge <- NULL

# write results to disk
write.table(res, file=paste0(prefix,"_ct_gsea_annotation.txt"), quote=F, row.names=T, col.names=T, sep="\t")
