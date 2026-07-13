###################################################################################################
## plot marker accessibility scores from Alex 20230620 for scifiATAC-seq maize 10x data
###################################################################################################

# load arguments
args <- commandArgs(trailingOnly=T)
if(length(args) != 7){stop("Rscript plot_marker_accessibility.R [meta] [gene_activity] [pcs.txt] [markers.bed] [threads] [prefix] [cluster_id]")}

#args
meta <- as.character(args[1])
geneact <- as.character(args[2])
pcs <- as.character(args[3])
mark <- as.character(args[4])
threads <- as.numeric(args[5])
prefix <- as.character(args[6])
cluster_id <- as.character(args[7])

# meta <- "Gm_atlas_Cotyledon_stage_seeds_var0.5.iNMF_metadata.txt"
# geneact <- "/scratch/xz24199/02Will82/scrna-seq/step04_annotation/Gm_atlas_Cotyledon_stage_seeds.rna.celltype.sparse.rds"
# pcs <- "Gm_atlas_Cotyledon_stage_seeds_var0.5.iNMF_embedding.txt"
# mark <- "/scratch/xz24199/02Will82/scatac-seq/markers/test_final_markers_C.txt"
# type <- "RNA"

# meta <- "/scratch/xz24199/02Will82/scrna-seq/step03_ziliang/Gm_atlas_Root.rna.Singlet.metadata.txt"
# geneact <- "/scratch/xz24199/02Will82/scrna-seq/step03_ziliang/Gm_atlas_Root.rna.sparse.rds"
# pcs <- "/scratch/xz24199/02Will82/scrna-seq/step03_ziliang/Gm_atlas_Root.rna.Harmony.txt"
# mark <- "/scratch/xz24199/02Will82/scatac-seq/markers/test_final_markers_C.txt"
# cluster_id <- "celltype"

# load functions
#source("functions.plot_marker_accessibility.R")
#source("functions.plot_marker_accessibility_xz.R")
#source("functions.plot_marker_accessibility_v2_xz.R")
source("functions.plot_marker_accessibility_v3_xz.R")

# load data
meta <- read.table(meta, comment.char="")
#meta$celltype <- paste0("RNA_",meta$celltype)
#meta <- subset(meta, type == type)
dat <- loadData(meta, pcs, geneact, mark, cluster_id)
b.meta <- dat$b
activity.all <- dat$activity
h.pcs1 <- dat$h.pcs
marker.info.dat <- dat$marker.info

# keep only top 250k cells
# b.meta <- head(b.meta[order(b.meta$total, decreasing=T),], n=100000)
activity.all <- activity.all[,rownames(b.meta)]
activity.all <- activity.all[Matrix::rowSums(activity.all)>0,]
h.pcs1 <- h.pcs1[rownames(b.meta),]
marker.info.dat <- marker.info.dat[rownames(marker.info.dat) %in% rownames(activity.all),]

# iterate over each major cluster
out <- runMajorPriori(b.meta, activity.all, h.pcs1, marker.info.dat, threads=threads, output=prefix, smooth.markers=T)
saveRDS(out,paste0(prefix,'_imputed_data.rds'))
