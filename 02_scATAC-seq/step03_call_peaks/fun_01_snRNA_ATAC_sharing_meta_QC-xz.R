###############################################################################
## isolate shared cell in ATAC and RNA
###############################################################################

 # module load R/4.4.2-gfbf-2024a
 # module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
 # ml R-bundle-CRAN/2024.11-foss-2024a

library(MASS)
library(viridis)
library(patchwork)
library(ggplot2)
library(Seurat)


# args
args <- commandArgs(T)
if(length(args)!=3){stop("Rscript snRNA_scATAC_meta_sharing_QC-xz.R [rna_meta] [atac_path] [name]")}

# load args
rna_meta <- as.character(args[1])
atac_path <- as.character(args[2])
#gene_sparse <- as.character(args[3])
prefix <- as.character(args[3])

# rna_meta <- "/scratch/xz24199/04_scifi_multiome/Zm_seedling/dmulti-RNA/step03-2_rm0_reclustering/Zm_All_seedling_integrated.rna.Singlet.metadata.txt"
# atac_path <- "/scratch/xz24199/04_scifi_multiome/Zm_seedling/dmulti-ATAC-chromap/02_QC"
# prefix <- "AH15AH17_eMulti_Zm_seedling"


#----check all clean cells----
meta.r <- read.table(rna_meta)
#meta.r$call[is.na(meta.r$call)] <- 0

#merge all atac metas
metas <- list.files(path=atac_path, pattern=".updated_metadata_v2.txt", full.names = TRUE)

meta.a <- lapply(metas,function(x){read.table(x)})
meta.a <- do.call(rbind,meta.a)
meta.a$library <- substr(meta.a$cellID, 23, nchar(meta.a$cellID))

#check overlapped cells
meta.r$tenx_CB <- substr(meta.r$cell_barcode, 1, 16)
meta.r$tenx_CB_lib <- paste(meta.r$library,meta.r$tenx_CB,sep=":")

meta.a$tenx_CB <- as.data.frame(do.call(rbind, strsplit(meta.a$cellID,"[:-]", perl=T)))[,3]
meta.a$tenx_CB_lib <- paste(meta.a$library,meta.a$tenx_CB,sep=":")

overlap <- intersect(meta.r$tenx_CB_lib, meta.a$tenx_CB_lib) #overlap# 7999

#update meta with overlap cells
meta.r$sharing <- ifelse(meta.r$tenx_CB_lib %in% overlap, "Shared","Unique")
meta.a$sharing <- ifelse(meta.a$tenx_CB_lib %in% overlap, "Shared","Unique")

meta.r$ATAC_cellID <- paste0("CB:Z:",meta.r$tenx_CB,"-",meta.r$library)

#---Plot Share vs Unique QC plot---
#Plot RNA
cell_counts <- table(meta.r$sharing)
meta.r$sharing <- factor(
  meta.r$sharing,
  levels = names(cell_counts),
  labels = paste0(names(cell_counts), "\n(N=", cell_counts, ")")
)
pdf(file=paste0(prefix, '_Sharing_QC-RNA.pdf'), width = 12, height = 3)
p1 <- ggplot(meta.r, aes(x = sharing, y = log10nUMI, fill = library)) +
  geom_violin(position = position_dodge(0.7), trim = FALSE) + 
  geom_boxplot(position = position_dodge(0.7), color = "white", 
               width = 0.05, show.legend = FALSE) + theme_classic()+xlab("") + ylab("Total nUMI (log10)") +
  theme(legend.position='top')
p2 <- ggplot(meta.r, aes(x = sharing, y = log10nGene, fill = library)) + 
  geom_violin(position = position_dodge(0.7), trim = FALSE) + 
  geom_boxplot(position = position_dodge(0.7), color = "white", 
               width = 0.05, show.legend = FALSE) + theme_classic()+xlab("") + ylab("Total nGene (log10)") +
  theme(legend.position='top')
p3 <- ggplot(meta.r, aes(x = sharing, y = pOrg, fill = library)) + 
  geom_violin(position = position_dodge(0.7), trim = FALSE) + 
  geom_boxplot(position = position_dodge(0.7), color = "white", 
               width = 0.05, show.legend = FALSE) + theme_classic()+xlab("") + ylab("Proportion of organelle") +
  theme(legend.position='top')
p4 <- ggplot(meta.r, aes(x = sharing, y = spliced_rate, fill = library)) + 
  geom_violin(position = position_dodge(0.7), trim = FALSE) + 
  geom_boxplot(position = position_dodge(0.7), color = "white", 
               width = 0.05, show.legend = FALSE) + theme_classic()+xlab("") + ylab("Proportion of spliced reads") +
  theme(legend.position='top')
wrap_plots(p1, p2, p3, p4, nrow=1)
dev.off()

meta.r$sharing <- ifelse(meta.r$tenx_CB_lib %in% overlap, "Shared","Unique") #rm the "\n"


#Plot ATAC
pdf(file=paste0(prefix, '_Sharing_QC-ATAC.pdf'), width = 12, height = 3)
p1 <- ggplot(meta.a, aes(x = sharing, y = log10nSites, fill = library)) +
  geom_violin(position = position_dodge(0.7), trim = FALSE) + 
  geom_boxplot(position = position_dodge(0.7), color = "white", 
               width = 0.05, show.legend = FALSE) + theme_classic()+xlab("") + ylab("Total nSite (log10)") +
  theme(legend.position='top')
p2 <- ggplot(meta.a, aes(x = sharing, y = pTSS, fill = library)) + 
  geom_violin(position = position_dodge(0.7), trim = FALSE) + 
  geom_boxplot(position = position_dodge(0.7), color = "white", 
               width = 0.05, show.legend = FALSE) + theme_classic()+xlab("") + ylab("Proportion in TSS") +
  theme(legend.position='top')
p3 <- ggplot(meta.a, aes(x = sharing, y = FRiP, fill = library)) +
  geom_violin(position = position_dodge(0.7), trim = FALSE) + 
  geom_boxplot(position = position_dodge(0.7), color = "white", 
               width = 0.05, show.legend = FALSE) + theme_classic()+xlab("") + ylab("FRiP") +
  theme(legend.position='top')
p4 <- ggplot(meta.a, aes(x = sharing, y = pOrg, fill = library)) + 
  geom_violin(position = position_dodge(0.7), trim = FALSE) + 
  geom_boxplot(position = position_dodge(0.7), color = "white", 
               width = 0.05, show.legend = FALSE) + theme_classic()+xlab("") + ylab("Proportion of organelle") +
  theme(legend.position='top')
wrap_plots(p1, p2, p3, p4, nrow=1)
dev.off()

#save updated meta data.
write.table(meta.r, file=paste(prefix,".LabelShared.cells.RNA.metadata.txt", sep = ""), sep="\t", quote=F, row.names=T, col.names=T)
write.table(meta.a, file=paste(prefix,".LabelShared.cells.ATAC.metadata.txt", sep = ""), sep="\t", quote=F, row.names=T, col.names=T)

#df <- meta.r[,c("ATAC_cellID","celltype")]
#write.table(df, file=paste0(prefix,"_ATACShared_cell_clusters.txt"), row.names=F, col.names = F, sep = "\t", quote = F)

df <- meta.r[grepl("Shared",meta.r$sharing),c("ATAC_cellID","celltype")]
write.table(df, file=paste0(prefix,"_ATACShared_cell_clusters.txt"), row.names=F, col.names = F, sep = "\t", quote = F)



# #save shared cell only
# #RNA
# meta.r <- meta.r[meta.r$sharing == "Share",]
# rds.r <- readRDS(file=paste0(rna_path,'/',prefix,".raw.seurat.rds"))
# count.r <- GetAssayData(rds.r, slot="counts", assay="RNA")
# count.r <- count.r[,colnames(count.r) %in% rownames(meta.r)]
# count.r <- count.r[,rownames(meta.r)]
# #update cell barcode id in meta and count.r
# rownames(meta.r) <- paste(meta.r$tenx_CB, meta.r$library, sep = "-")
# colnames(count.r) <- rownames(meta.r)
# saveRDS(count.r, file = paste0(prefix,".SharedOnly.cells.RNA.sparse.rds"))
# write.table(meta.r, file=paste(prefix,".SharedOnly.cells.RNA.metadata.txt", sep = ""), sep="\t", quote=F, row.names=T, col.names=T)

# #ATAC and gene activity
# meta.a <- meta.a[meta.a$sharing == "Share",]
# count.a <- readRDS(file=paste0(atac_path,'/',prefix,".raw.soc.rds"))
# count.a <- count.a$count
# gene.a <- readRDS(gene_sparse)
# count.a <- count.a[,colnames(count.a) %in% rownames(meta.a)]
# gene.a <- gene.a[,colnames(gene.a) %in% rownames(meta.a)]
# count.a <- count.a[,rownames(meta.a)]
# #gene.a <- gene.a[,rownames(meta.a)]

# #update cell barcode id in meta and count.r
# rownames(meta.a) <- paste(meta.a$tenx_CB, meta.a$library, sep = "-")
# meta.a$raw_cellID <- meta.a$cellID
# meta.a$cellID <- rownames(meta.a)
# #meta.a <- meta.a[rownames(meta.r),]
# colnames(count.a) <- rownames(meta.a)
# colnames(gene.a) <- gsub("CB:Z:","",colnames(gene.a))

# saveRDS(count.a, file = paste0(prefix,".SharedOnly.cells.ATAC.sparse.rds"))
# saveRDS(gene.a, file = paste0(prefix,".SharedOnly.cells.ATAC.gene.sparse.rds"))
# write.table(meta.a, file=paste(prefix,".SharedOnly.cells.ATAC.metadata.txt", sep = ""), sep="\t", quote=F, row.names=T, col.names=T)

#ATAC gene sparse.




# #update cell barcode id in meta and count.r
# rownames(meta.a) <- paste(meta.a$tenx_CB, meta.a$library, sep = "-")
# #meta.a <- meta.a[rownames(meta.r),]
# colnames(gene.a) <- rownames(meta.a)
# saveRDS(gene.a, file = paste0(prefix,".SharedOnly.cells.ATAC.sparse.rds"))
# write.table(meta.a, file=paste(prefix,".SharedOnly.cells.ATAC.metadata.txt", sep = ""), sep="\t", quote=F, row.names=T, col.names=T)
# 

