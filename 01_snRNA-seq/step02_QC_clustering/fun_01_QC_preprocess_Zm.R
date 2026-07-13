# module load R/4.3.1-foss-2022a
##########################################
#### QC the data & initial clustering ####
##########################################

##=============================== Load package ===============================
library(Seurat)
library(dplyr)
library(ggplot2)
library(monocle3)
#library(scCustomize)
#library(qs)   #quickly writing and reading any R object to and from disk.
library(psych)  ## for stats describeby()
library(DoubletFinder)
library(harmony)

args <- commandArgs(T)
if(length(args) != 2){stop("Rscript snRNA-seq_process_xz.R <path> <tissue>")}
path <- as.character(args[1])
tissue <- as.character(args[2])

#path <- "/scratch/xz24199/04_scifi_multiome/dmulti-RNA/step01-align"
#tissue <- "AH15_Zm_seedling-rep2"


#===============================Load the raw dataset===============================
# mv features.tsv genes.tsv
#data1 <- Read10X(data.dir = paste0(path,"/",tissue,"Solo.out/GeneFull/raw/"))
mtx <- paste0(path,"/",tissue,"Solo.out/GeneFull/raw/matrix.mtx") #used unique mapping reads.
cells <- paste0(path,"/",tissue,"Solo.out/GeneFull/raw/barcodes.tsv")
features <- paste0(path,"/",tissue,"Solo.out/GeneFull/raw/features.tsv")

counts <- ReadMtx(mtx = mtx, cells = cells, features = features)

# Initialize the Seurat object with the raw (non-normalized data).
obj <- CreateSeuratObject(counts = counts, project = tissue,  min.features = 10)#ziliang filter cells

# #=============================== QC===============================
# QC Note:
# important features:
# unique genes / Cell (too little -> broken; too many -> doublets)
# molecules / Cell
# Organelle reads%
# rule of thumb:
# A general rule of thumb when performing QC is to set thresholds for individual metrics to be as permissive as possible, and always consider the joint effects of these metrics.

## ===============================1 edit metadata for QC===============================
obj[["pOrg"]] <- PercentageFeatureSet(obj, pattern = "^Ze") / 100
obj[["pCp"]] <- PercentageFeatureSet(obj, pattern = "^ZemaC") / 100
obj[["pMt"]] <- PercentageFeatureSet(obj, pattern = "^ZeamM") / 100
#obj[["pNmt"]] <- PercentageFeatureSet(obj, pattern = "^Hetgly") / 100

# Add number of genes per UMI for each cell to metadata
# this socre shows the complexity of the RNA, Generally, we expect the novelty score to be above 0.80 for good quality cells.
#obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)#ziliang code might a bug
obj$GenesPerUMI <- obj$nFeature_RNA / obj$nCount_RNA

# Create metadata dataframe
mtdt <- obj@meta.data
mtdt$Tech <- "eMulti_RNA"
mtdt$cellID <- paste(rownames(mtdt),mtdt$orig.ident,mtdt$Tech,sep="-")
mtdt$cell_barcode <- rownames(mtdt)
#rownames(mtdt) <- mtdt$cellID 
# Rename columns
mtdt <- mtdt %>%
  dplyr::rename(library = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

mtdt$log10nGene <- log10(mtdt$nGene)
mtdt$log10nUMI <- log10(mtdt$nUMI)

#----Add splice rate to meta (Velocyto) -----
#calculate spliced count matix
mtx.am <- paste0(path,"/",tissue,"Solo.out/Velocyto/raw/ambiguous.mtx")
mtx.spliced <- paste0(path,"/",tissue,"Solo.out/Velocyto/raw/spliced.mtx")
mtx.unspliced <- paste0(path,"/",tissue,"Solo.out/Velocyto/raw/unspliced.mtx")
cells <- paste0(path,"/",tissue,"Solo.out/Velocyto/raw/barcodes.tsv")
features <- paste0(path,"/",tissue,"Solo.out/Velocyto/raw/features.tsv")

counts.am <- ReadMtx(mtx = mtx.am, cells = cells, features = features)
counts.spliced <- ReadMtx(mtx = mtx.spliced, cells = cells, features = features)
counts.unspliced <- ReadMtx(mtx = mtx.unspliced, cells = cells, features = features)

splice.res <- data.frame(cells=colnames(counts.am), ambiguous_n = colSums(counts.am), spliced_n = colSums(counts.spliced), unspliced_n =  colSums(counts.unspliced))
splice.res$total <- splice.res$ambiguous_n + splice.res$spliced_n + splice.res$unspliced_n
splice.res$spliced_rate <- splice.res$spliced_n/splice.res$total
splice.res$spliced_rate[is.nan(splice.res$spliced_rate)] <- 0

splice.res <- splice.res[splice.res$cells %in% rownames(mtdt),]
splice.res <- splice.res[rownames(mtdt),]
mtdt <- cbind(mtdt, splice.res)


#filter orgnelle genes.
#obj <- readRDS("AH12_Gm_Infected_root-rep1.EmptyDropFiltered.seurat.rds")
counts <- GetAssayData(obj, slot="counts", assay="RNA")   
counts <- counts[!grepl("^Ze",rownames(counts)),]
#counts <- counts[!grepl("^Hetgly",rownames(counts)),]

#update counts and meta
overlap <- intersect(rownames(mtdt), colnames(counts))
counts <- counts[,overlap]
mtdt <- mtdt[overlap,]
colnames(counts) <- mtdt$cellID
rownames(mtdt) <- mtdt$cellID

# create new Seurat object
obj.n <- CreateSeuratObject(counts=counts, project = tissue, meta.data=mtdt)
# Create .RData object to load at any time
# save(obj, file="objd_24sample_rawdata_sObj.RData")

saveRDS(obj.n, file=paste0(tissue,".raw.seurat.rds"))
write.table(mtdt, file=paste0(tissue,".raw.metadata.txt"), sep= "\t", quote=F)




