###################################################################################################
###                             chromVAR analysis of scATAC data                                ###
###################################################################################################
#module load R/4.3.1-foss-2022a
#History: add function to read meme files to scan.
#History: addmeme-v2 use TFBStools

# arguments
args <- commandArgs(TRUE)
if(length(args) != 7){stop("Rscript chromVAR.R <threads> <sparseMatrix.rds> <meta> <peaks.bed> <fai> <output_prefix>")}
threads <- as.numeric(args[1])
input.sp <- as.character(args[2])
metadata <- as.character(args[3])
peakFile <- as.character(args[4])
FAI <- as.character(args[5])
prefix <- as.character(args[6])
memeDir <- as.character(args[7])

#load libraries
library(chromVAR)
library(motifmatchr)
library(BiocParallel)
#library(BSgenome.Gmax.a4.v1) 
library(BSgenome.maizeV5)
#BiocManager::install("BSgenome.Gmax.NCBI.Gmv40")
#library(BSgenome.Gmax.NCBI.Gmv40) #genome length is not the same
library(Matrix)
library(SummarizedExperiment)
library(GenomicAlignments)
library(dplyr)
library(JASPAR2022)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(TFBSTools)
library(universalmotif)


# functions
loadPeaks <- function(x, y, peaks, extra_cols=4){

	# create ref
	fai <- lapply(as.character(y$V1), function(z){
		return(as.numeric(as.character(y$V2)[y$V1==z]))
	})
	names(fai) <- as.character(y$V1)

	# load bed
	bed <- as.data.frame(do.call(rbind, strsplit(rownames(x),"_")))	
	newbed <- read.table(peaks)
	
	rownames(newbed) <- paste(newbed$V1,newbed$V2,newbed$V3,sep="_")
	#newbed <- subset(newbed, newbed$V4 != "exons")
	rownames(bed) <- rownames(x)
	bed <- bed[rownames(newbed),]
	rownames(bed) <- seq(1:nrow(bed))

	# convert 2 GR
	colnames(bed) <- c("chr", "start", "end")
	bed$chr <- as.character(bed$chr)
	bed$start <- as.numeric(as.character(bed$start))
	bed$end <- as.numeric(as.character(bed$end))
	bed$keep <- ifelse(bed$start > fai[bed$chr] | bed$end > fai[bed$chr], 0, 1)
	x <- x[bed$keep > 0, ]
	bed <- bed[bed$keep > 0,]
	bed$keep <- NULL
	bed[, "start"] <- bed[, "start"]
	bed <- makeGRangesFromDataFrame(bed, keep.extra.columns = F)

	# sort
	sorted_bed <- sortSeqlevels(bed)
	sorted_bed <- sort(sorted_bed, ignore.strand = TRUE)
	sbeddf <- as.data.frame(sorted_bed)
	s.ids <- paste(sbeddf$seqnames,sbeddf$start,sbeddf$end,sep="_")
	shared <- intersect(s.ids, rownames(x))
	x <- x[shared,]
	sorted_bed <- subset(sorted_bed, c(s.ids %in% shared))
	
	return(list(bed=sorted_bed, cnts=x))
}
getJasparMotifs2 <- function(species = "Homo sapiens", collection = "CORE", ...){
    opts <- list()
    opts["species"] <- species
    opts["collection"] <- collection
    opts <- c(opts, list(...))
    out <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022, opts)
    if (!isTRUE(all.equal(TFBSTools::name(out), names(out))))
        names(out) <- paste(names(out), TFBSTools::name(out),
            sep = "_")
    return(out)
}

#dir <- memeDir

loadMEMEDir <- function(dir){
#dir <- memeDir
    files <- list.files(
        dir,
        pattern="\\.meme$",
        full.names=TRUE
    )

    message(length(files), " MEME files found.")

    pwm.list <- vector("list", length(files))

    for(i in seq_along(files)){

        message("[", i, "/", length(files), "] ",
                basename(files[i]))

        ## read one motif
        umotif <- universalmotif::read_meme(files[i])

        ## convert one motif -> one PWMatrix
        pwm.list[[i]] <- universalmotif::convert_motifs(
            umotif,
            class="TFBSTools-PWMatrix"
        )
    }

    ## make a PWMatrixList
    names(pwm.list) <- sapply(pwm.list, name)
	return(pwm.list)
    #TFBSTools::PWMatrixList(pwm.list)
}

pwMatrixList_to_integerPWMatrix <- function(pwm_list,
                                            bg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                            scale = 1000) {

  stopifnot(is.list(pwm_list))

  convert_one <- function(pwm) {

    # extract matrix
    mat <- tryCatch({
      pwm@profileMatrix
    }, error = function(e) {
      pwm@matrix
    })

    # enforce A/C/G/T order
    mat <- mat[c("A", "C", "G", "T"), , drop = FALSE]

    # PWM (log2 odds) -> probabilities
    prob <- sweep(2^mat, 1, bg, `*`)

    # normalize columns
    ppm <- apply(prob, 2, function(x) x / sum(x))

    # scale to integer matrix
    int_mat <- round(ppm * scale)
    int_mat <- apply(int_mat, c(1, 2), as.integer)

    # rebuild PWMatrix (KEEP CLASS)
    new_pwm <- TFBSTools::PWMatrix(
      ID = pwm@ID,
      name = pwm@name,
      profileMatrix = int_mat,
      pseudocounts = pwm@pseudocounts,
      tags = pwm@tags
    )

    return(new_pwm)
  }

  lapply(pwm_list, convert_one)
}

# set number of cores
register(MulticoreParam(threads))

# verbose
message("########################################")
message("########################################")
message("")
message("============================")
message("     running chromVAR       ")
message("============================")
message("")


###################################################################################################
### load and process data									   
###################################################################################################

# build counts matrix
message("Loading count matrix ...")
a <- readRDS(input.sp)

atac_id <- sub("CB:Z:","",colnames(a)) #convert cell id to RNA cell id
parts <- do.call(
  rbind,
  strsplit(atac_id,"-",fixed = TRUE)
)
rna_id <- paste0(parts[,1],"AGTGATTAGCA-",parts[,2],"-",parts[,3],"-eMulti_RNA")
colnames(a) <- rna_id

# input files
message("Loading peak information ...")
ref <- read.table(FAI)
obj <- loadPeaks(a, ref, peaks=peakFile)
peaks <- obj$bed
a <- obj$cnts


# load meta.data
message("Loading meta data ...")
meta <- read.table(metadata,comment.char = "")
meta <- meta[meta$sharing == "Shared",]
a <- a[,colnames(a) %in% rownames(meta)]
meta <- meta[colnames(a),]
meta$depth <- Matrix::colSums(a)
#meta$depth <- colSums(a)
message("cells = ",ncol(a), " | peaks = ", nrow(a))

rownames(a) <- NULL

# create frag counts object
message("Creating experiment object ...")
fragment_counts <- SummarizedExperiment(assays = list(counts = a),
                                        rowRanges = peaks,
                                        colData = meta)

# clean-up memory
rm(a)
rm(obj)

#seqnames(BSgenome.Gmax.NCBI.Gmv40) <- c("Gm01","Gm02","Gm03","Gm04","Gm05","Gm06","Gm07","Gm08","Gm09","Gm10","Gm11","Gm12","Gm13","Gm14","Gm15","Gm16","Gm17","Gm18","Gm19","Gm20")

# add GC data
message("Estimating GC bias ...")
fragment_counts <- addGCBias(fragment_counts, genome=BSgenome.maizeV5)

# filter cells
message("Filtering samples ...")
filtered_counts <- filterSamples(fragment_counts, min_depth=100, min_in_peaks=0.1, shiny=F)

# filter peaks
message("Filtering peaks ...")
filtered_counts <- filterPeaks(filtered_counts, non_overlapping=T, min_fragments_per_peak=10) #Here I give it the final filtered peaks so setting a small number for min_fragments_per_peak 

###############################################################################
## motif deviation
###############################################################################

# estimate deviations
message("Running motif analysis ...")
#jaspar2016_motifs[['first']] <- first #add one object to motif list.
jaspmotifs.at       <- getJasparMotifs2(species = "Arabidopsis thaliana")
jaspmotifs.zm       <- getJasparMotifs2(species = "Zea mays")
jaspmotifs.os       <- getJasparMotifs2(species = "Oryza sativa")
jaspmotifs.gm       <- getJasparMotifs2(species = "Glycine max")

message("Loading MEME motifs...")

pwms <- loadMEMEDir(memeDir)
pwms2 <- pwMatrixList_to_integerPWMatrix(pwms)


message(length(pwms)," MEME motifs loaded.")

jaspmotifs <- c(jaspmotifs.at, jaspmotifs.zm, jaspmotifs.os, jaspmotifs.gm, pwms2)

motif            <- matchMotifs(jaspmotifs, filtered_counts, genome = BSgenome.maizeV5)
dev.motif        <- computeDeviations(object = filtered_counts, annotations = motif)
dev.motif.scores <- deviationScores(dev.motif)
motif.devs       <- deviations(dev.motif)
saveRDS(motif, file=paste0(prefix,".motif_matches.rds"))
write.table(t(dev.motif.scores), file=paste0(prefix,".motif.scores.txt"), quote=F, row.names=T, col.names=T, sep="\t")
write.table(t(motif.devs), file=paste0(prefix,".motif.deviations.txt"), quote=F, row.names=T, col.names=T, sep="\t")

## background peaks
bbpeaks <- getBackgroundPeaks(filtered_counts)
write.table(bbpeaks, file=paste0(prefix,".backgroundPeaks.mat.txt"), quote=F, row.names=T, col.names=T, sep="\t")
message("--Finished--")
