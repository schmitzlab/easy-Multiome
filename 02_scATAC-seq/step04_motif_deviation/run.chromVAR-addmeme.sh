#!/bin/bash
#SBATCH --job-name=chromVAR               #Job name (testBowtie2)
#SBATCH --partition=iob_p         #Queue name (batch)
#SBATCH --nodes=1                 # Run all processes on a single node
#SBATCH --ntasks=1                #Run in a single task on a single node
#SBATCH --cpus-per-task=30         # Number of CPU cores per task (8)
#SBATCH --mem=500G                 # Job memory limit (10 GB)
#SBATCH --time=7-00:00:00            # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --export=ALL              # Do not load any users<U+0092> explicit environment variables
#SBATCH --output=%x_%j.out        # Standard output log, e.g., testBowtie2_1234.out
#SBATCH --error=%x_%j.err         # Standard error log, e.g., testBowtie2_1234.err
#SBATCH --mail-type=END,FAIL      # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=xz24199@uga.edu # Where to send mail
set -e
set -u
set -o pipefail
cd $SLURM_SUBMIT_DIR

# get modules
 module load R/4.4.2-gfbf-2024a
 module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
 ml R-bundle-CRAN/2024.11-foss-2024a

# variables
threads=30
fai=/path/to/Zm-B73-REFERENCE-NAM-5.0-mainChr_CpMt_chrom.sizes
datalist=maize_list.txt
meme='/path/to/motif/405_TFmeme'

cat $datalist |while read id
do
data=$(echo $id | tr -d "\n")
array=(${data/// })  #split with tab and space;
sparse=${array[0]}
meta=${array[1]}
peaks=${array[2]}
pca=${array[3]}

prefix=$(basename $peaks .bed)

# run

Rscript ./chromVAR_analysis_script.maize.addmeme-v2.R $threads $sparse $meta $peaks $fai $prefix $meme
Rscript ./project.motif.UMAP-v2.R $meta ${prefix}.motif.deviations.txt $pca ${prefix}.motif.deviations
Rscript ./project.motif.UMAP-v2.R $meta ${prefix}.motif.scores.txt $pca ${prefix}.motif.scores
done
