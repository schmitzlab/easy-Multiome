#!/bin/bash
#SBATCH --job-name=QC_update               #Job name (testBowtie2)
#SBATCH --partition=iob_p       #Queue name (batch)
#SBATCH --nodes=1                 # Run all processes on a single node
#SBATCH --ntasks=1                #Run in a single task on a single node
#SBATCH --cpus-per-task=10         # Number of CPU cores per task (8)
#SBATCH --mem=200G                 # Job memory limit (10 GB)
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
#cellranger-atac mkref cellRanger --config=./Arabidopsis_TAIR10_config.txt
 module load R/4.4.2-gfbf-2024a
 module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
 ml R-bundle-CRAN/2024.11-foss-2024a
 module load MACS2/2.2.9.1-foss-2023a


	#---maize---
	ref=/path/to/Zm-B73-REFERENCE-NAM-5.0-mainChr_CpMt_chrom.sizes
	ann=/path/to/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_mainChr_MtCp.RmNoStandGene.gtf
	
	# input
	bed=/path/to/01_mapping/AH15AH17_eMultiATAC_Zm_seedling-rep1.mq30.tn5.bed.gz
	name=AH15AH17_eMultiATAC_Zm_seedling-rep1
	
	echo " - running Socrates scATAC-seq QC analysis for $name ..."
	# run 
	Rscript ./QC_scATAC_data.v2_minuzx_20230203.R $bed $name $ann $ref
	Rscript ./fun_filter_low_qual_cells.v2_mihuzx_20230709.R $name.raw.soc.rds $name
    Rscript ./fun_scATAC_meta_QC-xz.R ${name}.updated_metadata.txt $name

