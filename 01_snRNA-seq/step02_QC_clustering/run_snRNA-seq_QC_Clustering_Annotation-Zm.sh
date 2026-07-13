#!/bin/bash
#SBATCH --job-name=QC_update               #Job name (testBowtie2)
#SBATCH --partition=iob_p       #Queue name (batch)
#SBATCH --nodes=1                 # Run all processes on a single node
#SBATCH --ntasks=1                #Run in a single task on a single node
#SBATCH --cpus-per-task=10         # Number of CPU cores per task (8)
#SBATCH --mem=100G                 # Job memory limit (10 GB)
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
# module load R/4.3.1-foss-2022a
# export LC_ALL=C

 module load R/4.4.2-gfbf-2024a
 module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
 ml R-bundle-CRAN/2024.11-foss-2024a

alin_path='/path/to/04_scifi_multiome/Zm_seedling/dmulti-RNA/step01-align'
marker_path='/path/to/04_scifi_multiome/markers'

alin_path='/scratch/xz24199/04_scifi_multiome/Zm_seedling/dmulti-RNA/step01-align'

for id in "${alin_path}"/*rep2Solo.out; do

 prefix=${id##*/}
 prefix=${prefix%Solo.out}
	
 #----run QC----------------
 echo " - running prprocess QC analysis for $prefix ..." 
 Rscript ./fun_01_QC_preprocess_Zm.R $alin_path $prefix

 echo " - running filter snRNA meta for $prefix ..."
 Rscript ./fun_03_snRNA_meta_QC-xz.R ${prefix}.raw.metadata.txt $prefix

 #----run clustering----------------
 echo " - running process and clustering for $prefix ..."
 Rscript ./fun_04_Clustering_RNA_Zm.R $prefix

 #------run annotation-----------------
 rna_count=${prefix}.rna.sparse.rds
 rna_meta=${prefix}.rna.Singlet.Cluster.metadata.txt
 rna_id=celltype
 rna_pca=${prefix}.rna.PCA.txt


 ls $marker_path/*.txt |while read marker
 do
  marker_id=$(basename $marker .txt)

  echo "Step05: plot RNA marker expression data in UMAP for tissue $prefix and marker $marker..."
  Rscript ./fun_05_plot_marker_accessibility.R $rna_meta $rna_count $rna_pca $marker 16 ${prefix}_${marker_id} $rna_id

  echo "Step06: plot imputed RNA marker expression data  heatmap for tissue $prefix and marker $marker..."
  Rscript ./fun_06_PlotClusterZscore_markers_imputed_heatmap.R $rna_meta ${prefix}_${marker_id}_imputed_data.rds $marker $rna_id ${prefix}_${marker_id}

 done

done
