#!/bin/bash
#SBATCH --job-name=snRNA_clustering               #Job name (testBowtie2)
#SBATCH --partition=iob_p       #Queue name (batch)
#SBATCH --nodes=1                 # Run all processes on a single node
#SBATCH --ntasks=1                #Run in a single task on a single node
#SBATCH --cpus-per-task=10         # Number of CPU cores per task (8)
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
#cellranger-atac mkref cellRanger --config=./Arabidopsis_TAIR10_config.txt
 module load R/4.4.2-gfbf-2024a
 module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
 ml R-bundle-CRAN/2024.11-foss-2024a


rds_path=/path/to/04_scifi_multiome/Zm_seedling/dmulti-RNA/step02_QC_and_clustering
meta_path=/path/to/04_scifi_multiome/Zm_seedling/dmulti-RNA/step02_QC_and_clustering
marker_path='/path/to/04_scifi_multiome/markers'
clean_meta="/path/to/04_scifi_multiome/Zm_seedling/dmulti-RNA/step03_integrate_clustering/Zm_All_seedling_integrated.rna.Singlet.metadata.txt"

#---For gsea annotation---
marker1="/path/to/markers/markers.filtered.v5.txt" #for DEG identification.
go="/path/to/GenomeRef/Zea_mays_v5/maize.B73.AGPv5.mappedGOterms.biologicalprocess.gmt"

tissue=Zm_seedling
prefix=Zm_All_seedling_integrated
 
  echo "Step01: clustering for tissue $tissue..."
  Rscript ./00_integrate_Clustering_RNA-rm0Recluster_Zm.R $meta_path $rds_path $tissue $clean_meta $prefix 
  
  rna_count=${prefix}.rna.sparse.rds
  rna_meta=${prefix}.rna.Singlet.metadata.txt
  rna_id=celltype
  rna_pca=${prefix}.rna.Harmony.txt

  ls $marker_path/*.txt |while read marker
  do
    marker_id=$(basename $marker .txt)
  
    rna_id="celltype"

    echo "Step02: plot RNA marker expression data in UMAP for tissue $prefix and marker $marker..."
    Rscript ./02_plot_marker_accessibility.R $rna_meta $rna_count $rna_pca $marker 16 ${prefix}_${marker_id} $rna_id

    echo "Step03: plot imputed RNA marker expression data  heatmap for tissue $prefix and marker $marker..."
    Rscript ./03_PlotClusterZscore_markers_imputed_heatmap.R $rna_meta ${prefix}_${marker_id}_imputed_data.rds $marker $rna_id ${prefix}_${marker_id}

  done
								
  echo "Step05: Identify DEGs for rna for $prefix..."
  Rscript ./05_pseudobulk_DAG_edgeR.R $rna_meta $rna_count $marker1 $rna_id ${prefix}.rna

  echo "Step06: Run rna GSEA for $prefix..."
  Rscript ./06_celltype_gsea.R ${prefix}.rna_DAG_pseudobulk.txt $go ${prefix}.rna 16

