#!/bin/bash
#SBATCH --job-name=call_peaks
#SBATCH --partition=iob_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=300G
#SBATCH --time=7-00:00:00
#SBATCH --export=ALL
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xz24199@uga.edu
set -e
set -u
set -o pipefail
cd $SLURM_SUBMIT_DIR

#load modules
module load R/4.4.2-gfbf-2024a
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
ml R-bundle-CRAN/2024.11-foss-2024a
module load MACS2/2.2.9.1-foss-2023a
module load parallel/20240722-GCCcore-13.3.0
#module load R/4.3.1-foss-2022a

###########input vars##########
#--------tissues--------

bam=/scratch/xz24199/04_scifi_multiome/Zm_seedling/dmulti-ATAC-chromap/bed/AH15AH17_eMultiATAC_Zm_seedling.tn5.mq30.sorted.bed.gz
rna_meta="/scratch/xz24199/04_scifi_multiome/Zm_seedling/dmulti-RNA/step03-2_rm0_reclustering/Zm_All_seedling_integrated.rna.Singlet.metadata.txt"
atac_path="/scratch/xz24199/04_scifi_multiome/Zm_seedling/dmulti-ATAC-chromap/02_QC"
tissue="AH15AH17_eMulti_Zm_seedling"

#------------------------
#Genome information.

ref=/path/to/Zm-B73-REFERENCE-NAM-5.0-mainChr_CpMt_sorted_chrom.sizes
genome_size=1200000000
threads=18

Rscript ./fun_01_snRNA_ATAC_sharing_meta_QC-xz.R $rna_meta $atac_path $tissue

#################################################
# function to merge cells given cluster and BAM #
#################################################
mergeCells(){
	
	# vars
	clust=$1
	bamf=$2
	clustf=$3
	tissue=$4


	# select cells from correct cluster
	awk -F'\t' -v cluster=$clust '$2==cluster' $clustf | cut -f1 - > ${tissue}.$clust.bc_IDs.txt
	zcat $bamf |grep -Ff ${tissue}.$clust.bc_IDs.txt - | cut -f1-4 - |grep -v 'ChrCp\|ChrMt\|USDA110' - > ${tissue}_$clust.bed

	rm ${tissue}.$clust.bc_IDs.txt

}
export -f mergeCells


#############################
# iterate over all clusters #
#############################

# function
iterateClusters(){

	# load parameters
	i=$1
	bam=$2
	clusters=$3
	tissue=$4
	ref=$5
	g_size=$6
#	q_value=$7
#	input=$6

	# make directory
	if [ ! -d $PWD/$tissue ]; then
		mkdir $PWD/$tissue
	fi

	# verbose
	echo "merging cells from cluster $i ..."
	mergeCells $i $bam $clusters $tissue

	# call peaks
	echo "calling peaks with macs2 cluster $i ..."
	readdepth=$( wc -l < ${tissue}_$i.bed)
	echo "total reads = $readdepth"
	if [ $readdepth -lt 10000000 ];then
		fdr=0.1
	else
		if [ $readdepth -lt 25000000 ];then
			fdr=0.05
		else
			if [ $readdepth -lt 50000000 ]; then
				fdr=0.025
			else
				if [ $readdepth -lt 100000000 ]; then
					fdr=0.01
				else
					fdr=0.001
				fi
			fi
		fi
	fi 
	echo "FDR set to $fdr for $tissue $i..."
	macs2 callpeak -t ${tissue}_$i.bed \
		-f BED \
        -g $g_size \
        --nomodel \
        --keep-dup all \
        --extsize 150 \
        --shift -75 \
		--qvalue $fdr \
        --outdir $PWD/$tissue \
		--bdg \
        -n ${tissue}.$i.macs2
	
	# un-corrected
	sort -k1,1 -k2,2n $PWD/$tissue/${tissue}.$i.macs2_treat_pileup.bdg |cleanBED.pl $ref $readdepth - > $PWD/$tissue/${tissue}.$i.macs2_treat_pileup.clean.bdg
	wigToBigWig -clip $PWD/$tissue/${tissue}.$i.macs2_treat_pileup.clean.bdg $ref $PWD/$tissue/${tissue}.$i.macs2.clean.bw

	# clean
	rm $PWD/$tissue/${tissue}.$i.macs2_treat_pileup.bdg
	rm $PWD/$tissue/${tissue}.$i.macs2_control_lambda.bdg
    
	# move bed files to tissue directory
	 mv ${tissue}_$i.bed $PWD/$tissue

}
export -f iterateClusters

# iterate over clusters
parallel -j $threads iterateClusters {1} $bam ${tissue}_ATACShared_cell_clusters.txt $tissue $ref $genome_size ::: $( cut -f2 ${tissue}_ATACShared_cell_clusters.txt | sort -k1,1n | uniq )


# create merged set of peaks
cd $tissue
adjustPeaks.sh $tissue
