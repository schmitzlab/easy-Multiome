#!/bin/bash
#SBATCH --job-name=STAR_mapping               #Job name (testBowtie2)
#SBATCH --partition=iob_p         #Queue name (batch)
#SBATCH --nodes=1                 # Run all processes on a single node
#SBATCH --ntasks=1                #Run in a single task on a single node
#SBATCH --cpus-per-task=30         # Number of CPU cores per task (8)
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


genome_index=/path/to/GenomeRef/Zea_mays_v5/V5_CpMt_gene/STAR_Zm_v5_scRNA-seq
gtf=/path/to/GenomeRef/Zea_mays_v5/V5_CpMt_gene/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_mainChr_MtCp.gtf
whitelist=/path/to/04_scifi_multiome/BC_while_list/737K-cratac-v1-96-well-barcode.txt 

file_base=/path/to/04_scifi_multiome/Zm_seedling/00_rawdata/rna/AH_17Xuan_eMultiome-S1-rna_S10_L004
prefix=AH15AH17_Zm_seedling-rep1


read1=${file_base}_R1_001.fastq.gz
read2=${file_base}_R3_001.fastq.gz
i5_idx=${file_base}_R2_001.fastq.gz


#reverse complement 10x barcode
module load SeqKit/2.5.1
i5_idx_rc=$(dirname $i5_idx)/$(basename $i5_idx .fastq.gz).RC.fastq.gz
seqkit seq -r -p $i5_idx -o $i5_idx_rc

#add barcode to read2
read2_BC=$(dirname $read2)/$(basename $read2 .fastq.gz).BC.fastq.gz
perl ./join_fastqs.pl $i5_idx_rc $read2 $read2_BC

# append 10x barcode to read name
module load UMI-tools/1.1.2-foss-2022a-Python-3.10.4

read1_bc_umi=$(dirname $read1)/$(basename $read1 .fastq.gz).BC.UMI.fastq.gz
read2_bc_umi=$(dirname $read2)/$(basename $read2 .fastq.gz).BC.UMI.fastq.gz
umi_tools extract --stdin $read1  --extract-method=regex \
		--read2-in=$read2_BC --bc-pattern2='(?P<cell_1>.{16})(?P<umi_1>.{10})(?P<discard_1>.{1})(?P<cell_2>.{11})(?P<discard_2>.{1})T{3}.*' \
		--stdout=$read1_bc_umi --read2-out=$read2_bc_umi

# add cell barcode and umi back to read2
read2_bc_umi_seq=$(dirname $read2)/$(basename $read2 .fastq.gz).BC.UMI.seq.fastq.gz
perl ./CB_to_read_name.pl $read2_bc_umi $read2_bc_umi_seq

#only trim reads1
module load cutadapt/4.5-GCCcore-11.3.0
read1_bc_umi_trim=$(dirname $read1)/$(basename $read1 .fastq.gz).BC.UMI.Trim.fastq.gz

cutadapt --nextseq-trim=1 --cores 16 \
-a "file:/home/xz24199/bin/Trimmomatic-0.39/adapters/NexteraPE-PE_TruSeq3-PE-2.fa" \
--output $read1_bc_umi_trim \
$read1_bc_umi
 
#trim both ends and output reads with minimum-length 20bp.
read1_bc_umi_trim2=$(dirname $read1)/$(basename $read1 .fastq.gz).BC.UMI.Trim2.fastq.gz
read2_bc_umi_seq_trim2=$(dirname $read2)/$(basename $read2 .fastq.gz).BC.UMI.seq.Trim2.fastq.gz

cutadapt --nextseq-trim=1 --cores 16 --minimum-length=20 \
-a "GGGGGGGG" \
-A "GGGGGGGG" \
--output $read1_bc_umi_trim2 --paired-output $read2_bc_umi_seq_trim2 \
$read1_bc_umi_trim $read2_bc_umi_seq


#maping with STARsolo
module load STAR/2.7.10b-GCC-11.3.0
echo "mapping sample $prefix"
STAR --runThreadN 30 --readFilesCommand zcat \
    --genomeDir $genome_index \
    --outFileNamePrefix $prefix --readFilesIn $read1_bc_umi_trim2 $read2_bc_umi_seq_trim2 \
    --sjdbGTFfile $gtf \
    --soloType CB_UMI_Simple \
	--soloCBstart 1 --soloCBlen 27 \
	--soloUMIstart 28 --soloUMIlen 10 \
	--soloStrand Forward \
	--soloCBwhitelist $whitelist \
	--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloUMIdedup 1MM_CR \
    --soloCellFilter EmptyDrops_CR \
    --soloFeatures Gene GeneFull SJ Velocyto \
    --soloMultiMappers PropUnique \
    --soloBarcodeReadLength 0 \
    --twopassMode Basic \
	--clip3pAdapterSeq AAAAAAAA \
	--clip3pAdapterMMp 0.1 \
    --outSAMtype BAM SortedByCoordinate \
	--outFilterMultimapNmax 1 \
    --limitBAMsortRAM 200000000000 \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM

