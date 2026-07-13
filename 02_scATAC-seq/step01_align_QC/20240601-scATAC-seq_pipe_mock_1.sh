#!/bin/bash

###################################
#parameter for cellranger-atac count
cellranger_ref=/path/to/index/Zm_cellranger_v2
fastq_path=/path/to/fastq/rep1
base=Zm_seedling-rep1

###################################
#parameter for process bam
chrom_info=/path/to/Zm_v5_chrom.sizes
threads=24
memory=500
qual=30

#run cell ranger
if [ ! -d "./01cellranger" ]; then
  mkdir ./01cellranger
fi

cd ./01cellranger
cellranger-atac count --id=$base \
                  --reference=$cellranger_ref \
                  --fastqs=$fastq_path \
                  --localcores=$threads \
                  --localmem=$memory
cd ..

#process bam
module load BWA/0.7.17-GCCcore-11.2.0
module load ucsc/443
module load picard/2.27.5-Java-15
module load SAMtools/1.16.1-GCC-11.3.0
module load Python/3.7.4-GCCcore-8.3.0 
module load HarfBuzz/4.2.1-GCCcore-11.3.0
module load FriBidi/1.0.12-GCCcore-11.3.0

if [ ! -d "./02process_bam" ]; then
  mkdir ./02process_bam
fi
cd 02process_bam

        echo "filtering bam based on mapQ..."
#		mkdir ${base}_bams
        echo "retaining only mapped reads ..."
	    samtools view -@ $threads -bhq $qual -f 3 ../01cellranger/${base}/outs/possorted_bam.bam > ${base}.mq${qual}.bam

        # run picard
        echo "removing dups - ${base} ..."
        java -Xmx120g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
        INPUT=${base}.mq${qual}.bam \
        OUTPUT=${base}.mq${qual}.rmdup.bam \
        METRICS_FILE=${base}.rmdup.metrics \
        BARCODE_TAG=CR \
        ASSUME_SORT_ORDER=coordinate \
        REMOVE_DUPLICATES=true \
        USE_JDK_DEFLATER=true \
        USE_JDK_INFLATER=true 
       # VALIDATION_STRINGENCY=LENIENT \
        
		#append sample id to barcode.
        echo "fixing barcode and filter Cp Mt reads..."
        perl ../fun_fixBC.pl ${base}.mq${qual}.rmdup.bam | samtools view -bhS - > ${base}.BC.mq${qual}.rmdup.bam
        samtools stats ${base}.BC.mq${qual}.rmdup.bam > ${base}.BC.mq${qual}.rmdup.bam.stats

        #make Tn5 bed files
        echo "making Tn5 bed files ..."
        perl ../fun_makeTn5bed.pl ${base}.BC.mq${qual}.rmdup.bam | sort -k1,1 -k2,2n - > ${base}.tn5.mq${qual}.bed
		#uniq ${base}.tn5.bed > ${base}.tn5.uniq.bed
		
		if [ -s ${base}.tn5.mq${qual}.bed ]
		then
		  rm ${base}.mq${qual}.bam
		  #rm ${base}.mq${qual}.hdr.bam
		  rm ${base}.mq${qual}.rmdup.bam
		  #rm ${base}.tn5.mq${qual}.50bpwindow.bed
		  #rm ${base}.tn5.mq${qual}.50bpwindow.bg
		fi
		
