#!/usr/bin/env bash
set -eu

export PATH="/Volumes/archive/deardenlab/asw-transcriptome/bin:${PATH}"

##Build index
bin/bowtie2-build output/trinity/Trinity.fasta \
	output/bowtie2/Trinity.fasta.index \
	--threads 50

##Perform Alignmnet
bin/bowtie2 \
	-p 10 -q \
	--threads 50 \
	-x output/bowtie2/Trinity.fasta.index \
	-1 output/bbduk_trim/abdo_r1.fq.gz,output/bbduk_trim/thorax_r1.fq.gz \
	-2 output/bbduk_trim/abdo_r2.fq.gz,output/bbduk_trim/thorax_r2.fq.gz \
	2>&1 1> /dev/null | tee output/bowtie2/align_stats.txt