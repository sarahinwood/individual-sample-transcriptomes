#!/usr/bin/env python3

import pathlib2
import pandas
import os

#############
# FUNCTIONS #
#############

def find_read_files(read_dir):
#Make list of files
	path_generator = os.walk(read_dir, followlinks = True)
	my_files = list((dirpath, filenames)
		for (dirpath, dirname, filenames)
		in path_generator)
#Make new dictionary & populate with files (flowcell = key)
	my_fastq_files = {}
	for dirpath, filenames in my_files:
		for filename in filenames:
			if filename.endswith('.fastq.gz'):
				my_flowcell = pathlib2.Path(dirpath).name
				my_fastq = str(pathlib2.Path(dirpath,filename))
				if my_flowcell in my_fastq_files:
					my_fastq_files[my_flowcell].append(my_fastq)
				else:
					my_fastq_files[my_flowcell]= []
					my_fastq_files[my_flowcell].append(my_fastq)
	return(my_fastq_files)

def sample_name_to_fastq(wildcards):
	sample_row = sample_key[sample_key['Sample_name'] == wildcards.sample]
	sample_id = sample_row.iloc[-1]['OGF_sample_ID']
	sample_flowcell = sample_row.iloc[-1]['Flow_cell']
	sample_all_fastq = [x for x in all_fastq[sample_flowcell]
						if '-{}-'.format(sample_id) in x]
	sample_r1 = sorted(list(x for x in sample_all_fastq
							if '_R1_' in os.path.basename(x)))
	sample_r2 = sorted(list(x for x in sample_all_fastq
							if '_R2_' in os.path.basename(x)))
	return({'r1': sample_r1, 'r2': sample_r2})

###########
# GLOBALS #
###########

read_dir = 'data/reads'

sample_key_file = 'data/sample_key.csv'

bbduk_adapters = '/adapters.fa'

#containers
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
busco_container = 'shub://TomHarrop/singularity-containers:busco_3.0.2'
trinity_container = 'shub://TomHarrop/singularity-containers:trinity_2.6.6'

#########
# SETUP #
#########

# generate name to filename dictionary
all_fastq = find_read_files(read_dir)

sample_key = pandas.read_csv(sample_key_file)

all_samples = sorted(set(sample_key['Sample_name']))

#########
# RULES #
#########

rule target:
	input:
		'output/trinity/Trinity.fasta'
		'output/fastqc'



rule run_Trinity:
	input:
		left = expand('output/bbmerge/{sample}_all_r1.fq.gz', sample=all_samples),
		right = expand('output/bbmerge/{sample}_unmerged_r2.fq.gz', sample=all_samples)
	output:
		'output/trinity/Trinity.fasta',
		'output/trinity/Trinity.fasta.gene_trans_map'
	params:
		outdir = 'output/trinity'
		left = lambda wildcards, input: ','.join(sorted(set(input.left))),
		right = lambda wildcards, input: ','.join(sorted(set(input.left)))
	singularity:
		trinity_container
	threads:
		20
	log:
		'output/logs/trinity.log'
	shell:
		'Trinity '
		'--SS_lib_type RF '
		'--max_memory 300G '
		'--CPU {threads} '
		'--output {params.outdir} '
		'--left {params.left} '
		'--right {params.right} '
		'--seqType fq'

rule merge_all_r1_reads
	input:
		r1 = 'output/bbmerge/{sample}_unmerged_r1.fq.gz',
		merged = 'output/bbmerge/{sample}_merged.fq.gz'
	output:
		joined = 'output/bbmerge/{sample}_all_r1.fq.gz'
	shell:
		'cat {input.r1} {input.merged} > {output.joined}'

rule bbmerge:
	input:
		r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
		r2 = 'output/bbduk_trim/{sample}_r2.fq.gz'
	output:
		merged = 'output/bbmerge/{sample}_merged.fq.gz',
		unm1 = 'output/bbmerge/{sample}_unmerged_r1.fq.gz',
		unm2 = 'output/bbmerge/{sample}_unmerged_r2.fq.gz',
		ihist = 'output/bbmerge/{sample}_ihist.txt'
	params:
		adapters = bbduk_adapters
	log:
		'output/logs/bbduk_merge/{sample}.log'
	singularity:
		bbduk_container
	threads:
		20
	shell:
		'bbmerge.sh '
		'threads={threads} '
		'in={input.r1} '
		'in2={input.r2} '
		'out={output.merged} '
		'outu1={output.unm1} '
		'outu2={output.unm2} '
		'ihist={output.ihist} '
		'verystrict=t '
		'adapters={params.adapters} '
		'&> {log}'

rule fastqc:
	input:
		expand('output/bbduk_trim/{sample}_r{n}.fq.gz',
			sample=all_samples, n=[1,2])
	output:
		directory('output/fastqc')
	shell:
		'mkdir -p {output} ; '
		'fastqc --outdir {output} {input}'


rule bbduk_trim:
	input:
		r1 = 'output/joined/{sample}_r1.fq.gz',
		r2 = 'output/joined/{sample}_r2.fq.gz'
	output:
		r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
		r2 = 'output/bbduk_trim/{sample}_r2.fq.gz'
	params:
		adapters = bbduk_adapters
	log:
		'output/logs/bbduk_trim/{sample}.log'
	threads:
		20
	singularity:
		bbduk_container
	shell:
		'bbduk.sh '
		'in={input.r1} '
		'in2={input.r2} '
		'out={output.r1} '
		'out2={output.r2} '
		'ref={params.adapters} '
		'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 '
		'&> {log}'

rule cat_reads:
	input:
		unpack(sample_name_to_fastq)
	output:	
		r1 = temp('output/joined/{sample}_r1.fq.gz'),
		r2 = temp('output/joined/{sample}_r2.fq.gz')
	threads:
		1
	shell:
		'cat {input.r1} > {output.r1} & '
		'cat {input.r2} > {output.r2} & '
		'wait'







