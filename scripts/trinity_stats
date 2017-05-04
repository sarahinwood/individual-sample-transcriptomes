#!/usr/bin/env bash

bin/trinity/util/TrinityStats.pl \
	output/trinity/Trinity.fasta \
	> output/trinity_stats/stats.txt

bin/trinity/util/misc/contig_ExN50_statistic.pl \
	output/trinity_abundance/RSEM.TPM.not_cross_norm \
	output/trinity/Trinity.fasta \
	>output/trinity_stats/xn50.out.txt \
	2>output/trinity_stats/xn50.err.txt