#!/usr/bin/env bash

/Volumes/BiocArchive/archive/deardenlab/asw-transcriptome/bin/trinity/util/abundance_estimates_to_matrix.pl \
 	--est_method RSEM \
 	--cross_sample_norm none \
 	--out_prefix output/trinity_abundance/RSEM \
 	output/trinity_abundance/RSEM.isoforms.results