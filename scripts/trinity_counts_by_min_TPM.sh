#!/usr/bin/env bash

bin/trinity/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl \
output/trinity_abundance/RSEM.TPM.not_cross_norm \
>output/trinity_abundance/RSEM.TPM.not_cross_norm.counts_by_min_TPM