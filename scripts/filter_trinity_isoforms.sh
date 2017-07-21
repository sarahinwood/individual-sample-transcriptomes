#!/usr/bin/env bash

set -eu

bin/bbmap/filterbyname.sh \
in=output/trinity/Trinity.fasta \
include=t \
names=output/trinity_abundance/isoforms_by_expression.txt \
out=output/trinity_filtered_isoforms/isoforms_by_expression.fasta

bin/bbmap/filterbyname.sh \
in=output/trinity/Trinity.fasta \
include=t \
names=output/trinity_abundance/isoforms_by_length.txt \
out=output/trinity_filtered_isoforms/isoforms_by_length.fasta 