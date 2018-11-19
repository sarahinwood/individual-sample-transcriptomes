#!/usr/bin/env bash

set -eu

bin/bbmap/filterbyname.sh \
in=output/trinity/Trinity.fasta \
include=t \
names=output/length_filtered/transcript_ids_above_500.txt \
out=output/length_filtered/longest_transcripts_above_500.fasta