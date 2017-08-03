#!/usr/bin/env bash

set -eu

bbmerge="bin/bbmap/bbmerge.sh"

thorax_r1="output/bbduk_trim/thorax_r1.fq.gz"
thorax_r2="output/bbduk_trim/thorax_r2.fq.gz"

abdo_r1="output/bbduk_trim/abdo_r1.fq.gz"
abdo_r2="output/bbduk_trim/abdo_r2.fq.gz"

"${bbmerge}" in="${thorax_r1}" in2="${thorax_r2}" \
	out=output/bbmerge/thorax_merged.fq.gz \
	outu1=output/bbmerge/thorax_r1_unmerged.fq.gz \
	outu2=output/bbmerge/thorax_r2_unmerged.fq.gz \
	ihist=output/bbmerge/thorax_ihist.txt \
	verystrict=t \
	adapters=bin/bbmap/resources/adapters.fa \
	&> output/bbmerge/thorax.log.txt &

"${bbmerge}" in="${abdo_r1}" in2="${abdo_r2}" \
	out=output/bbmerge/abdo_merged.fq.gz \
	outu1=output/bbmerge/abdo_r1_unmerged.fq.gz \
	outu2=output/bbmerge/abdo_r2_unmerged.fq.gz \
	ihist=output/bbmerge/abdo_ihist.txt \
	verystrict=t \
	adapters=bin/bbmap/resources/adapters.fa \
	&> output/bbmerge/abdo.log.txt &

wait

echo "done"

exit 0