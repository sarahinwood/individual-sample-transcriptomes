#!/usr/bin/env bash

bbduk="bin/bbmap/bbduk.sh"

outdir="output/bbduk_trim"

thorax_r1="data/thorax_r1.fastq.gz"
thorax_r2="data/thorax_r2.fastq.gz"

abdo_r1="data/abdo_r1.fastq.gz"
abdo_r2="data/abdo_r2.fastq.gz"

"${bbduk}" in="${thorax_r1}" in2="${thorax_r2}" \
    out=output/bbduk_trim/thorax_r1.fq.gz \
    out2=output/bbduk_trim/thorax_r2.fq.gz \
    ref=bin/bbmap/resources/adapters.fa \
    ktrim=r k=23 mink=11 hdist=1 tpe tbo \
    qtrim=r trimq=15 \
    &> "${outdir}"/thorax.log.txt &

"${bbduk}" in="${abdo_r1}" in2="${abdo_r2}" \
    out=output/bbduk_trim/abdo_r1.fq.gz \
    out2=output/bbduk_trim/abdo_r2.fq.gz \
    ref=bin/bbmap/resources/adapters.fa \
    ktrim=r k=23 mink=11 hdist=1 tpe tbo \
    qtrim=r trimq=15 \
    &> "${outdir}"/abdo.log.txt &

wait

echo "done"

exit 0
