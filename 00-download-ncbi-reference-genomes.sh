#!/usr/bin/env bash

# mamba install ncbi-datasets
outdir="ncbi"
outfile="test.zip"
filename="${outdir}/${outfile}"

mkdir -p $outdir

accessions=(GCF_001922985.1 GCF_014211975.1 GCF_014217335.1 GCF_014218805.1 GCF_002157225.2 GCF_902141825.1)

datasets download genome accession \
    ${accessions[@]} \
    --include cds,protein,genome,seq-report \
    --filename $filename

unzip $filename -d ${filename/.zip/}
