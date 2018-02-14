#!/usr/bin/env bash

contigs=$1
bins=$2
binname=$3
outpath=$4


seqtk subseq ${contigs} <(grep "\<${binname}\>" ${bins} | cut -f 1) > ${outpath}


