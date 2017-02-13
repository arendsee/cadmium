#!/usr/bin/env bash

set -o nounset
set -o pipefail

species=$1

source config
source shell-utils.sh

parse_script=$PWD/parse-gff.py

# Prepare FASTA file of genes, regions potentially include UTRs and introns
input_gff=$INPUT/gff/$species.gff
input_fna=$INPUT/fna/$species.fna
output_fna=$INPUT/gene/$species.gene.fna

safe-mkdir $INPUT/gene

check-read $input_gff $0
check-read $input_fna $0

# select mRNA and reduce 9th column to Name
$parse_script -s mRNA -r Name -d $input_gff |
    rename_for_bedtools 9 |
    bedtools getfasta   \
        -fi $input_fna  \
        -name           \
        -bed /dev/stdin \
        -fo /dev/stdout |
    sed '/>/s/::.*//' |
    $smof clean > $output_fna
