#!/usr/bin/env bash

set -o nounset
set -o pipefail

source config
source shell-utils.sh

species=$(cat $INPUT/species)

usage(){
cat << EOF >&2
Arguments
 -h print this help message

Extract
 1. protein sequences for each species
 2. nucleotide sequences for each search interval
 3. all ORF GFF file
 4. all ORF translations
 5. all transcript ORFs
 6. all transcript ORF translations

REQUIRES:
    bedtools v2.23.0
    emboss transseq
    emboss getorf
    smof
EOF
    exit 0
}

while getopts "h" opt; do
    case $opt in
        h)
            usage ;;
    esac
done

parallel --will-cite "./2a_get_transcript_orfs.sh {} " ::: $species
parallel --will-cite "./2b_get_genes.sh           {} " ::: $species
parallel --will-cite "./2c_get_proteins.sh        {} " ::: $species
parallel --will-cite "./2d_get_all_orfs.sh        {} " ::: $species
