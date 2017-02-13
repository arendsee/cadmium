#!/usr/bin/env bash

set -o nounset
set -o pipefail

species=$1

source config
source shell-utils.sh

parse_script=$PWD/parse-gff.py

input_gff=$INPUT/gff/$species.gff
input_fna=$INPUT/fna/$species.fna
output_faa=$INPUT/trans-orf/$species.faa

echo "Writing $species to $output_faa" >&2

mkdir -p $INPUT/trans-orf

check-read $input_gff $0
check-read $input_fna $0

tmpdir=/tmp/2a-$species
kill-mkdir $tmpdir

nexons=$(cut -f3 $input_gff | grep -cx exon)
nnames=$(cut -f9 $input_gff | grep -cP '(Name=|ID=)')
nparent=$(cut -f9 $input_gff | grep -cP 'Parent=')
echo hi >&2

if [[ $nexons -eq 0 ]]
then
    echo "Error: GFF file for '$species' must contain exon entries" >&2
fi

if [[ $nnames -eq 0 ]]
then
    echo "Error: GFF file for '$species' must contain Name or ID tags in 9th column" >&2
fi

if [[ $nparent -eq 0 ]]
then
    echo "Error: GFF file for '$species' must contain Parent tags in 9th column" >&2
fi

$parse_script -s exon -r Name Parent -dmp $input_gff |
    sort -k10 -k4n |
    rename_for_bedtools 10 |
    cut -f1-9 |
    bedtools getfasta   \
        -fi $input_fna  \
        -bed /dev/stdin \
        -fo /dev/stdout \
        -name |
    sed 's/::.*//' |
    awk '
        $1 ~ /^>/ && $1 in a { next }
        {a[$1]++; print}
    ' |
    getorf -filter -find 1 -minsize 30 |
    $smof clean -s > $output_faa
