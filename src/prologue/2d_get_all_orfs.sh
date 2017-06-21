#!/usr/bin/env bash

set -o nounset
set -o pipefail

species=$1

smof=smof
output_faa=z.faa
output_gff=z.gff
input_fna=$1

# source config
# source shell-utils.sh
#
# # Get open reading frames from each input genome
# input_fna=$INPUT/fna/$1.fna
# output_faa=$INPUT/orf-faa/$1.faa
# output_gff=$INPUT/orf-gff/$1.gff
#
# safe-mkdir $INPUT/orf-faa
# safe-mkdir $INPUT/orf-gff
#
# check-read $input_fna $0

$smof clean --reduce-header $input_fna |
   # Find all START STOP bound ORFs with 10+ AA
    getorf -filter -find 1 -minsize 30 |
    # Filter out all ORFs with unknown residues
    $smof grep -v -q X |
    # Pipe the protein sequence to a protein fasta file
    tee  $output_faa |
    # Select only header lines
    grep '>' |
    # Parse a header such as:
    # >scaffold_1_432765 [258 - 70] (REVERSE SENSE)
    perl -pe 's/>(\S+)_(\d+) \[(\d+) - (\d+)\](.*)/$1 $3 $4 $1_$2 $5/' |
    # Prepare GFF
    awk '
        BEGIN{OFS="\t"}
        {
            seq_name = $1
            uid = $4
            if($5 ~ /REVERSE/){
                strand = "-"
            } else {
                strand = "+"
            }
            if($2 < $3){
                start = $2
                stop  = $3
            } else {
                start = $3
                stop  = $2
            }
        }
        { print seq_name, ".", "ORF", start, stop, ".", strand, ".", uid }
    ' > $output_gff

# # Remove the extra gunk getorf appends to headers
# $smof clean -s $output_faa > $output_faa~; mv $output_faa~ $output_faa
