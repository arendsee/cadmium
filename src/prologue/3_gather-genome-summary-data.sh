#!/usr/bin/env bash

set -o nounset
set -o errexit
set -o pipefail

source config
source shell-utils.sh

usage (){
cat << EOF
Generate summary files for a set of genomes

Creates the following files in the input/stat folder

1. scaffold-lengths.tab - records the length of each assembled scaffold
   [ species | scaffold | length ]
2. nstrings.tab - records positions of strings of N in assembly
   [ species | scaffold | start | stop ]
3. kb-composition.tab - records base composition of 10000 random 1000-base samples
   [ species | scaffold | strand | A | C | G | T | N ]  # order may vary

OPTIONAL ARGUMENTS
  -i Input directory for genome sequences (default: input/fna)
  -o Output directory (default: input/stat)
  -h print this help message

REQUIREMENTS
  $smof
  parallel
  awk
  sed
EOF
    exit 0
}


idir=$INPUT/fna
odir=$INPUT/stat

while getopts "hi:o:" opt; do
    case $opt in
        h)
            usage ;;
        i) 
            idir=${OPTARG%/} ;;
        o) 
            odir=${OPTARG%/} ;;
        ?)
            exit 1 ;;
    esac 
done

check-dir $idir $0
safe-mkdir $odir

make-header () {
    echo -e "$(tr ' ' '\t' <<< $@)"
}

scaflen=$odir/scaffold-lengths.tab
nstring=$odir/nstrings.tab
nuccomp=$odir/kb-composition.tab

# Find lengths of all scaffolds in the genome file (for a fully assembled
# genome, scaffolds will correspond to chromosomes)
#
# OUTPUT COLUMNS
# 1. species name
# 2. scaffold name
# 3. scaffold length 
write-scaffold-lengths () {
    make-header 'species scaffold length'
    ls $idir/*fna | parallel "$smof stat -q {} > $odir/{/}.tab "
    for j in $odir/*fna.tab
    do
        s=${j%.fna.tab}
        s=${s##*/}
        sed "s/^/$s\t/" $j
        rm $j
    done
}

# Find positions of runs of unknown bases
#
# OUTPUT COLUMNS
# 1. species name
# 2. scaffold name
# 3. n-string start
# 4. n-string stop
write-nstrings () {
    make-header 'species scaffold start stop'
    ls $idir/*fna | parallel "$smof grep -Poq --gff --gff-type {/.} 'N+' {}" |
        awk 'BEGIN{FS="\t"; OFS="\t"} {print $3, $1, $4, $5}'
}

write-nucleotide-composition () {
    [[ -r $scaflen ]] || write-scaffold-lengths
    for j in $idir/*.fna
    do
        s=${j%.fna}
        s=${s##*/}
        bedtools random  \
            -l 1000      \
            -n 10000     \
            -seed 123456 \
            -g <(awk -v s=$s -v OFS="\t" '$1 == s {print $2,$3}' $scaflen) |
        bedtools getfasta   \
            -s              \
            -fi $j          \
            -bed /dev/stdin \
            -fo /dev/stdout |
        sed -r "s/>([^:]+):([0-9]+)-[0-9]+\((.)\)/>$s:\1:\2:\3/"
    done                     |
        $smof clean -xru -t n |
        $smof stat -qc        |
        tr ':' "\t"          |
        # Expand header to accomadate the added columns
        sed "1s/seqid/species\tscaffold\tstart\tstrand/"
}

write-scaffold-lengths       > $scaflen
write-nstrings               > $nstring
write-nucleotide-composition > $nuccomp

exit 0
