#!/usr/bin/env bash

set -o nounset
set -o pipefail

source config
source shell-utils.sh

exit_status=0

usage (){
cat << EOF >&2
Run focal species and target synteny maps through Synder to get search intervals

Output is written into the folder input/maps

Files are formated as: [query species].vs.[target species].map.tab

REQUIRES:
    synder
EOF
    exit $exit_status
}

species=$(cat $INPUT/species)
syndir=$INPUT/syn
gffdir=$INPUT/gff
mapdir=$INPUT/maps
glenfile=$INPUT/stat/scaffold-lengths.tab

safe-mkdir $mapdir

check-dir  $syndir   $0
check-dir  $gffdir   $0
check-read $glenfile $0

while getopts "h" opt; do
    case $opt in
        h)
            usage ;;
    esac 
done

# Get genome lengths file
genlen () {
    awk -v s=$1 'BEGIN{OFS="\t"} NR > 1 && $1 == s {print $2, $3}' $glenfile
}

status-check () {
    if [[ $1 != 0 ]]; then
        echo -e "\e[1;31m$2\e[0m"
        exit_status=1
    fi
}

for s in $species
do
    if [[ $s != $FOCAL_SPECIES ]]
    then
        map=$mapdir/$FOCAL_SPECIES.vs.$s.map.tab
        log=/tmp/$s.log
        echo $s

        synfile="$syndir/$FOCAL_SPECIES.vs.$s.syn"
        tmpque=/tmp/que$RANDOM
        tmptar=/tmp/tar$RANDOM

        check-read  $synfile $0

        genlen       $s             > $tmptar
        genlen       $FOCAL_SPECIES > $tmpque

        time $synder search         \
            -s $synfile             \
            -i $INPUT/search.gff    \
            -t $tmptar              \
            -q $tmpque              \
            -b $synder_search_bases \
            -k $synder_k            \
            -x d > $map 2> $log
        status-check $? "  synder failed"

        echo "Log written to '$s'"

        rm $tmpque $tmptar
    fi
done

exit $exit_status
