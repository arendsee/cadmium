#!/usr/bin/env bash
set -u

usage (){
    echo "Filter best hits from a blast output file"
    echo "REQUIRED ARGUMENTS"
    echo " -i Input blast results file (output of blastp-vs-others.sh)"
    echo " -o Output rows from input with max hits for each query gene against each target species" 
    exit 0
}

# print help with no arguments
[[ $# -eq 0 ]] && usage

while getopts "hi:o:" opt; do
    case $opt in
        h)
            usage ;;
        i) 
            blastresults=$OPTARG ;;
        o)
            maxhits=$OPTARG ;;
    esac 
done

awk '
    BEGIN{OFS="\t"}
    {k=$2"\t"$3; if(a[k] < $6){a[k] = $6; s[k]=$0}}
    END{for (k in a) {print s[k]}}
' $blastresults > $maxhits

awk '
    !$3 in s {s[$3] = 0}
    $5 < 1e-5 {s[$3]++}
    END{for (k in s) {print k, s[k]}}
    ' blast-results-max.tab |
        awk '$2 < 3 {print $1}'
