#!/usr/bin/env bash
set -u

usage (){
    echo "Blast focal species proteins against a database for each of the other species"
    echo "REQUIRED ARGUMENTS"
    echo "  -h display this help message"
    echo "  -i input directory which conains *.faa files"
    echo "  -o output TAB-delimited file (STDOUT by default)"
    echo "  -f focal species"
    echo "  -n number of cores to use"
    exit 0
}

# print help with no arguments
[[ $# -eq 0 ]] && usage

num_threads=1
outfile=/dev/stdout
while getopts "hi:o:f:n:" opt; do
    case $opt in
        h)
            usage ;;
        i) 
            idir=${OPTARG%%/} ;;
        o) 
            outfile=$OPTARG ;;
        f)
            focal_species=$OPTARG ;;
        n)
            num_threads=$OPTARG ;;
    esac 
done

blastdb=$idir/blastdb
[[ -d $blastdb ]] || mkdir $blastdb

echo -e "focal_species\ttarget_species\tqseqid\tsseqid\tevalue\tbitscore" > $outfile
for f in $idir/*.faa
do
    n=`basename $f`
    n=${n%%.faa}
    if [[ ! -r $blastdb/$n.pin ]]
    then
        makeblastdb -dbtype prot -in $f -title $n -out $blastdb/$n
    fi
    blastp \
        -query $idir/$focal_species.faa \
        -db $blastdb/$n \
        -evalue '0.1' \
        -outfmt '6 qseqid sseqid evalue bitscore' \
        -num_threads $num_threads |
        awk -v focal=$focal_species -v target=$n 'BEGIN{FS="\t";OFS="\t"}{print focal, target, $0}'
done >> $outfile
