set -u

usage (){
cat << EOF
Given a target interval, get the full DNA sequence and between-STOP aa sequences
REQUIRED ARGUMENTS
  -g A GFF file telling where the target interval(s) are
  -f A DNA fasta file containing whole target genome (must be synced with GFF)
  -q The query protein sequence
EOF
    exit 0
}

# print help with no arguments
[[ $# -eq 0 ]] && usage

while getopts "hq:g:f:" opt; do
    case $opt in
        h)
            usage ;;
        q)
            quefaa=$OPTARG ;;
        g) 
            gff=$OPTARG ;;
        f)
            fna=$OPTARG ;;
    esac 
done

seqbase=`basename $quefaa`
seqbase=$(sed -r 's/\..*//' <<< $seqbase)

i=0
while read line
do
    i=$(( i + 1 ))
    base=$PWD/${seqbase}_$i
    echo $line | tr ' ' '\t' > $base.gff
    bedtools getfasta -fi $fna -bed $base.gff -s -fo $base.fna
    getorf -filter -find 0 -methionine N < $base.fna  > $base.faa

    mkdir tmp-blast
    cd tmp-blast
    makeblastdb -in $base.faa -dbtype prot -out orf -title orf 
    cd ..

    blastp \
        -query $quefaa \
        -db tmp-blast/orf \
        -outfmt '6 qseqid sseqid evalue bitscore qstart qend sstart send' \
        > $base.blast 

    rm -rf tmp-blast
done < $gff
