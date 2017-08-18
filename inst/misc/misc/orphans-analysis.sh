#!/usr/bin/env bash
set -u

# Start with a GFF file (for now, just A. thaliana, will generalize later)
gff=input/gff/Arabidopsis_thaliana.gff

# Search just the
dbs=input/maps/db

get-clean-gff () {
cat $gff |
    # Extract chromosome entries (dropping organelle genomes)
    grep -E '^Chr'  |
    # Extract the mRNA entries
    grep mRNA |
    # Remove Chr, leaving the chr number, needed since GFF input to synder
    # requires id, not name, of scaffolds
    sed 's/^Chr//'
}

get-orphans () {
    orp=~/research/DATASETS/strata-loci.tab
    grep thaliana $orp | cut -f2
}

get-clean-gff | grep -f <(get-orphans) > input/orphans.gff

mkdir -p input/si

for s in Arabidopsis_lyrata
do
    out=input/si/Arabidopsis_thaliana.vs.$s.si.txt
    gff=input/orphans.gff
    db=$dbs/Arabidopsis_thaliana_$s.txt
    synder -i $gff -s $db -c contig > $out
done
