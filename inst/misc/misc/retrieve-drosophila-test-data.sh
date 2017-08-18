#!/usr/bin/env bash
set -u

species='ananassae erecta grimshawi melanogaster mojavensis persimilis pseudoobscura sechellia simulans virilis willistoni yakuba'
base=ftp://ftp.flybase.net/genomes/12_species_analysis/genomes

mkdir drosophila
cd drosophila
for s in $species
do
    sbase=$base/Drosophila_$s/current
    gff=$s.gff
    faa=$s.faa
    fna=$s.fna
    wget -O $faa.gz $sbase/fasta/*all-translation*
    wget -O $gff.gz $sbase/gff/*all-filtered*
    gunzip *gz
    # These Drosophila GFF files contain the FASTA reference after the feature rows
    # So here I will extract this FASTA data into a dedicated *.fna file and remove
    # the sequence content from the GFF file.
    sed -n '/^>/,$p' $gff > $fna
    sed '/FASTA/,$d' $gff | sed '/###/d' > z
    mv z $gff
done
cd ..
