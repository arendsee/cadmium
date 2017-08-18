#!/bin/bash

# -------------------------------------------------------------------
# Prepare a PBS file for submission to a cluster
#
# Example:
# $ bash build-pbs.sh melanogaster simulans
# Where the genomic fasta files data/melanogaster.fna and data/simulans.fna
# must exist.
# -------------------------------------------------------------------

set -u

cat satsuma.pbs |
    sed "s/QUERY/$1/" |
    sed "s/TARGET/$2/" > $1.vs.$2.pbs
