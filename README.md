[![Travis-CI Build Status](https://travis-ci.org/arendsee/fagin.svg?branch=dev)](https://travis-ci.org/arendsee/fagin)
[![Coverage Status](https://img.shields.io/codecov/c/github/arendsee/fagin/master.svg)](https://codecov.io/github/arendsee/fagin?branch=dev)

# Fagin

A pipeline for the classification of orphans into origin classes using a syntenic filter.

## Funding

This work is funded by the National Science Foundation grant:

[NSF-IOS 1546858](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1546858)
Orphan Genes: An Untapped Genetic Reservoir of Novel Traits

## Installation

```R
devtools::install_github("arendsee/fagin")
library(fagin)
```

## Dependencies

Currently `fagin` has no dependencies outside of R. It makes heavy use of
bioconductor (e.g. `Biostring`, `GenomicRanges`, and `GenomicFeatures`). It
also uses the rather experimental packages 'synder' and 'rmonad'.

## Input

The following is required

 - Phylogeny for all included species
 - Name of the focal species
 - Synteny map for the focal species versus each other species
 - For each species
   - GFF file (must at least include gene models)
   - Full genome (GFF reference)

## Documentation

Go [here](https://github.com/arendsee/fagin-case-studies) to see working case
studies that you can adapt for your own projects.

You can also check out the (under construction) wiki [here](https://github.com/arendsee/fagin/wiki).

## Configuration

To run and configure `fagin`, you need to set paths to your data in
a configuration object. The default configuration can be generated

```R
config()
```

This will need to be tailored to your specific needs. To run the full fagin analysis, call

```R
# Where con is your configuration object
run_fagin(con)
```

## Pipeline

 - Identify target genes that overlap the search space.
 - Search the query protein against the overlapping target gene's ORFs
 - Search the query gene DNA against the search interval DNA sequences
 - Predict ancestor states
