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

All results of the analysis will be archived (in the folder `ARCHIVE` by
default, this can be customized in the config). The archive contains the
following Rdata files:

 1. `d1.Rda` - gathering and checking raw data

 2. `d2.Rda` - comparisons of target species to the focal species

 3. `d3.Rda` - building a feature table from the target/focal comparisons

 4. `d4.Rda` - labels and trees built from the feature table

 5. `d5.Rda` - results of superimposing labels on the phylogenetic tree

 6. `final_obj.Rda` -  the final rmonad object


## Pipeline

 - Identify target genes that overlap the search space.
 - Search the query protein against the overlapping target gene's ORFs
 - Search the query gene DNA against the search interval DNA sequences
 - Predict ancestor states


## TODO

Content
 - [x] generalize from 'orphan' to 'query'
 - [ ] add orthology statistics
 - [ ] add RNA-seq input
 - [ ] implement non-quadratic alignment (e.g. HMMER or BLAST)
 - [ ] account for phylostratigraphic bias (alla Moyer)
 - [ ] get proper ORF finder
 - [ ] additional classes
 - [ ] * transposon footprints (Ux)
 - [ ] * overprinting (Ox)
 - [ ] * transposon association (Ox)
 - [ ] * frameshift (Ox)

Implementation
 - [x] merge with R synder version
 - [x] * ab initio refactor as pure R package
 - [x] * replace all shellscripts
 - [x] * add test suite
 - [ ] parallelize everything (divide-analyze-recombine)
 - [x] toss knitr, modularize for interactive exploration

Final Destination
 - [ ] Handle phylostratigraphy
 - [ ] Handle synteny map creation (MUMmer?)
 - [ ] Incorporate BLAST and manage its results (or something better than BLAST)
 - [ ] Incorporate RNAseq data
 - [ ] Visualize the results genewise with Trelliscope
 - [ ] integrate `phylostratr`
 - [x] toss knitr, modularize for interactive exploration
 - [ ] integrate with Trelliscope and datadr
 - [ ] parallelize everything (divide-analyze-recombine)
 - [ ] add logs (this will need some `rmonad` modification)

Report
 - [ ] Better order for UNO secondary classes
 - [ ] Distinguish between scaffolds and chromosomes
 - [ ] Retrieve cytology info?
 - [ ] Add alignment stats - overlaps, length of target and query
 - [ ] prepare report for each query
 - [ ] create gene pages
