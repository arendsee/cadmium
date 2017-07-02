[![Travis-CI Build Status](https://travis-ci.org/arendsee/fagin.svg?branch=master)](https://travis-ci.org/arendsee/fagin)
[![Coverage Status](https://img.shields.io/codecov/c/github/arendsee/fagin/master.svg)](https://codecov.io/github/arendsee/fagin?branch=master)

# Fagin

A pipeline for the classification of orphans into origin classes using a syntenic filter.

# Installation

```R
devtools::install_github("arendsee/fagin")
library(fagin)
```

# Dependencies

Currently `fagin` has no dependencies outside of R. It makes heavy use of the
core libraries of bioconductor (`Biostring` and `GenomicRanges`). The only
really experimental (read unstable) dependency is `synder`.

# Input

The following is required

 - Phylogeny for all included species
 - Name of the focal species
 - Synteny map for the focal species versus each other species
 - For each species
   - GFF file (must at least include gene models)
   - Full genome (GFF reference)

# Configuration

To run and configure `fagin`, you need to set paths to your data in
a configuration object. The default configuration can be generated

```R
config()
```

and taylored to specific needs.

# Pipeline

 - Identify target genes that overlap the search space.
 - Search the query protein against the overlapping target gene's ORFs
 - Search the query gene DNA against the search interval DNA sequences
 - Predict ancestor states


# TODO

Content
 - [x] generalize from 'orphan' to 'query'
 - [ ] prepare report for each query
 - [ ] create gene pages
 - [ ] add orthology statistics
 - [ ] print classification tree in report
 - [ ] add RNA-seq input
 - [ ] implement non-quadratic alignment (e.g. HMMER or BLAST)

Implementation
 - [x] Merge with R synder version
 - [x] * ab initio refactor as pure R package
 - [x] * replace all shellscripts
 - [x] * add test suite
 - [ ] all data in RSQLite databases (constant memory)
 - [ ] parallelize everything (divide-analyze-recombine)
 - [ ] integrate BLAST orphan identification
 - [ ] integrate phylostratigraphy
 - [ ] toss knitr, modularize for interactive exploration
 - [ ] integrate with Trelliscope and datadr

Final Destination
 - [ ] Fagin will include everything needed for orphan analysis
 - [ ] Incorporate BLAST and manage its results (or something better than BLAST)
 - [ ] Phylostratigraphic contextualization using NCBI common tree
 - [ ] Syntenic map creation
 - [ ] Visualize the results genewise with Trelliscope
