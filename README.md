**This program is under heavy development, stability coming soon**

# Fagin

A pipeline for the classification of orphans into origin classes using a syntenic filter.

# Input

 The following is required

 - Phylogeny for all included species
 - Name of the focal species
 - Synteny map for the focal species versus each other species
 - For each species
   - GFF file (must at least include gene models)
   - Full genome (GFF reference)

# Assumptions about the input

 - The Parent tag in all GFFs for the CDS type maps uniquely to a protein. It
   doesn't have to be a protein ID, it may be a model id or transcript id. If
   multiple proteins have the same parent, the proteins will be merged
   incorrectly into one protein.
 - As above, an exon Parent tag must map uniquely to a transcript. This is less
   likely to be a problem.

# Pipeline

 - Identify target genes that overlap the search space.
 - Search the query protein against the overlapping target gene's ORFs
 - Search the query gene DNA against the search interval DNA sequences
 - Predict ancestor states

# Output

   A PDF file describing the results of the run.

 - The phylogenetic tree of all species in the pipeline
 - Genome summary for each species 
 - Overall statistics for classifications
 - Visualizations of overall statistics
 - Origin page for each orphan gene

# Running Fagin

 1. `./configure.sh` 
 2. Fill in generated config scripts
 3. `make load && make test && make run`

## 1. Run configure script

```
./configure.sh
```

This script will
  - check all dependencies, installing locally necessary
  - check all required R libraries, and install if possible
  - build basic config files

## 2. Setup config files

Fill in missing fields in `preconfig.sh` and `runconfig.R`. Details for each
required field are in the config files.

## 3. Run the analysis

```
make load && make test && make run
```


# TODO

Content
 - [ ] generalize from 'orphan' to 'query'
 - [ ] prepare report for each query
 - [ ] create gene pages
 - [ ] add orthology statistics
 - [ ] print feature table
 - [ ] print label table
 - [ ] print all intermediate data
 - [ ] only run functions specified by the function tree
 - [ ] print classification tree in report
 - [ ] add RNA-seq input
 - [ ] implement non-quadratic alignment

Implementation
 - [ ] Merge with R synder version
 - [ ] * ab initio refactor as pure R package
 - [ ] * replace all shellscripts
 - [ ] * full test coverage
 - [ ] all data in RSQLite databases (constant memory)
 - [ ] parallelize everything (divide-analyze-recombine)
 - [ ] integrate BLAST orphan identification
 - [ ] integrate phylostratigraphy
 - [ ] toss knitr, modularize for interactive exploration
 - [ ] integrate with Trelliscope and datadr

Final Destination
 - [ ] Fagin will include everything needed for orphan analysis
 - [ ] Identification of initial orphans using BLAST
 - [ ] Manage blast results (not a trivial manner)
 - [ ] Phylostratigraphic contextualization using NCBI common tree
 - [ ] Syntenic map creation
 - [ ] Synder interface
 - [ ] Search within syntenic context
 - [ ] Visualize the results genewise with Trelliscope
