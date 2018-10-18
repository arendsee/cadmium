# v0.11.0

 * Update yeast query and control genes. Now the query genes are the
   Saccharomyces-specific genes (from `phylostratr`) and the controls are an
   equal number of random non-Saccharomyves-specific genes.
 * Update vignette to include more summary info
 * Fix missing newlines in warning messages
 * Standardize naming conventions (snakecase for functions)
 * Fix ORF identification with the ORFik package

# v0.9.0

 * Add complete yeast data set
 * Replace bisquares regression with L1 regression 
 * Increase rmonad granularity
 * Add battery of GFF tests and transforms before conversion to TxDb object
