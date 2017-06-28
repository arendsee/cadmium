#!/usr/bin/env R
# CONFIGURATION TEMPLATE - DO NOT CHANGE

# This script will be sourced as an R file. So everything written in it must be
# legal R code. Strings must be quoted. You may add any R code you like here.



# ===================================================================
# REQUIRED INPUT
# ===================================================================

# The focal species for this analysis. This should be the query sequence used
# in building the synteny maps. The name MUST be quoted, e.g.
# "Arabidopsis_thaliana"
FOCAL_SPECIES = NULL




# ===================================================================
# The defaults below are usable for a first pass
# ===================================================================

# Base p-value cutoffs (these will be adjusted for multiple testing
# query protein versus target gene alignments
PROT2PROT_PVAL = 0.05

# query protein versus all SI translated ORFs
PROT2ALLORF_PVAL = 0.05

# query protein versus translated ORFs from spliced transcripts
PROT2TRANSORF_PVAL = 0.05

# query genes versus entire SI (nucleotide match)
DNA2DNA_PVAL = 0.05

# Number of simulations
PROT2PROT_NSIMS     = 1e3
PROT2ALLORF_NSIMS   = 1e3
PROT2TRANSORF_NSIMS = 1e3

# Maximum value of m*n that will be searched
DNA2DNA_MAXSPACE = 1e8

# Ratio of search interval to query interval below which an indel is called
INDEL_THRESHOLD = 0.25




# ===================================================================
# You shouldn't change anything below
# The correct configuration will automatically be produced
# ===================================================================

HOME          = "FAGIN_HOME"
FAA_DIR       = "FAGIN_HOME/input/faa"
GFF_DIR       = "FAGIN_HOME/input/gff"
SYN_DIR       = "FAGIN_HOME/input/syn"
GENE_DIR      = "FAGIN_HOME/input/gene"
GENOME_DIR    = "FAGIN_HOME/input/fna"
SI_DIR        = "FAGIN_HOME/input/maps"
SPECIES_FILE  = "FAGIN_HOME/input/species"
SCAFLEN       = 'FAGIN_HOME/input/stat/scaffold-lengths.tab'
NSTRINGS      = 'FAGIN_HOME/input/stat/nstrings.tab'
ORPHAN_LIST   = "FAGIN_HOME/input/orphan-list.txt"
KBCOMP        = 'FAGIN_HOME/input/stat/kb-composition.tab'
ORFGFF        = "FAGIN_HOME/input/orf-gff"
ORFFAA        = "FAGIN_HOME/input/orf-faa"
TRANS_ORF     = "FAGIN_HOME/input/trans-orf"
CACHE         = "FAGIN_HOME/cache"
DECISION_TREE = "FAGIN_HOME/etc/cds.yaml"
TREE          = "FAGIN_HOME/input/tree"
