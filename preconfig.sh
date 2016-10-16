
# ===================================================================
# REQUIRED INPUT                                                  
# ===================================================================

# -----------------------------------------------------------------------------
# GFF_DIR is a directory containing a GFF file for each species used in the
# pipeline. This GFF file must contain at minimum mRNA and coding sequence
# (CDS) features. The last column must contain a unique id for the specific
# gene model (mRNA). All start and stop positions must be relative to the
# reference genomes in FNA_DIR (see argument -n).
# 
# Chr1   .   mRNA   3631   5899   .   +   .   AT1G01010.1 Chr1   .   CDS 3760
# 3913   .   +   .   AT1G01010.1 Chr1   .   CDS    3996   4276   .   + .
# AT1G01010.1
# -----------------------------------------------------------------------------
GFF_DIR=~/db/fagin-input/gff


# -----------------------------------------------------------------------------
# FNA_DIR is a directory containing a single genome sequence file for each
# species used in the pipeline. The files must be in FASTA format.
# -----------------------------------------------------------------------------
FNA_DIR=~/db/fagin-input/fna


# -----------------------------------------------------------------------------
# syn_dir should be the name of directory containing one synteny map for each
# species that will be compared. each synteny map should consist of a single
# file named according to the pattern "<query>.vs.<target>.tab", for example,
# "arabidopsis_thaliana.vs.arabidopsis_lyrata.tab". these files should contain
# the following columns:
#
# 1. query contig name (e.g. chromosome or scaffold)
# 2. query start position
# 3. query stop position
# 4. target contig name
# 5. target start position
# 6. target stop position
# 7. score (not necessarily used)
# 8. strand relative to query
#
# here is an example:
#
# chr2	193631	201899	tchr2	193631	201899	100	+
# chr2	225899	235899	tchr2	201999	202999	100	+
# chr1	5999	6099	tchr1	6099	6199	100	+
# chr1	5999	6099	tchr1	8099	8199	100	+
# chr1	17714	18714	tchr2	17714	18714	100	+
# chr2	325899	335899	tchr2	301999	302999	100	+
#
# a synteny map like this can be created using a whole genome synteny program,
# such as satsuma (highly recommended). building a single synteny map requires
# hundreds of cpu hours, so it is best done on a cluster. an example pbs script
# is provided, see src/satsuma.pbs.
# -----------------------------------------------------------------------------
SYN_DIR=~/db/fagin-input/syn


# -----------------------------------------------------------------------------
# TREE is a newick format file specifying the topology of the species tree. It
# must contain all species used in the pipeline AND NO OTHERS (I may relax this
# restriction later).
#
# NOTE: There must be no spaces in the species names.
#
# Here is an example tree:
#
# (Brassica_rapa,(Capsella_rubella,(Arabidopsis_lyrata,Arabidopsis_thaliana)));
# -----------------------------------------------------------------------------
TREE=~/db/fagin-input/brassicaceae.tree


# -----------------------------------------------------------------------------
# FOCAL_SPECIES is the name of the one species whose orphans will be
# characterized by this pipeline (e.g. Arabidopsis_thaliana). The name must be
# consistent with the species names used elsewhere in the pipeline.
#
# For now, there can be only one focal species. Future releases may contain an
# all-vs-all option.
# -----------------------------------------------------------------------------
FOCAL_SPECIES=Arabidopsis_thaliana


# -----------------------------------------------------------------------------
# A list of the genes that will be analyzed. All the genes in the list must be
# represented by gene models of the same name in the focal species GFF file.
# -----------------------------------------------------------------------------
ORPHAN_LIST=~/db/fagin-input/orphan-list.txt




# ===================================================================
# The defaults below are usable for a first pass
# ===================================================================

# -----------------------------------------------------------------------------
# Minimum length of an interval in the synteny map (target side) --- This value
# can be adjusted to change the resolution of the map. Higher values reduce the
# number of syntenic links that will be included, thus generally increase the
# number queries contained entirely within the search intervals.  The cost is
# increased search space.
# -----------------------------------------------------------------------------
MINLEN=0

# -----------------------------------------------------------------------------
# --- Input start and stop bases for building the synder database
# This is a 4 character string of 0s and 1s. For building the database, only
# the first two matter. Options are ...
#   1. 0000 - input start and stop are 0-based
#   2. 1000 - input start is 1-based and stop is 0-based
#   3. 0100 - input start is 0-based and stop is 1-based
#   4. 1100 - input start and stop are both 1-based
# The correct choice dependes on the synteny builder you are using. Satsuma,
# the tool I use, seems to be 0-based relative to start positions and 1-based
# relative to stops. I have not been able to find details on their output, but
# the lowest stop positions equal 0 and the highest stops equal contig length
# (n, rather than n-1).
# -----------------------------------------------------------------------------
synder_db_bases=0100


# -----------------------------------------------------------------------------
# Stop and start base for the input and output to synder (after building the
# database). A properly formatted GFF files should be 1-based, according to the
# specs[1]. Fagin assumes 1-based output (required by Bioconductor and R in
# general). So the value should usually be the 1111.
#
# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
# -----------------------------------------------------------------------------
synder_search_bases=1111


# -----------------------------------------------------------------------------
# Synder v0.12.2+ can merge doubly-overlapping syntenic blocks and also skip a
# specified number of off-scaffold intervals that interrupt an otherwise valid
# contiguous set.
# -----------------------------------------------------------------------------
synder_k=10




# ===================================================================
# Automatically generated variables, do not change
# ===================================================================

# -----------------------------------------------------------------------------
# Fagin home directory containing
#  .
#  ├── bin
#  ├── dbg
#  ├── doc
#  ├── etc
#  ├── src
#  ├── configure.sh
#  ├── Makefile
#  ├── README.md
#  └── VERSION
# -----------------------------------------------------------------------------
HOME_DIR="/home/shoggoth/src/git/fagin"

# -----------------------------------------------------------------------------
# File to hold links to input data and intermediate files. This directory
# should be on a filesystem with adequate storage (+10Gb, depending on project
# size).
# -----------------------------------------------------------------------------
INPUT="/home/shoggoth/src/git/fagin/input"

# -----------------------------------------------------------------------------
# Default paths for automatically buildable dependencies. If the program
# already exists in $PATH, that value will take precedence and overwrite the
# value in /bin. More accurately, no program will be built in /bin if the
# program exists in $PATH already.
# -----------------------------------------------------------------------------

synder='/home/shoggoth/src/git/fagin/bin/synder'
smof='/home/shoggoth/src/git/fagin/bin/smof'
