#!/usr/bin/Rscript --vanilla

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("ape"))

parser <- ArgumentParser(
  formatter_class='argparse.RawTextHelpFormatter',
  description='Get species names (leaf labels) from a phylogenetic tree',
  usage='get-species-from-tree.R <filename>')

parser$add_argument(
  'tree',
  help='Nexus format tree'
)

args <- parser$parse_args()

cat(read.tree(args$tree)$tip.label, sep="\n")
