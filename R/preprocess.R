# #' Load a nexus format tree
# #'
# #' @param treefile A filename for a nexus tree
# #' @return a phylo object
# #' @export
# load_tree_io <- function(treefile){
#   ape::read.tree(treefile)
# }
#
# #' Get species from tree
# #'
# #' @param tree phylo object representing all under inspection
# #' @return character vector of unique species
# #' @export
# get_species_from_tree <- function(tree){
#   tree$tip.label
# }

# write.species_data <- function(){
#
# }
#
# build_derived_inputs <- function(config){
#   tree <- load_tree_io(config@input@tree)
#   species <- get_species_from_tree(tree)
# }
