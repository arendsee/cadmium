#' Get the target species from a config object
#'
#' @export
get_targets <- function(con){
  setdiff(
    ape::read.tree(con@input@tree)$tip.label,
    con@input@focal_species
  ) %>% gsub(pattern=" ", replacement="_")
}

#' Get all species in the study
#'
#' @export
get_species <- function(con){
  ape::read.tree(con@input@tree)$tip.label %>%
    gsub(pattern=" ", replacement="_")
}
