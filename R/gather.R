load_species <- function(species_name, input){

  dna_ <- rmonad::funnel(species_name, dir=input@fna_dir) %*>%
    get_genome_filename %>>%
    load_dna

  gff_ <- rmonad::funnel(species_name, dir=input@gff_dir) %*>%
    get_gff_filename %>>%
    load_gene_models

  nstrings_    <- dna_ %>>% derive_nstring
  orfgff_      <- dna_ %>>% derive_orfgff
  orffaa_      <- list(dna_, orfgff_) %*>% derive_orffaa
  aa_          <- list(dna_, gff_) %*>% derive_aa
  trans_       <- list(dna_, gff_) %*>% derive_trans
  transorfgff_ <- trans_ %>>% derive_transorfgff

  transorffaa_ <-
  list(
    trans       = trans_,
    transorfgff = transorfgff_
  ) %*>% derive_transorffaa

  specsum_ <- rmonad::funnel(
    gff.summary         = gff_         %>>% summarize_gff,
    dna.summary         = dna_         %>>% summarize_dna,
    aa.summary          = aa_          %>>% summarize_faa,
    trans.summary       = trans_       %>>% summarize_dna,
    orfgff.summary      = orfgff_      %>>% summarize_granges,
    orffaa.summary      = orffaa_      %>>% summarize_faa,
    transorfgff.summary = transorfgff_ %>>% summarize_granges,
    transorffaa.summary = transorffaa_ %>>% summarize_faa,
    nstring.summary     = nstrings_    %>>% summarize_nstring
  ) %*>%
    new(Class="species_summaries")

  specfile_ <- rmonad::funnel(
    gff.file         = gff_         %>>% to_cache( label="gff"         , group=species_name ),
    dna.file         = dna_         %>>% to_cache( label="dna"         , group=species_name ),
    aa.file          = aa_          %>>% to_cache( label="aa"          , group=species_name ),
    trans.file       = trans_       %>>% to_cache( label="trans"       , group=species_name ),
    orfgff.file      = orfgff_      %>>% to_cache( label="orfgff"      , group=species_name ),
    orffaa.file      = orffaa_      %>>% to_cache( label="orffaa"      , group=species_name ),
    transorfgff.file = transorfgff_ %>>% to_cache( label="transorfgff" , group=species_name ),
    transorffaa.file = transorffaa_ %>>% to_cache( label="transorffaa" , group=species_name ),
    nstring.file     = nstrings_    %>>% to_cache( label="nstring"     , group=species_name )
  ) %*>%
    new(Class="species_data_files")

  rmonad::funnel(
    files = specfile_,
    summaries = specsum_
  ) %*>%
    new(Class="species_meta")

}


load_species_with_cache <- function(species_name, ...){
  # # Retrieve the cached data if present
  # result <- from_cache(species_name, "species_meta")
  # if(is.null(result)){
  #   result <- load_species(species_name, ...)
  #   to_cache(result, species_name, "species_meta")
  # }
  # result

  load_species(species_name, ...)
}


load_synmaps <- function(target_species, focal_species, syndir){
  rmonad::funnel(
    focal_name  = focal_species,
    target_name = target_species,
    dir         = syndir
  ) %*>%
    get_synmap_filename %v>%
    load_synmap %>%
    {
      rmonad::funnel(
        synmap.file = . %>>% to_cache(label="synmap", group=target_species),
        synmap.summary = . %>>% summarize_syn
      )
    } %*>%
      new(Class="synteny_meta")
}


#' Derive secondary data from the minimal required inputs
#'
#' The required inputs to Fagin are genomes fasta files, GFFs, synteny maps, a
#' species tree, and the focal species name. This function derives all required
#' secondary data, but does no analysis specific to the query intervals (e.g.
#' the orphan genes candidates).
#'
#' This function does not handle the full data validation. All individual
#' pieces are validated, but not the relationships between them. This is the
#' role of the \code{validate_derived_inputs} function.
#'
#' @param con The config object that provides paths to required data
#' @return derived_input object
#' @export
load_data <- function(con){

  # set as monad here, so other functions will uniquely map to it
  con_ <- rmonad::as_monad(con)

  tree_ <- con_ %>>% { load_tree(.@input@tree) }

  species_names_ <- tree_ %>>% {

    "Extract species list from the species tree. The tree, rather than the
    input files, determines which species are used in the analysis."

    .$tip.label
  }

  focal_species_ = con_ %>>% { .@input@focal_species } %>>% {

    "The input focal species may have gaps, but for the internal one, we remove
    gaps."

    gsub(" ", "_", .)
                            
  } %>% rmonad::funnel(specs=species_names_) %*>% {
    
    "Assert that the focal species is in the phylogenetic tree"

    if(! (. %in% specs)){
      msg <- "Focal species '%s' is not in the species tree [%s]"
      stop(sprintf(msg, ., paste(specs, collapse=", ")))
    }

    .

  }

  species_meta_list_ <- con_ %>>% { con@input } %>%
    rmonad::funnel(specs=species_names_) %*>%
    {

      "For each species, collect and cache all required data"

      ss <- BiocParallel::bplapply(specs, load_species_with_cache, .)
      names(ss) <- specs

      rmonad::combine(ss)

    }

  synmap_meta_list_ <- con_ %>>%
    { .@input@syn_dir } %>%
    rmonad::funnel(specs=species_names_, focal=focal_species_) %*>%
    {

      "Load synteny map for focal species against each target species"

      ss <- lapply(
        setdiff(specs, focal),
        function(x) {
          load_synmaps(
            target_species = x,
            focal_species = focal,
            syndir = .
          )
        }
      )
      names(ss) <- setdiff(specs, focal)

      rmonad::combine(ss)
    }

  queries_ <- con_ %>>% { load_queries(.@input@query_gene_list) }

  rmonad::funnel(
    tree          = tree_,
    focal_species = focal_species_,
    queries       = queries_,
    species       = species_meta_list_,
    synmaps       = synmap_meta_list_
  ) %*>%
    new(Class="derived_input")

}
