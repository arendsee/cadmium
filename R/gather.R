load_species <- function(species_name, input){

  inputs_ <- funnel(
    spec=species_name,
    conf=input
  )

  dna_ <- inputs_ %*>%
    (function(s,i){
       gene_genome_filename(s, i@fna_dir)
    }) %>>% load_dna

  gff_ <- inputs_ %*>%
    (function(s,i){
       get_gff_filename(s, i@gff_dir)
    }) %>>% load_gff

  nstrings_    <- dna_ %>>% derive_nstring
  orfgff_      <- dna_ %>>% derive_orfgff
  orffaa_      <- list(dna_, orfgff_) %*>% derive_orffaa
  aa_          <- list(dna_, gff_) %*>% derive_aa
  trans_       <- list(dna_, gff_) %*>% derive_trans
  transorfgff_ <- trans_ %>>% derive_transorfgff
  transorffaa_ <- list(trans_, transorfgff_) %>>% derive_transorffaa

  specsum_ <- funnel(
    gff.summary         = gff_         %>>% summarize_gff,
    dna.summary         = dna_         %>>% summarize_dna,
    aa.summary          = aa_          %>>% summarize_faa,
    trans.summary       = trans_       %>>% summarize_dna,
    orfgff.summary      = orfgff_      %>>% summarize_granges,
    orffaa.summary      = orffaa_      %>>% summarize_faa,
    transorfgff.summary = transorfgff_ %>>% summarize_granges,
    transorffaa.summary = transorffaa_ %>>% summarize_faa,
    nstring.summary     = nstring_     %>>% summarize_nstring
  ) %*>%
    new("species_summaries")

  specfile_ <- funnel(
    gff.file         = gff_         %>>% to_cache( label="gff"         , group=species_name ),
    dna.file         = dna_         %>>% to_cache( label="dna"         , group=species_name ),
    aa.file          = aa_          %>>% to_cache( label="aa"          , group=species_name ),
    trans.file       = trans_       %>>% to_cache( label="trans"       , group=species_name ),
    orfgff.file      = orfgff_      %>>% to_cache( label="orfgff"      , group=species_name ),
    orffaa.file      = orffaa_      %>>% to_cache( label="orffaa"      , group=species_name ),
    transorfgff.file = transorfgff_ %>>% to_cache( label="transorfgff" , group=species_name ),
    transorffaa.file = transorffaa_ %>>% to_cache( label="transorffaa" , group=species_name ),
    nstring.file     = nstring_     %>>% to_cache( label="nstring"     , group=species_name )
  ) %*>%
    new("species_data_files")

  funnel(
    files = specfile_,
    summaries = specsum_
  ) %*>%
    new("species_meta")

}


load_species_with_cache <- function(species_name, ...){
  # Retrieve the cached data if present
  result <- from_cache(species_name, "species_meta")  
  if(is.null(result)){
    result <- load_species(species_name, ...) 
    to_cache(result, species_name, "species_meta")  
  }
  result
}


load_synmaps <- function(target_species, focal_species, syndir){
  funnel(
    target_species,
    focal_species,
    syndir
  ) %*>%
    get_synmap_filename %>>%
    load_synmap %>%
    funnel(
      synmap.file = . %>>% to_cache(label="synmap", group=target_species),
      synmap.summary = . %>>% summarize_syn
    ) %*>%
      new("synteny_meta")
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

  con_ <- as_monad(con)

  tree_ <- con_ %>>% load_tree(.@input@tree)

  species_names_ <- tree_ %>>% {
    
    "Extract species list from the species tree. The tree, rather than the
    input files, determines which species are used in the analysis."

    .$tip.label
  }

  focal_species_ = con_ %>>% {

    "The input focal species may have gaps, but for the internal one, we remove
    gaps."
    
    gsub(" ", "_", .@input@focal_species)
                            
  }

  species_meta_list_ <- funnel(con=con_, specs=species_names_) %*>%
    function(con, specs){

      "For each species, collect and cache all required data"

      ss <- BiocParallel::bplapply(specs, load_species_with_cache, con@input)
      names(ss) <- specs
      ss
    }

  synmap_meta_list_ <- funnel(con=con_, specs=species_names_, focal=focal_species_) %*>%
    function(con, specs, focal){

      "Load synteny map for focal species against each target species"

      ss <- lapply(
        setdiff(specs, focal_species),
        load_synmaps,
        focal_species = focal,
        syndir        = con@input@syn_dir
      )
      names(ss) <- setdiff(specs, focal)
      ss
    }

  queries_ <- con_ %>>% { load_queries(.@input@query_gene_list) }

  funnel(
    tree = tree_,
    focal_species = focal_species_, 
    queries = queries_,
    species = species_meta_list_, 
    synnaps = synmap_meta_list_
  ) %*>%
    new("derived_input")

}
