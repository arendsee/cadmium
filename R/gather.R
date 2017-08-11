# A wrapper for scanFa that sets the seqname to the first word in the header
# This is probably the behaviour the function should have, since this is done
# when subsetting a sequence using a GRanges object.
scanFa_trw <- function(x, ...){

  "Load an XStringSet object from an indexed FASTA file. Use the first word in
  the header as the sequence name"

  seq <- Rsamtools::scanFa(x, ...)
  names(seq) <- sub(' .*', '', names(seq)) 
  seq
}


load_species <- function(species_name, input){

  "Generate, summarize and merge all derived data for one species. The only
  inputs are a genome and a GFF file of gene models."

  txdb_ <-
    get_gff_filename(species_name, dir=input@gff_dir) %v>%
    {
      specname <- gsub(pattern='_', replacement=' ', x=species_name)
      GenomicFeatures::makeTxDbFromGFF(., organism=specname)
    }

  dna_ <-
    get_genome_filename(species_name, dir=input@fna_dir) %v>%
    load_dna

  seqinfo_ <- dna_ %>>% scanFa_trw %>>% {
    rmonad::funnel(
      seqnames   = names(.),
      seqlengths = IRanges::width(.),
      isCircular = NA,
      genome     = species_name
    )
  } %*>% GenomeInfoDb::Seqinfo


  nstrings_ <- dna_ %>>% scanFa_trw %>>% derive_nstring

  orfgff_ <- dna_ %>>% scanFa_trw %>>% derive_orfgff

  orffaa_ <-
    rmonad::funnel(
      dna = dna_,
      gff = orfgff_
    ) %*>%
    extractWithComplements %>>%
    Biostrings::translate(if.fuzzy.codon="solve")

  aa_ <-
    rmonad::funnel(
      x = dna_,
      transcripts =
        txdb_ %>>% GenomicFeatures::cdsBy(by="tx", use.names=TRUE)
    ) %*>%
    GenomicFeatures::extractTranscriptSeqs %>>%
    Biostrings::translate(if.fuzzy.codon="solve")

  trans_ <-
    rmonad::funnel(
      x = dna_,
      transcripts = txdb_ %>>% GenomicFeatures::exonsBy(by="tx", use.names=TRUE)
    ) %*>%
    GenomicFeatures::extractTranscriptSeqs %>>%
    {
      filepath <- paste0(".", species_name, "_trans.fna")
      Biostrings::writeXStringSet(., filepath=filepath)   
      Rsamtools::indexFa(filepath)
      Rsamtools::FaFile(filepath)
    }

  transorfgff_ <- trans_ %>>% scanFa_trw %>>% derive_orfgff

  transorffaa_ <-
    rmonad::funnel(
      dna = trans_,
      gff = transorfgff_
    ) %*>%
    extractWithComplements %>>%
    Biostrings::translate(if.fuzzy.codon="solve")

  specsum_ <- rmonad::funnel(
    gff.summary         = txdb_        %>>% summarize_gff,
    dna.summary         = dna_         %>>% summarize_dna,
    aa.summary          = aa_          %>>% summarize_faa,
    trans.summary       = trans_       %>>% summarize_dna,
    orfgff.summary      = orfgff_      %>>% summarize_granges,
    orffaa.summary      = orffaa_      %>>% summarize_faa,
    transorfgff.summary = transorfgff_ %>>% summarize_granges,
    transorffaa.summary = transorffaa_ %>>% summarize_faa,
    nstring.summary     = nstrings_    %>>% summarize_nstring
  ) %*>%
    new(Class="species_summaries") %>_%
  {

    "Warn if any proteins have internal stop codons"

    stops <- .@aa.summary@has_internal_stop
    tab <- .@aa.summary@table

    if(any(stops)){
      n <- sum(stops)
      total <- length(stops)
      offenders <- tab[stops, ]$seqids %>% as.character %>%
        paste0(collapse=", ")
      msg <- "%s of %s proteins have internal stops, offending genes: %s"
      warning(sprintf(msg, n, total, offenders))
    }

  } %>_% {

    "Warn if names used in the protein file do not match those of the transcripts"

    aaids <- .@aa.summary@table$seqids
    trids <- .@trans.summary@table$seqids

    if(! setequal(aaids, trids)){

      not_in_tr <- setdiff(aaids, trids)
      not_in_aa <- setdiff(trids, aaids)

      msg <- paste(
        "Protein and transcript names do not match:",
        "transcript ids missing in proteins: %s",
        "protein ids missing in transcript: %s",
        sep="\n"
      )

      warning(sprintf(
        msg,
        paste(not_in_aa, collapse=", "),
        paste(not_in_tr, collapse=", ")
      ))
    }

  }

  specfile_ <- rmonad::funnel(
    gff.file         = txdb_        %>>% to_cache( label="gff"         , group=species_name ),
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
    files     = specfile_,
    summaries = specsum_,
    seqinfo   = seqinfo_
  ) %*>%
    new(Class="species_meta")

}


load_synmap_meta <- function(tspec, fspec, syndir){

  "Load the synteny map for the given focal and target species. Then cache the
  synteny map and return the cached filename and a summary of the map as a
  `synteny_meta` object."

  target_species <- GenomeInfoDb::genome(tspec@seqinfo) %>% unique
  focal_species <- GenomeInfoDb::genome(fspec@seqinfo) %>% unique

  if(length(target_species) != 1  ||  length(focal_species) != 1){
    stop("More than one species found in Seqinfo file, this should not happen.")
  }

  if(any(is.na(target_species)) || any(is.na(focal_species))){
    stop("Species must be set in each Seqinfo object (must not be NA)")
  }

  rmonad::funnel(
    focal_name  = focal_species,
    target_name = target_species,
    dir         = syndir
  ) %*>%
    get_synmap_filename %v>%
    synder::read_synmap(
      seqinfo_a = fspec@seqinfo,
      seqinfo_b = tspec@seqinfo
    ) %>%
    {
      rmonad::funnel(
        synmap.file = . %>>% to_cache(label="synmap", group=target_species),
        synmap.summary = . %>>% summarize_syn
      )
    } %*>%
      new(Class="synteny_meta")
}


#' Load primary data from the minimal required inputs
#'
#' @param con The config object that provides paths to required data
#' @return derived_input object
#' @export
primary_data <- function(...) {

  # TODO: this is a temporary caching function, repeal and replace.

  if(file.exists('d.Rda')){
    load('d.Rda')
  } else {
    d <- .primary_data(...)
  }

  save(d, file='d.Rda')

  d

}

.primary_data <- function(con){

  "
  Derive secondary data from the minimal required inputs
  
  The required inputs to Fagin are genomes fasta files, GFFs, synteny maps, a
  species tree, and the focal species name. This function derives all required
  secondary data, but does no analysis specific to the query intervals (e.g.
  the orphan genes candidates).
  
  This function does not handle the full data validation. All individual
  pieces are validated, but not the relationships between them. This is the
  role of the `validate_derived_inputs` function.
  "

  # set as monad here, so other functions will uniquely map to it
  con_ <- rmonad::as_monad(con)

  tree_ <- con_ %>>% { load_tree(.@input@tree) }

  species_names_ <- tree_ %>>% {

    "Extract species list from the species tree. The tree, rather than the
    input files, determines which species are used in the analysis."

    .$tip.label

  } %>>% {

    "Remove spaces in the species names."

    gsub(" ", "_", .)

  }

  focal_species_ = con_ %>>% { .@input@focal_species } %>%
  rmonad::funnel(specs=species_names_) %*>% {
    
    "Assert that the focal species is in the phylogenetic tree"

    if(! (. %in% specs)){
      msg <- "Focal species '%s' is not in the species tree [%s]"
      stop(sprintf(msg, ., paste(specs, collapse=", ")))
    }

    .

  }

  species_meta_list_ <- con_ %>>% { .@input } %>%
    rmonad::funnel(specs=species_names_) %*>%
    {

      "For each species, collect and cache all required data"

      ss <- lapply(specs, load_species, .)
      names(ss) <- specs

      rmonad::combine(ss)

    }

  synmap_meta_list_ <-
    rmonad::funnel(
      syndir = con@input@syn_dir,
      fspec  = species_meta_list_ %>>% { .[[con@input@focal_species]] },
      tspecs = species_meta_list_ %>>% {.[names(.) != con@input@focal_species]}
    ) %*>%
    {

      "Load synteny map for focal species against each target species"

      ss <- lapply(
        tspecs,
        load_synmap_meta,
        fspec  = fspec,
        syndir = syndir
      )

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
