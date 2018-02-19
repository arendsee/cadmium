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

get_trans <- function(m){
  rmonad::funnel(
    x = m,
    transcripts = view(m, 'transcript') %>>%
                  GenomicFeatures::exonsBy(by="tx", use.names=TRUE)
  ) %*>%
  GenomicFeatures::extractTranscriptSeqs %>>%
  {

    "Ensure all transcripts are named. If any of the names are missing (NA),
    then the Biostrings::writeXStringSet will report an error as a note,
    which won't stop processing. By checking here I can stop analysis at the
    right time.
    
    If any transcripts do have missing names, then these are removed with a
    warning.
    "

    na_indices <- which(is.na(names(.)))

    if(length(na_indices) > 0){ msg <-
"%s of %s transcripts had missing names. This is not good. You should look into
the problem. For now, the offending transcripts have been removed."
      warning(sprintf(msg, length(na_indices), length(.)))
      . <- .[-na_indices] 
    }

    .

  } %>>% {

    "Print the transcripts to a temporary file"

    filepath <- file.path(con@archive, paste0(".", species_name, "_trans.fna"))
    Biostrings::writeXStringSet(., filepath=filepath)
    filepath

  } %>>% {

    "Build an indexed version of the genome"

    Rsamtools::indexFa(.)

    .

  } %>>% {

    "Get a reference to the genome"

    Rsamtools::FaFile(.)

  }
}

load_species <- function(species_name, con){

  "Generate, summarize and merge all derived data for one species. The only
  inputs are a genome and a GFF file of gene models."

  format_translation_warning <- make_format_translation_warning(species_name)

  # local cacher for current target species
  .cacher <- function(m, tag){
    cacher(m, c(tag, species_name))
  }

  .view <- function(m, tag){
    rmonad::view(m, c(tag, species_name))
  }

  .tag <- function(m, tag){
    rmonad::tag(m, c(tag, species_name))
  }

  # load genome
  get_genome_filename(species_name, dir=con@input@fna_dir) %>%
      .cacher('genome_filename') %>>%
      load_dna %>%
        .cacher('genome') %>>%
        summarize_dna %>%
        .cacher('summary_genome') %>%

    # make Seqinfo
    .view('genome') %>>%
      scanFa_trw %>>%
      {GenomeInfoDb::Seqinfo(
        seqnames   = names(.),
        seqlengths = IRanges::width(.),
        isCircular = NA,
        genome     = species_name
      )} %>%
      .cacher('seqinfo') %>%

    # nstring
    .view('genome') %>>%
      scanFa_trw %>>% derive_nstring %>%
      .cacher('nstring') %>>%
      summarize_nstring %>%
      .cacher('summary_nstring') %>%
    .view('genome') %__%

    # GFF
    get_gff_filename(species_name, dir=con@input@gff_dir) %>%
      .tag('gff_filename') %>%
      rmonad::funnel(
        filename = .,
        seqinfo_ = .view(., 'seqinfo')
      ) %*>% load_gene_models %>% .cacher("gff") %>>%
      summarize_gff %>% .cacher("summary_gff") %>%

    # ngORF GFF
    .view('genome') %>>%
      scanFa_trw %>>%
      derive_orfgff %>%
      .cacher('orfgff') %>>%
      summarize_granges %>%
      .cacher('summary_orfgff') %>%

    # ngORF faa
    .view('genome') %>%
      rmonad::funnel(
        dna = .,
        gff = .view(., 'orfgff')
      ) %*>%
      extractWithComplements %>>%
      {
        list(format_warnings=format_translation_warning)
        Biostrings::translate(., if.fuzzy.codon="solve")
      } %>%
      .cacher('orffaa') %>>%
      summarize_faa %>% .cacher("summary_orffaa")

    # .view('gff') %>>%
    #   GenomicFeatures::cdsBy(by="tx", use.names=TRUE) %>%
    #   .cacher('transcripts') %>>% {
    #
    #     "Find the phase of the first CDS in each model. This will be zero if the
    #     model is complete."
    #
    #     # exon_rank unfortunately is what is says it is: exon rank, not CDS rank.
    #     # So the code:
    #     # R> unlist(.)[mcols(unlist(.))$exon_rank == 1]$cds_name %>% as.integer
    #     # does not work. So instead I have to do the following, which is vastly slower:
    #     sapply(., function(x) GenomicRanges::mcols(x)$cds_name[1]) %>% as.integer
    #
    #   } %>_% {
    #
    #     "Warn if there are any incomplete models"
    #
    #     incomplete <- sum(. != 0)
    #     total <- length(.)
    #
    #     if(incomplete > 0){
    #       warning(sprintf(
    #       "%s of %s gene models in %s are incomplete (the phase of the first CDS is
    #       not equal to 0). This does not include partial gene models that happen to
    #       begin in-phase. Details are recorded in the species_summary field
    #       model_phases.",
    #       incomplete, total, species_name))
    #     }
    #
    #   } %>% .cacher('aa_model_phase') %>%
    # .view('transcript') %>>% {
    #
    #     "Abominable hack. To get around the GenomicFeatures issue reported
    #     [here](https://support.bioconductor.org/p/101245/), I encode the phase
    #     information in the CDS name. Here I extract that data and offset the
    #     reading frame as needed."
    #
    #     ori <- .
    #
    #     unlist(ori) %>%
    #       {
    #         meta <- GenomicRanges::mcols(.)
    #         meta$phase <- as.integer(meta$cds_name)
    #         GenomicRanges::start(.) <- ifelse(
    #           meta$exon_rank == 1,
    #           GenomicRanges::start(.) + meta$phase,
    #           GenomicRanges::start(.)
    #         )
    #         .
    #       } %>%
    #       relist(ori)
    #
    #   } %>%
    #   .cacher('transcript') %>%

    # # aa
    # .view('genome') %>%
    #   rmonad::funnel(
    #     x = .,
    #     transcripts = .view('transcript')
    #   ) %*>%
    # GenomicFeatures::extractTranscriptSeqs %>>%
    # {
    #   list(format_warnings=format_translation_warning)
    #   Biostrings::translate(., if.fuzzy.codon="solve")
    # } %>%
    # .cacher('aa') %>>%
    # summarize_faa %>% .cacher('summary_aa') %>%
    #
    # # aa model phase
    # .view('aa') %>%
    #   rmonad::funnel(phases = view(., 'aa_model_phase', aa = .) %*>% {
    #     # This should always be true. If it is not, there is a logical problem
    #     # in the code, not the data.
    #     stopifnot(length(aa) == length(phases))
    #     list(
    #       summary = factor(phases) %>% summary,
    #       incomplete_models = names(aa)[phases != 0]
    #     )
    #   } %>% .cacher('aa_model_phase') %>%
    #
    # # trans
    # .view('genome') %>% get_trans %>% cacher('trans')

  ### TODO: revive this
  # transorffaa_ <-
  #   rmonad::funnel(
  #     dna = trans_,
  #     gff = transorfgff_
  #   ) %*>%
  #   extractWithComplements %>>%
  #   {
  #     list(format_warnings=format_translation_warning)
  #     Biostrings::translate(., if.fuzzy.codon="solve")
  #   }

  ##### TODO: add these warnings in the appropriate places
  # {
  #
  #   "Warn if any proteins have internal stop codons"
  #
  #   stops <- .@aa.summary@has_internal_stop
  #   tab <- .@aa.summary@table
  #
  #   if(any(stops)){
  #     n <- sum(stops)
  #     total <- length(stops)
  #     offenders <- tab[stops, ]$seqids %>% as.character %>%
  #       paste0(collapse=", ")
  #     msg <- "%s of %s proteins in %s have internal stops. These amino acid sequences were constructed from the mRNA models specified in the GFF file: %s"
  #     warning(sprintf(msg, n, total, species_name, offenders))
  #   }
  #
  # } %>_% {
  #
  #   "Warn if names used in the protein file do not match those of the transcripts"
  #
  #   aaids <- .@aa.summary@table$seqids
  #   trids <- .@trans.summary@table$seqids
  #
  #   if(! setequal(aaids, trids)){
  #
  #     not_in_tr <- setdiff(aaids, trids)
  #     not_in_aa <- setdiff(trids, aaids)
  #
  #     msg <- paste(
  #       "Protein and transcript names in %s do not match:",
  #       "transcript ids missing in proteins: %s",
  #       "protein ids missing in transcript: %s",
  #       sep="\n"
  #     )
  #
  #     warning(sprintf(
  #       msg,
  #       species_name,
  #       paste(not_in_aa, collapse=", "),
  #       paste(not_in_tr, collapse=", ")
  #     ))
  #   }
  #
  # }

}


# load_synmap_meta <- function(tspec, fspec, syndir){
#
#   "Load the synteny map for the given focal and target species. Then cache the
#   synteny map and return the cached filename and a summary of the map as a
#   `synteny_meta` object."
#
#   target_species <- GenomeInfoDb::genome(tspec@seqinfo) %>% unique
#   focal_species <- GenomeInfoDb::genome(fspec@seqinfo) %>% unique
#
#   if(length(target_species) != 1  ||  length(focal_species) != 1){
#     stop("More than one species found in Seqinfo file, this should not happen.")
#   }
#
#   if(any(is.na(target_species)) || any(is.na(focal_species))){
#     stop("Species must be set in each Seqinfo object (must not be NA)")
#   }
#
#   synmap <- rmonad::funnel(
#     focal_name  = focal_species,
#     target_name = target_species,
#     dir         = syndir
#   ) %*>%
#     get_synmap_filename %v>%
#     synder::read_synmap(
#       seqinfo_a = fspec@seqinfo,
#       seqinfo_b = tspec@seqinfo
#     )
#
#   synmap %>>% summarize_syn %>% cacher('synmap', tspec) %__%
#     synmap %>% cacher('synmap_summary', tspec)
# }


#' Load primary data
#'
#' @export
primary_data <- function(con){

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

  m <- load_tree(con@input@tree) %>% cacher('tree') %>>%
    {

      "Extract species list from the species tree. The tree, rather than the
      input files, determines which species are used in the analysis."

      .$tip.label

    } %>_% {
      targets <- .
      con@input@focal_species %>% cacher('focal_species') %>>% {

        "Assert that the focal species is in the phylogenetic tree"
    
        if(! (. %in% targets)){
          msg <- "Focal species '%s' is not in the species tree [%s]"
          stop(sprintf(msg, ., paste(targets, collapse=", ")))
        }

        NULL
      }
    } %>>% {
  
      "Remove spaces in the species names."
  
      gsub(" ", "_", .)
  
    } %>% cacher('target_species') %>>% {

      "For each species, collect and cache all required data"

      ss <- lapply(., load_species, con=con)
      names(ss) <- .
      rmonad::combine(ss)
      'yolo'
    } %__%
    con@input@query_gene_list %>>% load_gene_list   %>% cacher("query_genes") %__%
    con@input@control_gene_list %>>% load_gene_list %>% cacher("control_genes")

  m

  # rmonad::funnel(
  #   syndir = con@input@syn_dir,
  #   fspec  = view(m, 'focal_species'),
  #   tspecs = view(m, 'target_species')
  # ) %*>% {
  #   "Load synteny map for focal species against each target species"
  #
  #   ss <- lapply(
  #     tspecs,
  #     load_synmap_meta,
  #     fspec  = fspec,
  #     syndir = syndir
  #   )
  #
  #   rmonad::combine(ss)
  # } %>>%
  #   # pass this on as an Rmonad object
  #   as_monad(lossy=FALSE)
}
