translate <- function(dna, species_name){

  "_ :: DNAStringSet -> AAStringSet   -- Species for warning annotation

  This function handles ambiguous bases, translating them into amino acids when
  possible, otherwise leaving them as X. 

  Warnings will be raised if the CDS is not a multiple of 3.
  "

  list(format_warnings=make_format_translation_warning(species_name))

  Biostrings::translate(dna, if.fuzzy.codon="solve")
}


filter_with_warnings__zero_length_proteins <- function(aa, species){

  "_ :: AAStringSet -> AAStringSet   -- Species for warning annotation
  
  Check for zero length proteins. If any are found, raise and warning and
  then remove them."

  if(any(GenomicRanges::width(aa) == 0)){
    zero_width_models <- names(aa)[GenomicRanges::width(aa) == 0]
    bad <- length(zero_width_models)
    total <- length(aa)
    msg <- paste("%s or %s mRNAs code for a 0 length protein. This is bad.",
                 "The following mRNA IDs are being removed: %s")
    warning(sprintf(msg, bad, total, paste(zero_width_models, collapse=", ")))
  }

  aa[GenomicRanges::width(aa) > 0]
}


trim_CDS_with_non_zero_phase <- function(grlist){

  "_ :: GRangesList -> GRangesList
  
  Offset start based on phase.
  
  WARNING: Abominable hack. To get around the GenomicFeatures issue
  reported [here](https://support.bioconductor.org/p/101245/), I encode
  the phase information in the CDS name. Here I extract that data and
  offset the reading frame as needed."

  gr <- unlist(grlist)

  meta <- GenomicRanges::mcols(gr)

  # Assert that the CDS name is in the set {0,1,2}
  if(! setequal(meta$cds_name, 0:2)){
    stop("Expected the `cds_name` field to be overloaded with phase info",
         "[0,1,2]. But it is not. This is a bug in the code.")
  }

  meta$phase <- as.integer(meta$cds_name)
  meta <- meta %>% as.data.frame %>%
    dplyr::group_by(cds_id) %>%
    dplyr::mutate(last_exon = exon_rank == max(exon_rank)) %>%
    as.data.frame

  GenomicRanges::start(gr) <- ifelse(
    meta$exon_rank == 1 & GenomicRanges::strand(.) == "+",
    GenomicRanges::start(gr) + meta$phase,
    GenomicRanges::start(gr)
  )
  GenomicRanges::end(gr) <- ifelse(
    meta$last_exon &  GenomicRanges::strand(gr) == "-",
    GenomicRanges::end(gr) - meta$phase,
    GenomicRanges::end(gr)
  )
  GenomicRanges::mcols(gr)$type <- "exon" # actually CDS, but
                                         # `extractTranscriptSeqs`
                                         # wants exons

  relist(gr, grlist)
}


convert_GRanges_to_SynderGFF <- function(gr){

  "GRanges -> synder::GFF"

  synder::as_gff(gr, id=GenomicRanges::mcols(gr)$seqid)
}


check_for_incomplete_models <- function(phases, species_name){

  "Warn if there are any incomplete models"

  if(! setequal(phases, 0:2)){
    stop("Expected phase information, with values from [0,1,2]")
  }

  incomplete <- sum(phases != 0)
  total <- length(phases)

  if(incomplete > 0){
    warning(sprintf(
    "%s of %s gene models in %s are incomplete (the phase of the first CDS is
    not equal to 0). This does not include partial gene models that happen to
    begin in-phase. Details are recorded in the species_summary field
    model_phases.",
    incomplete, total, species_name))
  }
}


extract_phase_from_GRangesList <- function(grlist){

  "GRangesList -> Phase  -- relies on dirty hack
  
  Find the phase of the first CDS in each model. This will be zero if the
  model is complete.

  WARNING: in load_gene_models, I set the 'cds_name' to 'phase'.  exon_rank
  unfortunately is what is says it is: exon rank, not CDS rank.
  "

  # NOTE: in load_gene_models, I set the 'cds_name' to 'phase'.
  # exon_rank unfortunately is what is says it is: exon rank, not CDS rank.
  # So the code:
  # R> unlist(.)[mcols(unlist(.))$exon_rank == 1]$cds_name %>% as.integer
  # does not work. So instead I have to do the following, which is vastly slower:
  phase <- sapply(., function(x) GenomicRanges::mcols(x)$cds_name[1]) %>% as.integer

  if(! setequal(phase, 0:2)){
    stop("Expected phase information, with values from [0,1,2]")
  }

  phase
}


check_for_internal_stops <- function(aa_summary){

  "Warn if any proteins have internal stop codons"

  stops <- aa_summary@has_internal_stop
  tab <- aa_summary@table

  if(any(stops)){
    n <- sum(stops)
    total <- length(stops)
    offenders <- tab[stops, ]$seqids %>% as.character %>%
      paste0(collapse=", ")
    msg <- "%s of %s proteins in %s have internal stops. These amino acid sequences were constructed from the mRNA models specified in the GFF file: %s"
    warning(sprintf(msg, n, total, species_name, offenders))
  }
}


check_protein_transcript_match <- function(aa_summary, trans_summary){

  "Warn if names used in the protein file do not match those of the transcripts"

  aaids <- aa_summary@table$seqids
  trids <- trans_summary@table$seqids

  if(! setequal(aaids, trids)){

    not_in_tr <- setdiff(aaids, trids)
    not_in_aa <- setdiff(trids, aaids)

    msg <- paste(
      "Protein and transcript names in %s do not match:",
      "transcript ids missing in proteins: %s",
      "protein ids missing in transcript: %s",
      sep="\n"
    )

    warning(sprintf(
      msg,
      species_name,
      paste(not_in_aa, collapse=", "),
      paste(not_in_tr, collapse=", ")
    ))
  }
}


m_get_proteins <- function(gffDB, genomeDB, species_name){

  "
  _ :: GFF -> Genome -> m ProteinSeq
  _ :: TxDb -> FaFile -> m AAStringSet
  
  Take a GFF database object (TxDb) and a genomeDB (FaFile) and build an
  rmonad graph containing the tagged fields 'faa' and 'summary_faa', the
  annotated protein sequences and the protein summary object, respectively.
  "

  .view <- function(m, tag) rmonad::view(m, tag, species_name)
  .tag  <- function(m, tag) rmonad::tag( m, tag, species_name)

  #- GeneModels ->
  #- GeneModelList {
  #-   cds_id    = Int,
  #-   cds_name  = String,
  #-   exon_rank = Int
  #- }
  gffDB %>>%
    GenomicFeatures::cdsBy(by="tx", use.names=TRUE) %>%
    .tag("cdsRangeList") %>>%
    #- GeneModelList -> [Phase]
    extract_phase_from_GRangesList %>_%
    #- [Phase] -> *Warning
    check_for_incomplete_models %>% .tag("aa_model_phase") %>%

  .view("cdsRangeList") %>>%
    #- GeneModelList -> GeneModelList
    trim_CDS_with_non_zero_phase %>>%
    #- GeneModelList -> GeneModelList
    filter_with_warning__unnamed_entries %>% .tag("cdsRanges") %>%
    #- Genome -> GeneModelList -> CDS
    rmonad::funnel(
      x = genomeDB,
      transcripts = .
    ) %*>%
    GenomicFeatures::extractTranscriptSeqs %>>%
    #- CDS -> AASeqs
    translate(species_name) %>% .tag("faa") %>>% 
    filter_with_warnings__zero_length_proteins(species_name) %>>%
    #- AASeqs -> AASummary
    summarize_faa %>% .tag("summary_aa") %>_%
    #- AASummary -> *Warning
    check_for_internal_stops
}
