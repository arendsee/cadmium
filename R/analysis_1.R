# A wrapper for scanFa that sets the seqname to the first word in the header
# This is probably the behaviour the function should have, since this is done
# when subsetting a sequence using a GRanges object.
scanFa_trw <- function(x, ...){

  "Load an XStringSet object from an indexed FASTA file. Use the first word in
  the header as the sequence name"

  seq <- Rsamtools::scanFa(x, ...)
  names(seq) <- sub(" .*", "", names(seq)) 
  seq
}

get_trans_dna <- function(x, transcripts){
  GenomicFeatures::extractTranscriptSeqs(x, transcripts) %>>%
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

#- _ :: Genome -> GenomeSeq  -- with seqinfo
m_db2seq <- function(m, species_name){
  m %>>% scanFa_trw %>%
  {
    rmonad::funnel(
      seq = .,
      seqinfo_ = . %>>%
        #- _ :: Seqid -> SeqLength -> isCircular -> SpeciesName -> GenomeSeqinfo
        (function(seq){
          GenomeInfoDb::Seqinfo(
            seqnames   = names(seq),
            seqlengths = IRanges::width(seq),
            isCircular = NA,
            # The bioconductor people, in all their wisdom, decide to ignore
            # whatever name I give the species. Apparently, if it doesn't exist
            # in their web browser, it isn't worth storing.
            genome = species_name
          )
        })
    )
  } %*>% { GenomeInfoDb::seqinfo(seq) <- seqinfo_; seq }
}

# Generate, summarize and merge all derived data for one species. The only
# inputs are a genome and a GFF file of gene models.
load_species <- function(species_name, con){

  # Type descriptions
  #- SpeciesName :: character
  #- FolderPath :: character
  #- GenomeDB :: FaFile
  #- GenomeSeq :: DNAStringSet
  #- TranscriptomeDB :: FaFile
  #- TranscriptomeSeq :: DNAStringSet
  #- GenomeSeq :: DNAStringSet
  #- CDS :: DNAStringSet
  #- GenomeSeqinfo :: Seqinfo
  #- Length :: integer | numeric
  #- Count :: integer
  #- NString :: GRanges
  #- ORFRanges :: GRanges
  #- GeneModelsDB :: TxDb
  #- GeneModelList :: GeneModelList
  #- NStringSummary :: [Length]
  #- Density :: density
  #- AA :: A | C | D | E | ... | X  -- enum of amino acids 
  #- Phase :: 0 | 1 | 2

  #- DNASummary :: {
  #-   n_triple      = [Count],
  #-   initial_codon = [Count],
  #-   final_codon   = [Count]
  #- }
  #- AASummary :: {
  #-   initial_residue   = PairList (AA, Count),
  #-   final_residue     = PairList (AA, Count),
  #-   has_internal_stop = Bool,
  #-   comp              = Matrix Count,
  #-   table             = Table {
  #-     seqids = Seqid, 
  #-     length = Length
  #-   }
  #- }
  #- GFFSummary :: {
  #-   seqstats    = data.frame,
  #-   mRNA_length = NumericSummary,
  #-   CDS_length  = NumericSummary,
  #-   exon_length = NumericSummary
  #- }
  #- NumericSummary :: {
  #-   min     = [Num],
  #-   q25     = [Num],
  #-   median  = [Num],
  #-   q75     = [Num],
  #-   max     = [Num],
  #-   mean    = [Num],
  #-   sd      = [Num],
  #-   n       = Int,
  #-   density = Density
  #- }
  #- GRangesSummary :: {
  #-   seqstats = Table {
  #-     seqid = Seqid,
  #-     min   = Int,
  #-     max   = Int
  #-   },
  #-   width = NumericSummary
  #- }

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

  # genome and summary_genome
  #- _ :: SpeciesName -> FolderPath -> FilePath
  get_genome_filename(species_name, dir=con@input@fna_dir) %>%
      .cacher("genome_filename") %>>%
      #- _ :: FilePath -> Genome
      load_dna %>%
        .cacher("genomeDB") %>>%
        #- _ :: Genome -> DNASummary
        summarize_dna %>%
        .cacher("summary_genome") %>%

    # seqinfo
    .view("genomeDB") %>%
      m_db2seq(species_name) %>% .cacher('genomeSeq') %>>%
      GenomeInfoDb::seqinfo() %>% .cacher('seqinfo') %>%

    # nstring and summary_nstring
    .view("genomeSeq") %>>%
      #- _ :: GenomeSeq -> NString
      derive_nstring %>%
      .cacher("nstring") %>>%
      #- _ :: NString -> NStringSummary
      summarize_nstring %>%
      .cacher("summary_nstring") %>%

    .view("genomeDB") %__%

    # gff and summary_gff
    #- _ :: SpeciesName -> FolderPath -> FilePath
    get_gff_filename(species_name, dir=con@input@gff_dir) %>%
      .tag("gff_filename") %>%
      #- _ :: FilePath -> GenomeSeqinfo -> m GeneModels
      rmonad::funnel(
        file = .,
        seqinfo_ = .view(., "seqinfo")
      ) %*>% load_gene_models %>% .cacher("gffDB") %>>%
      #- _ :: GeneModels -> GFFSummary
      summarize_gff %>% .cacher("summary_gff") %>%

    # orfgff and summary_orfgff
    .view("genomeSeq") %>>%
      #- _ :: GenomeSeq -> ORFRanges 
      derive_orfgff %>>%
      { synder::as_gff(., id=GenomicRanges::mcols(.)$seqid) } %>%
      .cacher("orfgff") %>>%
      #- _ :: GRanges -> GRangesSummary
      summarize_granges %>%
      .cacher("summary_orfgff") %>%

    # orffaa and summary_orffaa
    .view("genomeDB") %>%
      #- Genome -> ORFRanges -> DNASeqs
      rmonad::funnel(
        dna = .,
        gff = .view(., "orfgff")
      ) %*>%
      extractWithComplements %>>%
      #- DNASeqs -> AASeqs
      {
        list(format_warnings=format_translation_warning)
        Biostrings::translate(., if.fuzzy.codon="solve")
      } %>%
      .cacher("orffaa") %>>%
      #- AASeqs -> AASummary
      summarize_faa %>% .cacher("summary_orffaa") %>%

    # transgff
    .view("gffDB") %>>%
      #- GeneModels ->
      #- GeneModelList {
      #-   cds_id    = Int,
      #-   cds_name  = String,
      #-   exon_rank = Int
      #- }
      GenomicFeatures::cdsBy(by="tx", use.names=TRUE) %>%
      .cacher("pre_transcript") %>>%
      #- GeneModelList -> [Phase]
      {

        "Find the phase of the first CDS in each model. This will be zero if the
        model is complete."


        # NOTE: in load_gene_models, I set the 'cds_name' to 'phase'.
        # exon_rank unfortunately is what is says it is: exon rank, not CDS rank.
        # So the code:
        # R> unlist(.)[mcols(unlist(.))$exon_rank == 1]$cds_name %>% as.integer
        # does not work. So instead I have to do the following, which is vastly slower:
        sapply(., function(x) GenomicRanges::mcols(x)$cds_name[1]) %>% as.integer

      } %>_%
      #- [Phase] -> *Warning
      {

        "Warn if there are any incomplete models"

        incomplete <- sum(. != 0)
        total <- length(.)

        if(incomplete > 0){
          warning(sprintf(
          "%s of %s gene models in %s are incomplete (the phase of the first CDS is
          not equal to 0). This does not include partial gene models that happen to
          begin in-phase. Details are recorded in the species_summary field
          model_phases.",
          incomplete, total, species_name))
        }

      } %>% .cacher("aa_model_phase") %>%
    .view("pre_transcript") %>>%
      #- GeneModelList -> GeneModelList  -- trim starts
      {

        "Offset start based on phase.
        
        WARNING: Abominable hack. To get around the GenomicFeatures issue
        reported [here](https://support.bioconductor.org/p/101245/), I encode
        the phase information in the CDS name. Here I extract that data and
        offset the reading frame as needed."

        ori <- .

        unlist(ori) %>%
          {
            meta <- GenomicRanges::mcols(.)
            meta$phase <- as.integer(meta$cds_name)
            GenomicRanges::start(.) <- ifelse(
              meta$exon_rank == 1,
              GenomicRanges::start(.) + meta$phase,
              GenomicRanges::start(.)
            )
            .
          } %>%
          relist(ori)

      } %>%
      .cacher("transgffDB") %>%

    # aa and summary_aa
    .view("genomeDB") %>%
      #- Genome -> GeneModelList -> CDS
      rmonad::funnel(
        x = .,
        transcripts = .view(., "transgffDB")
      ) %*>%
      GenomicFeatures::extractTranscriptSeqs %>>%
      #- CDS -> AASeqs
      {
        list(format_warnings=format_translation_warning)
        Biostrings::translate(., if.fuzzy.codon="solve")
      } %>%
      .cacher("faa") %>>% 
      #- AASeqs -> AASummary
      summarize_faa %>% .cacher("summary_aa") %>_%
      #- AASummary -> *Warning
      check_for_internal_stops %>%

    # aa_model_phase
    .view("faa") %>%
      #- [Phase] -> AASeqs -> {
      #-   summary = NamedList (0 = Count, 1 = Count, 2 = Count),
      #-   incomplete_models = [Seqid]
      #- }
      rmonad::funnel(phases = .view(., "aa_model_phase"), aa = .) %*>% {
        # This should always be true. If it is not, there is a logical problem
        # in the code, not the data.
        stopifnot(length(aa) == length(phases))
        list(
          summary = factor(phases) %>% table,
          incomplete_models = names(aa)[phases != 0]
        )
      } %>% .cacher("aa_model_phase") %>%

    # transfna
    .view("gffDB") %>>%
      #- GeneModels -> GenomicRanges
      GenomicFeatures::exonsBy(by="tx", use.names=TRUE) %>%
      #- Genome -> GenomicRanges -> Transcriptome
      rmonad::funnel(
        x = .view(., "genomeDB"),
        transcripts = .
      ) %*>%
      get_trans_dna %>% .cacher("transcriptomeDB") %>>%
      #- Transcriptome -> DNASummary
      summarize_dna %>% .cacher("summary_transfna") %>%

    # transorfgff
    .view("transcriptomeDB") %>>%
      #- Transcriptome -> TranscriptomeSeq
      scanFa_trw %>>%
      #- TranscriptomeSeq -> ORFRange
      derive_orfgff %>>%
      { synder::as_gff(., id=GenomicRanges::mcols(.)$seqid) } %>%
      .cacher("transorfgff") %>>%
      #- ORFRange -> GRangesSummary
      summarize_granges %>%
      .cacher("summary_transorfgff") %>%
      #- AASummary -> GRangesSummary -> *Warning
      rmonad::funnel(
        aa_summary = .view(., "summary_aa"),
        trans_summary = .,
      ) %*>%
      check_protein_transcript_match %>%

    .view("transcriptomeDB") %>%
      m_db2seq(species_name) %>% .cacher("transcriptomeSeq") %>%

    # transorfaa and summary_transorfaa
    .view("transcriptomeDB") %>%
      #- Transcriptome -> ORFRange -> CDS
      rmonad::funnel(
        dna = .,
        gff = .view(., "transorfgff")
      ) %*>%
      extractWithComplements %>>%
      #- CDS -> AASeqs
      {
        list(format_warnings=format_translation_warning)
        Biostrings::translate(., if.fuzzy.codon="solve")
      } %>% .cacher("transorfaa") %>>%
      #- AASeqs -> AASummary
      summarize_faa %>% .cacher("summary_transorfaa")

}


load_synmap_meta <- function(tseqinfo, fseqinfo, target_species, focal_species, syndir){

  "Load the synteny map for the given focal and target species. Then cache the
  synteny map and return the cached filename and a summary of the map as a
  `synteny_meta` object."

  if(length(target_species) != 1  ||  length(focal_species) != 1){
    stop("More than one species found in Seqinfo file, this should not happen.")
  }

  if(any(is.na(target_species)) || any(is.na(focal_species))){
    stop("Species must be set in each Seqinfo object (must not be NA)")
  }

  get_synmap_filename(
    focal_name  = focal_species,
    target_name = target_species,
    dir         = syndir
  ) %>% cacher(c("synmap_file", target_species)) %>>%
    synder::read_synmap(
      seqinfo_a = fseqinfo,
      seqinfo_b = tseqinfo
    ) %>% cacher(c("synmap", target_species)) %>>%
    summarize_syn %>% cacher(c("synmap_summary", target_species))
}


#' Load primary data
#'
#' @export
primary_data <- function(con){

  "
  Derive secondary data from the minimal required inputs
  
  The required inputs to Fagin are genome fasta files, GFFs, synteny maps, a
  species tree, and the focal species name. This function derives all required
  secondary data, but does no analysis specific to the query intervals (e.g.
  the orphan genes candidates).
  
  This function does not handle the full data validation. All individual
  pieces are validated, but not the relationships between them. This is the
  role of the `validate_derived_inputs` function.
  "

  load_tree(con@input@tree) %>% cacher("tree") %>>%
    {

      "Extract species list from the species tree. The tree, rather than the
      input files, determines which species are used in the analysis."

      .$tip.label

    } %>_% {
      "Assert that the focal species is in the phylogenetic tree"

      if(! (con@input@focal_species %in% .)){
        msg <- "Focal species '%s' is not in the species tree [%s]"
        stop(sprintf(msg, ., paste(targets, collapse=", ")))
      }
    } %>>% {

      "Remove spaces in the species names."

      gsub(" ", "_", .)

    } %>% cacher("all_species") %>>% {
      setdiff(., con@input@focal_species)
    } %>%
      cacher('target_species') %>%
      view('all_species') %>%

    rmonad::loop(
      FUN = load_species,
      con = con
    ) %__%

    con@input@query_gene_list %>>%
      load_gene_list %>% cacher("query_genes") %__%

    con@input@control_gene_list %>>%
      load_gene_list %>% cacher("control_genes") %>%
      
    rmonad::view("target_species") -> m

    # FIXME: this is too ugly, maybe can add some better handling in rmonad?
    # Loop doesn't quite work here, because I need two things from inside the
    # Rmonad.
    # I want something simple like this:
    #   rmonad::foo(
    #     load_synmap_meta,
    #     tseqinfo = view(c("seqinfo", target)) 
    #     fseqinfo = view(c("seqinfo", focal))
    #     fseqinfo = con@input@syn_dir
    #   )
    if(get_OK(m, m@head)){
      focal <- con@input@focal_species
      for(target in rmonad::get_value(m, tag='target_species')[[1]]){
        m <-
          rmonad::funnel(
            tseqinfo = rmonad::view(m, c("seqinfo", target)),
            fseqinfo = rmonad::view(m, c("seqinfo", focal))
          ) %*>%
          load_synmap_meta(
            target_species = target,
            focal_species = focal,
            syndir=con@input@syn_dir
          )
      }
    }

    m

}
