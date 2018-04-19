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
  #-   table       = data.frame,
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
  #-   table = Table {
  #-     seqid = Seqid,
  #-     min   = Int,
  #-     max   = Int
  #-   },
  #-   width = NumericSummary
  #- }

  .view <- function(m, tag){
    rmonad::view(m, tag, species_name)
  }

  .tag <- function(m, tag){
    rmonad::tag(m, tag, species_name)
  }

  # genome and summary_genome
  #- _ :: SpeciesName -> FolderPath -> FilePath
  get_readable_filename(
    species_name,
    dir = con@input@fna_dir, 
    ext = c("fna", "fa", "fasta") # allowed extensions genome FASTA file
  )  %>>%
    #- _ :: FilePath -> Genome
    load_dna %>%
      .tag("genomeDB") %>>%
      #- _ :: Genome -> DNASummary
      summarize_dna %>%
      .tag("summary_genome") %>%

    # seqinfo
    .view("genomeDB") %>>%
      convert_FaFile_to_XStringSet %>% .tag('genomeSeq') %>>%
      GenomeInfoDb::seqinfo() %>% .tag('seqinfo') %>%

    # nstring and summary_nstring
    .view("genomeSeq") %>>%
      #- _ :: GenomeSeq -> NString
      derive_nstring %>%
      .tag("nstring") %>>%
      #- _ :: NString -> NStringSummary
      summarize_nstring %>%
      .tag("summary_nstring") %>%

    .view("genomeDB") %__%

    # gff and summary_gff
    #- _ :: SpeciesName -> FolderPath -> FilePath
    get_readable_filename(
      species_name,
      dir = con@input@gff_dir,
      ext = c("gff3","gff")
    ) %>%
      #- _ :: FilePath -> GenomeSeqinfo -> m GeneModels
      rmonad::funnel(
        file = .,
        seqinfo_ = .view(., "seqinfo")
      ) %*>% load_gene_models %>% .tag("gffDB") %>>%
      #- _ :: GeneModels -> GFFSummary
      summarize_gff %>% .tag("summary_gff") %>%

    # orfgff and summary_orfgff
    .view("genomeSeq") %>>%
      #- _ :: GenomeSeq -> ORFRanges 
      derive_orfgff %>>%
      convert_GRanges_to_SynderGFF %>%
      .tag("orfgff") %>>%
      #- _ :: GRanges -> GRangesSummary
      summarize_granges %>%
      .tag("summary_orfgff") %>%

    #- Genome -> ORFRanges -> DNASeqs
    {rmonad::funnel(
      dna = .view(., "genomeDB"),
      gff = .view(., "orfgff")
    )} %*>%
      extractWithComplements %>>%
      translate(label=species_name) %>%
      .tag("orffaa") %>>%
      #- AASeqs -> AASummary
      summarize_faa %>% .tag("summary_orffaa") %>%

    {rmonad::funnel(
      gffDB = .view(., 'gffDB'),
      genomeDB = .view(., "genomeDB")
    )} %*>% m_get_proteins(species_name=species_name) %>%

    # aa_model_phase
    .view("faa") %>%
      #- [Phase] -> AASeqs -> {
      #-   summary = NamedList (0 = Count, 1 = Count, 2 = Count),
      #-   incomplete_models = [Seqid]
      #- }
      rmonad::funnel(phases = .view(., "aa_model_phase"), aa = .) %*>%
      summarize_phase %>% .tag("summary_phase") %>%

    # transfna
    .view("gffDB") %>>%
      #- GeneModels -> GenomicRanges
      GenomicFeatures::exonsBy(by="tx", use.names=TRUE) %>%
      #- Genome -> GenomicRanges -> Transcriptome
      rmonad::funnel(
        x = .view(., "genomeDB"),
        transcripts = .
      ) %*>%
      get_trans_dna %>% .tag("transcriptomeDB") %>>%
      #- Transcriptome -> DNASummary
      summarize_dna %>% .tag("summary_transfna") %>%

    # transorfgff
    .view("transcriptomeDB") %>>%
      #- Transcriptome -> TranscriptomeSeq
      convert_FaFile_to_XStringSet %>>%
      #- TranscriptomeSeq -> ORFRange
      derive_orfgff %>>%
      #- ORFRange -> ORFRange  -- slightly different format, some checking
      convert_GRanges_to_SynderGFF %>% .tag("transorfgff") %>>%
      #- ORFRange -> GRangesSummary
      summarize_granges %>% .tag("summary_transorfgff") %>%
      #- AASummary -> GRangesSummary -> *Warning
      {rmonad::funnel(
        aa_summary = .view(., "summary_aa"),
        trans_summary = .view(., "summary_transorfgff")
      )} %*>%
      check_protein_transcript_match %>%

    .view("transcriptomeDB") %>>%
      convert_FaFile_to_XStringSet %>% .tag("transcriptomeSeq") %>%

    #- Transcriptome -> ORFRange -> CDS
    {rmonad::funnel(
      dna = .view(., "transcriptomeDB"),
      gff = .view(., "transorfgff")
    )} %*>%
      extractWithComplements %>>%
      #- CDS -> AASeqs
      translate(label=species_name) %>% .tag("transorfaa") %>>%
      #- AASeqs -> AASummary
      summarize_faa %>% .tag("summary_transorfaa")

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
  ) %>% rmonad::tag(c("synmap_file", target_species)) %>>%
    synder::read_synmap(
      seqinfo_a = fseqinfo,
      seqinfo_b = tseqinfo
    ) %>% rmonad::tag(c("synmap", target_species)) %>>%
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
