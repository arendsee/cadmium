# Generate, summarize and merge all derived data for one species. The only
# inputs are a genome and a GFF file of gene models.
load_species <- function(species_name, con){

  message("Loading ", species_name)

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

  get_gff_from_txdb <- function(x, ...){
    synder::as_gff(
      x,
      id=GenomicRanges::mcols(x)$tx_name,
      ...
    )
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
      make_seqinfo(species_name) %>% .tag('seqinfo') %>%

    # nstring and summary_nstring
    .view("genomeSeq") %>>%
      #- _ :: GenomeSeq -> NString
      derive_nstring %>% .tag("nstring") %>>%
      #- _ :: NString -> NStringSummary
      summarize_nstring %>% .tag("summary_nstring") %>%

    .view("genomeDB") %>>%

    # gff and summary_gff
    #- _ :: SpeciesName -> FolderPath -> FilePath
    {get_readable_filename(
      species_name,
      dir = con@input@gff_dir,
      ext = c("gff3","gff")
    )} %>%
      #- _ :: FilePath -> GenomeSeqinfo -> m GeneModels
      rmonad::funnel(
        file = .,
        seqinfo_ = .view(., "seqinfo")
      ) %*>% load_gene_models %>% .tag("gffDB") %>>%
      #- _ :: GeneModels -> GFFSummary
      summarize_gff %>% .tag("summary_gff") %>%

    # GeneModels -> mRNA
    .view("gffDB") %>>%
      GenomicFeatures::transcripts() %>%
      rmonad::funnel(
        x = .,
        seqinfo_ = .view(., "seqinfo")
      ) %*>%
      get_gff_from_txdb(type = "mRNA") %>% .tag("mRNA") %>%

    # GeneModels -> CDS
    .view("gffDB") %>>%
      GenomicFeatures::cds() %>%
      rmonad::funnel(
        x = .,
        seqinfo_ = .view(., "seqinfo")
      ) %*>%
      get_gff_from_txdb(type = "CDS") %>% .tag("CDS") %>%

    # GeneModels -> Exon
    .view("gffDB") %>>%
      GenomicFeatures::exons() %>%
      rmonad::funnel(
        x = .,
        seqinfo_ = .view(., "seqinfo")
      ) %*>%
      get_gff_from_txdb(type = "exon") %>% .tag("exon") %>%

    # orfgff and summary_orfgff
    .view("genomeSeq") %>>%
      #- _ :: GenomeSeq -> ORFRanges
      derive_genomic_ORFs(con) %>>%
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
      extract_with_complements %>>%
      fuzzy_translate(label=species_name) %>% .tag("orffaa") %>>%
      #- AASeqs -> AASummary
      summarize_faa %>% .tag("summary_orffaa") %>%

    {rmonad::funnel(
      gffDB    = .view(., "gffDB"),
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
      get_trans_dna(species_name=species_name) %>% .tag("transcriptomeDB") %>>%
      #- Transcriptome -> DNASummary
      summarize_dna %>% .tag("summary_transfna") %>%

    # transorfgff
    .view("transcriptomeDB") %>>%
      #- Transcriptome -> TranscriptomeSeq
      convert_FaFile_to_XStringSet %>>%
      #- TranscriptomeSeq -> ORFRange
      derive_transcript_ORFs(con) %>>%
      #- ORFRange -> ORFRange  -- slightly different format, some checking
      convert_GRanges_to_SynderGFF %>% .tag("transorfgff") %>>%
      #- ORFRange -> GRangesSummary
      summarize_granges %>% .tag("summary_transorfgff") %>%
      #- AASummary -> GRangesSummary -> *Warning
      {rmonad::funnel(
        aa_summary    = .view(., "summary_aa"),
        trans_summary = .view(., "summary_transfna")
      )} %*>%
      check_protein_transcript_match %>%

    .view("transcriptomeDB") %>>%
      convert_FaFile_to_XStringSet %>% .tag("transcriptomeSeq") %>%

    #- Transcriptome -> ORFRange -> CDS
    {rmonad::funnel(
      dna = .view(., "transcriptomeDB"),
      gff = .view(., "transorfgff")
    )} %*>%
      extract_with_complements %>>%
      #- CDS -> AASeqs
      fuzzy_translate(label=species_name) %>% .tag("transorfaa") %>>%
      #- AASeqs -> AASummary
      summarize_faa %>% .tag("summary_transorfaa")
}

load_synmap_meta <- function(m, focal, target, syndir){
  m %>>%
    {

      "Load synteny maps"

      get_synmap_filename(
        focal_name  = focal,
        target_name = target,
        dir         = syndir
      )
    } %>% rmonad::tag("synmap_file", target) %>%
    funnel(
      seqinfo_a = rmonad::view(., "seqinfo", focal),
      seqinfo_b = rmonad::view(., "seqinfo", target) 
    ) %*>%
      synder::read_synmap %>% rmonad::tag("synmap", target) %>>%
      summarize_syn %>% rmonad::tag("synmap_summary", target)
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

  {
    as_monad(get_species(con)) %>%
    rmonad::loop(
      FUN = load_species,
      con = con
    ) %>>%

    {
      load_gene_list(con@input@query_gene_list)
    } %>% rmonad::tag("query_genes") %>>%

    {
      load_gene_list(con@input@control_gene_list)
    } %>% rmonad::tag("control_genes")

  } -> m

  for(target in get_targets(con)){
    m <- load_synmap_meta(
      m,
      focal  = con@input@focal_species,
      target = target,
      syndir = con@input@syn_dir
    )
  }

  m
}
