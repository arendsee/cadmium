# Build table of binary features
buildFeatureTable <- function(m, species_name, con, group, gene_tag){

  .view_target <- function(m, tag){
    rmonad::view(m, tag, species_name)
  }

  .view_focal <- function(m, tag){
    rmonad::view(m, tag, con@input@focal_species)
  }

  .view <- function(m, tag){
    rmonad::view(m, tag, group, species_name)
  }


  rmonad::funnel(
    genes = rmonad::view(m, gene_tag),
    gene2genome = rmonad::view(m, "gene2genome.pval", group),
    gaps = .view(m, "gaps"),
    indels = .view(m, "indels"),
    aa2aa = rmonad::view(m, "aa2aa.pval", group),
    aa2orf = rmonad::view(m, "aa2orf.pval", group),
    aa2transorf = rmonad::view(m, "aa2transorf.pval", group),
    synder_summary = .view(m, "summary_synder_flags"),
    synder_out = .view(m, "synder_out"),
    aln = .view(m, "gene2genome_aln"),
    genome_hits = .view(m, "gene2genome")
  ) %*>% {

    "Set p-value cutoffs for each set of alignments (aa-vs-aa, aa-vs-orf,
    aa-vs-mRNA, gene-vs-genome). The p-value threshold set in the configuration
    is the overall desired p-value (defaulting to 0.05)."

    p2p_cutoff <- con@alignment@thresholds@prot2prot
    p2a_cutoff <- con@alignment@thresholds@prot2allorf
    d2d_cutoff <- con@alignment@thresholds@dna2dna
    p2t_cutoff <- con@alignment@thresholds@prot2transorf

    met <- merge(
      subset(gene2genome, species == species_name),
      genome_hits$map[, c("query", "cds_match", "exon_match", "mrna_match")],
      by="query"
    ) %>%
      dplyr::select(-target, -species, -pval) %>%
      dplyr::filter(pval.adj < d2d_cutoff)

    data.frame(
      seqid = genes,
      # matches somewhere in at least one search interval
      nuc = genes %in% met$query,
      # matches CDS in at least one search interval
      cds = genes %in% subset(met, cds_match)$query,
      # matches transcript in at least one search interval
      rna = genes %in% subset(met, mrna_match)$query,
      # matches exon in at least one search interval
      exo = genes %in% subset(met, exon_match)$query,
      # at least search interval overlaps a N-string
      nst = genes %in% gaps$query,
      # number of confirmed indels (based on search interval size)
      ind = genes %in% indels,
      # the query has an ortholog in the target
      gen = genes %in% subset(aa2aa, pval.adj < p2p_cutoff & species == species_name)$query,
      # ORF match in SI
      orf = genes %in% subset(aa2orf, pval.adj < p2a_cutoff & species == species_name)$query,
      # ORF match to spliced transcript (possibly multi-exonic)
      trn = genes %in% subset(aa2transorf, pval.adj < p2t_cutoff & species == species_name)$query,
      # synteny is scrambled
      scr = genes %in% subset(synder_summary, incoherent)$attr,
      # at least one search interval maps off scaffold
      una = genes %in% subset(synder_summary, unassembled)$attr,
      # search interval was not processed for technical reasons (e.g. too big)
      tec = genes %in% aln$skipped,
      stringsAsFactors = FALSE
    )
  } %>% rmonad::tag("feature_table", group, species_name)
}

buildLabelsTree <- function(feats, con){

  "Merge all features into the feature decision tree. Each node contains a
  logical vector specifying which queries are members of the node. The names
  in the decision tree must match the names in the feature table."

  root <- con@decision_tree %>%
    data.tree::as.Node(replaceUnderscores=FALSE)
  classify <- function(node, membership=logical(0)){
    if(length(membership) == 0){
      membership <- rep(TRUE, nrow(feats))
    }
    node$membership <- membership
    node$N <- sum(membership)
    if(node$name %in% names(feats)){
      kids <- node$children
      if(length(kids) == 2){
        yes <-  feats[[node$name]] & membership
        no  <- !feats[[node$name]] & membership
        classify(kids[[1]], yes)
        classify(kids[[2]], no)
      } else if (length(kids) != 0) {
        stop("Nodes have either 0 or 2 children")
      }
    } else if(node$isRoot){
      classify(node$children[[1]], membership)
    }
  }
  classify(root)
  root
}

labelTreeToTable <- function(root, feats){

  "
  Build a dataframe for each node
  For example:
          seqid primary secondary
  1 AT1G29418.1   ORFic        O1
  2 AT3G15909.1   ORFic        O1

  Then bind the rows of the dataframes for each node.
  "

  toTable <- function(node) {
    if(!is.null(node$N) && node$N > 0){
      seqid     <- feats$seqid[node$membership]
      primary   <- node$primary
      secondary <- node$secondary
    # If there are no members in the class
    } else {
      seqid     <- character(0)
      primary   <- character(0)
      secondary <- character(0)
    }
    data.frame(
      seqid     = seqid,
      primary   = primary,
      secondary = secondary,
      stringsAsFactors=FALSE
    )
  }
  # Get a table for each class
  # --- NOTE: Without simplify=FALSE something truly dreadful happens.
  # The proper output (where a,b,c ... are dataframe columns):
  #   <Node_1> a b c
  #   <Node_2> d e f
  #   <Node_3> g h i
  #   <Node_4> j k l
  # The automagically borken output
  #   <Node_1> a
  #   <Node_2> b
  #   <Node_3> c
  #   <Node_4> d
  # The columns are recast as vectors, fed to the wrong children, and the
  # leftovers are tossed.
  root$Get(toTable, filterFun = data.tree::isLeaf, simplify=FALSE) %>%
    # Bind all tables into one
    do.call(what=rbind) %>%
    # Remove rownames
    set_rownames(NULL)
}


#' Merge labels for all species
#' @export
determine_labels <- function(m, con){

  message("Determining labels")

  # TODO: fix rmonad::views so that it returns a named list
  query_features <- list()
  control_features <- list()
  for(target in get_targets(con)){
    query_features[[target]]   <- rmonad::view(m, "feature_table", "query",   target)
    control_features[[target]] <- rmonad::view(m, "feature_table", "control", target)
  }

  rmonad::combine(query_features) %>>%
    .determine_labels(con) %>% rmonad::tag("query_labels") %__%
  rmonad::combine(control_features) %>>%
    .determine_labels(con) %>% rmonad::tag("control_labels")
}
.determine_labels <- function(features, con){

  descriptions <- c(
    O1 = 'genic: known gene',
    O2 = 'genic: unknown ORF on known mRNA',
    O3 = 'genic: unknown ORF off known mRNA',
    U1 = 'unknown: maps off the scaffold',
    U2 = 'unknown: possible indel',
    U3 = 'unknown: possible N-string',
    # U4 = 'unknown: possible resized',
    U5 = 'unknown: syntenically scrambled',
    U6 = 'unknown: seriously unknown',
    U7 = 'unknown: skipped for technical reasons',
    N1 = 'non-genic: DNA match to CDS',
    N2 = 'non-genic: DNA match to exon',
    N3 = 'non-genic: DNA match to mRNA',
    N4 = 'non-genic: No DNA match to known gene'
  )

  labelTrees_ <- lapply(features, buildLabelsTree, con)

  labels_ <- lapply(
      seq_along(labelTrees_),
      function(i) labelTrees_[[i]] %>% labelTreeToTable(features[[i]])
    ) %>% magrittr::set_names(names(features))

  # Summarize labels as a table with the columns:
  #
  #  1. secondary
  #  2. species
  #  3. count
  #  4. description
  label_summary_ <-
    lapply(labels_, dplyr::count, primary, secondary)       %>%
    reshape2::melt(id.vars=c('primary', 'secondary'))       %>%
    dplyr::rename(species=L1, count=value)                  %>%
    dplyr::select(secondary, species, count)                %>%
    tidyr::complete(secondary, species, fill=list(count=0)) %>%
    dplyr::mutate(count = as.integer(count))                %>%
    dplyr::mutate(description = descriptions[secondary])    %>%
    dplyr::arrange(description, secondary, species, count)  %>%
    as.data.frame

  list(
    trees   = labelTrees_,
    labels  = labels_,
    summary = label_summary_
  )
}


plotDecisionTree <- function(root){
  GetNodeLabel <- function(node) { sprintf('%s\n%s', node$name, node$N) }

  GetEdgeLabel <- function(node) {
    if(!node$isRoot){
      if(node$name == node$parent$children[[1]]$name){
        'yes'
      } else {
        'no'
      }
    } else {
      NULL
    }
  }

  GetNodeShape <- function(node) {
    if(node$isLeaf){
      'circle'
    } else {
      'square'
    }
  }

  data.tree::SetEdgeStyle(root, label=GetEdgeLabel)
  data.tree::SetNodeStyle(root, label=GetNodeLabel, shape=GetNodeShape)
  data.tree::SetGraphStyle(root, rankdir="LR")

  plot(root)
}


#' Load tertiary data
#'
#' For each species, build a feature table
#'
#' @export
tertiary_data <- function(m, con){

  message("Building feature tables")

  for(s in get_targets(con)){
    m <- buildFeatureTable(m, con, species_name=s, group = "query",   gene_tag = "query_genes")
    m <- buildFeatureTable(m, con, species_name=s, group = "control", gene_tag = "control_genes")
  }

  m
}
