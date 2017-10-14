#' Build table of binary features
buildFeatureTable <- function(d, con){

  n <- length(d$queries)

  { 

    "Set p-value cutoffs for each set of alignments (aa-vs-aa, aa-vs-orf,
    aa-vs-mRNA, gene-vs-genome). The p-value threshold set in the configuration
    is the overall desired p-value (defualting to 0.05). Here we need to
    account for multiple testing. For now we use a Bonferoni correction."

    list(
      p2p_cutoff = con@alignment@thresholds@prot2prot     / length(d$queries),
      p2a_cutoff = con@alignment@thresholds@prot2allorf   / length(d$queries),
      d2d_cutoff = con@alignment@thresholds@dna2dna       / length(d$queries),
      p2t_cutoff = con@alignment@thresholds@prot2transorf / length(d$queries)
    )

  } %*>% {

    "From the raw data, distill a binary feature table"

    met <- GenomicRanges::mcols(d$gene2genome$map) %>%
           as.data.frame %>%
           subset(pval < d2d_cutoff)

    rmonad::funnel(
      seqid = d$queries,
      # matches somewhere in at least one search interval
      nuc = d$queries %in% met$query,
      # matches CDS in at least one search interval
      cds = d$queries %in% subset(met, cds_match)$query,
      # matches transcript in at least one search interval
      rna = d$queries %in% subset(met, mrna_match)$query,
      # matches exon in at least one search interval
      exo = d$queries %in% subset(met, exon_match)$query,
      # at least search interval overlaps a N-string
      nst = d$queries %in% as.character(d$gapped$query),
      # number of confirmed indels (based on search interval size)
      ind = d$queries %in% d$indels,
      # the query has an ortholog in the target
      gen = d$queries %in% subset(d$aa2aa$map, pval < p2p_cutoff)$query ,
      # ORF match in SI
      orf = d$queries %in% subset(d$aa2orf$map, pval < p2a_cutoff)$query ,
      # ORF match to spliced transcript (possibly multi-exonic)
      trn = d$queries %in% subset(d$aa2transorf$map, pval < p2t_cutoff)$query ,
      # synteny is scrambled
      scr = d$queries %in% subset(d$synder_summary, incoherent)$attr,
      # at least one search interval maps off scaffold
      una = d$queries %in% subset(d$synder_summary, unassembled)$attr,
      # search interval was not processed for technical reasons (e.g. too big)
      tec = d$queries %in% d$gene2genome$skipped
    )

  } %>_% {

    "Assert all feature vectors are of the correct length"

    stopifnot(sapply(., length) == n)

  } %>>% as.data.frame(stringsAsFactors=FALSE)
}

#' Load tertiary data
#'
#' @export
tertiary_data <- function(secondary_data, con){

  "For each species, build a feature table"

  buildFeatureTables <- function(d){
    ss <- lapply(d, buildFeatureTable, con=con)
    names(ss) <- names(secondary_data)
    rmonad::combine(ss)
  }

  qss <- buildFeaturesTables(secondary_data$query_results)
  css <- buildFeaturesTables(secondary_data$control_results)
  
  rmonad::funnel(query=qss, control=css)
  
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
determine_labels <- function(features, con){

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

  labelTrees_ <- lapply(features, function(feat) feat %>>% buildLabelsTree(con) ) %>% rmonad::combine()
  
  labels_ <- labelTrees_ %>>%
    {
      lapply(
        seq_along(.),
        function(i) .[[i]] %>>% labelTreeToTable(features[[i]])
      ) %>% magrittr::set_names(names(features))
    } %>>% rmonad::combine()

  label.summary_ <- labels_ %>>%
  {

    "Summarize labels as a table with the columns:

      1. secondary
      2. species
      3. count
      4. description
    "

    lapply(., dplyr::count, primary, secondary)             %>%
    reshape2::melt(id.vars=c('primary', 'secondary'))       %>%
    dplyr::rename(species=L1, count=value)                  %>%
    dplyr::select(secondary, species, count)                %>%
    tidyr::complete(secondary, species, fill=list(count=0)) %>%
    dplyr::mutate(count = as.integer(count))                %>%
    dplyr::mutate(description = descriptions[secondary])    %>%
    dplyr::arrange(description, secondary, species, count)  %>%
    as.data.frame
  }

  rmonad::funnel(
    trees   = labelTrees_,
    labels  = labels_,
    summary = label.summary_
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
