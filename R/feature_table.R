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
      # Synteny is scrambled
      scr = d$queries %in% d$scrambled,
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
      # at least one search interval maps off scaffold
      una = d$queries %in% d$unassembled,
      # search interval was not processed for technical reasons (e.g. too big)
      tec = d$queries %in% d$gene2genome$skipped
    )

  } %>_% {

    "Assert all feature vectors are of the correct length"

    stopifnot(sapply(., length) == n)

  } %>>% as.data.frame(stringsAsFactors=FALSE)
}

tertiary_data <- function(secondary_data, config){

  "For each species, build a feature table"

  ss <- lapply(secondary_data, buildFeatureTable, config)
  names(ss) <- names(secondary_data)

  rmonad::combine(ss)
  
}



#' Merge labels for all species
determineLabels <- function(features, con){



  buildLabelsTree <- function(feats, con){
    root <- con@decision_tree %>%
      data.tree::as.Node(replaceUnderscores=FALSE)
    classify <- function(node, membership=NULL){
      if(is.null(membership)){
        membership <- rep(TRUE, nrow(feats))
      }
      node$membership <- membership
      node$N <- sum(membership)
      if(node$name %in% names(feats)){
        if(length(node$children) == 2){
          yes <-  feats[[node$name]] & membership
          no  <- !feats[[node$name]] & membership
          classify(node$children[[1]], yes)
          classify(node$children[[2]], no)
        }
      } else if(node$isRoot){
        classify(node$children[[1]], membership)
      }
    }
    classify(root)
    root
  }



  labelTreeToTable <- function(root, feats){
    toTable <- function(node) {
      # Build a dataframe for each node
      # For example:
      #         seqid primary secondary
      # 1 AT1G29418.1   ORFic        O1
      # 2 AT3G15909.1   ORFic        O1
      if(node$N > 0){
        d <- data.frame(seqid = feats$seqid[node$membership])
        d$primary <- node$primary
        d$secondary <- node$secondary
        d
      # If there are no members in the class
      } else {
        NULL
      }
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
      # Remove any NULL elements (so rbind doesn't crash)
      lapply(function(x) if(!is.null(x)) { x })           %>%
      # Bind all tables into one
      do.call(what=rbind)                                 %>%
      # Remove missing elements ... TODO: is this necessary?
      filter(!is.na(seqid))                               %>%
      # Remove rownames
      set_rownames(NULL)
  }



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
    N4 = 'non-genic: No DNA match'
  )

  labelTrees <- lapply(features, buildLabelsTree, con)

  labels <- lapply(
    names(labelTrees),
    function(x) labelTreeToTable(labelTrees[[x]], features[[x]])
  ) %>% set_names(names(labelTrees))


  lapply(names(labelTrees), function(x) lapply(features[[x]][2:11], sum) %>% unlist) %>%
    set_names(names(labelTrees))

  label.summary <- labels                                   %>% 
    lapply(count, primary, secondary)                       %>% 
    melt(id.vars=c('primary', 'secondary'))                 %>% 
    dplyr::rename(species=L1, count=value)                  %>% 
    dplyr::select(secondary, species, count)                %>% 
    tidyr::complete(secondary, species, fill=list(count=0)) %>% 
    dplyr::mutate(count = as.integer(count))                %>% 
    dplyr::mutate(description = descriptions[secondary])    %>% 
    dplyr::arrange(description, secondary, species, count)  %>% 
    as.data.frame

  list(
    trees   = labelTrees,
    labels  = labels,
    summary = label.summary
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

  SetEdgeStyle(root, label=GetEdgeLabel)
  SetNodeStyle(root, label=GetNodeLabel, shape=GetNodeShape)
  SetGraphStyle(root, rankdir="LR")

  plot(root)
}
