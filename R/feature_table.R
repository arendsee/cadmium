#' Build table of binary features
buildFeatureTable <- function(result, query, config){
  orphans <- query$orphans 

  # Perform bonferoni corrections on all pvalue cutoffs
  p2p_cutoff <- config$prot2prot_pval     / length(orphans)
  p2a_cutoff <- config$prot2allorf_pval   / length(orphans)
  d2d_cutoff <- config$dna2dna_pval       / length(orphans)
  p2t_cutoff <- config$prot2transorf_pval / length(orphans)

  getGaps <- function(q2g){
    if(is.null(q2g)){
      q2g
    } else {
      q2g$query 
    }
  }

  # Synteny is scrambled
  scr <- orphans %in% result$scrambled
  # at least one search interval overlaps a target CDS
  cds <- orphans %in% (result$features$CDS$query %>% unique)
  # at least one search interval overlaps a target mRNA
  rna <- orphans %in% (result$features$mRNA$query %>% unique)
  # at least search interval overlaps a N-string
  nst <- orphans %in% getGaps(result$query2gap)
  # number of confirmed indels (based on search interval size)
  ind <- orphans %in% result$ind.stats$indeled.queries
  # number of confirmed resized (based on search interval size)
  res <- orphans %in% result$ind.stats$resized.queries
  # the query has an ortholog in the target
  gen <- orphans %in% (result$prot2prot$map %>% subset(pval < p2p_cutoff) %$% query)
  # ORF match in SI
  orf <- orphans %in% (result$prot2allorf$map %>% subset(pval < p2a_cutoff) %$% query)
  # has nucleotide match in SI
  nuc <- orphans %in% (result$dna2dna$map %>% subset(pval < d2d_cutoff) %$% query)
  # ORF match to spliced transcript (possibly multi-exonic)
  trn <- orphans %in% (result$prot2transorf$map %>% subset(pval < p2t_cutoff) %$% query)
  # at least one search interval maps off scaffold (flag==4 or flag==5)
  una <- orphans %in% result$unassembled
  # search interval was not processed for technical reasons (e.g. too big)
  tec <- orphans %in% result$skipped

  data.frame(
    seqid = orphans,
    scr   = scr,
    cds   = cds,
    rna   = rna,
    gen   = gen,
    nst   = nst,
    ind   = ind,
    res   = res,
    orf   = orf,
    nuc   = nuc,
    trn   = trn,
    una   = una,
    tec   = tec
  )
}

buildLabelsTree <- function(feats, config){

  require(yaml)
  root <- yaml.load_file(config$f_decision_tree) %>%
    as.Node(replaceUnderscores=FALSE)

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
  root$Get(toTable, filterFun = isLeaf, simplify=FALSE) %>%
    # Remove any NULL elements (so rbind doesn't crash)
    lapply(function(x) if(!is.null(x)) { x })           %>%
    # Bind all tables into one
    do.call(what=rbind)                                 %>%
    # Remove missing elements ... TODO: is this necessary?
    filter(!is.na(seqid))                               %>%
    # Remove rownames
    set_rownames(NULL)
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

#' Merge labels for all species
determineLabels <- function(query, results, config){

  descriptions <- c(
    O1 = 'genic: known gene',
    O2 = 'genic: unknown ORF on known mRNA',
    O3 = 'genic: unknown ORF off known mRNA',
    U1 = 'unknown: maps off the scaffold',
    U2 = 'unknown: possible indel',
    U3 = 'unknown: possible N-string',
    U4 = 'unknown: possible resized',
    U5 = 'unknown: syntenically scrambled',
    U6 = 'unknown: seriously unknown',
    U7 = 'unknown: skipped for technical reasons',
    N1 = 'non-genic: no gene in SI',
    N2 = 'non-genic: CDS in SI',
    N3 = 'non-genic: mRNA but not CDS in SI'
  )

  features <- lapply(results, buildFeatureTable, query, config)

  labelTrees <- lapply(features, buildLabelsTree, config)

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

#' Given a species tree and a set of labels, determine orphan origin
determineOrigins <- function(labels, config){

  require(ape)
  require(data.tree)

  root <- read.tree(config$f_tree) %>% as.Node(replaceUnderscores=FALSE)
  classify <- function(node){
      if(node$isLeaf){
        if(node$name == config$focal_species){
          node$cls <- NULL
        } else {
          node_labels <- labels$labels[[node$name]] %>% dplyr::arrange(seqid)
          node$cls <- node_labels$primary
          names(node$cls) <- node_labels$seqid
        }
      } else {
        child_cls <- lapply(node$children, classify)
        if(length(child_cls) != 2){
          warning('Species tree must be bifurcating')
        }
        a <- child_cls[[1]]
        b <- child_cls[[2]]
        if(!is.null(a) && !is.null(b)){
          stopifnot(names(node$cls) == names(a))
          stopifnot(names(node$cls) == names(b))
          node$cls <- ifelse(a == "Unknown" & b == "Non-ORFic", "Non-ORFic", "Unknown")
          node$cls <- ifelse(a == "ORFic" | b == "ORFic", "ORFic", node$cls)
        } else if(!is.null(a)){
          stopifnot(names(node$cls) == names(a))
          node$cls <- a
        } else if(!is.null(b)){
          stopifnot(names(node$cls) == names(b))
          node$cls <- b
        }
      }
      node$gen <- sum(node$cls == 'ORFic')
      node$non <- sum(node$cls == 'Non-ORFic')
      node$unk <- sum(node$cls == 'Unknown')
      invisible(node$cls)
  }

  findFocalSpecies <- function(node){
    if(node$name == config$focal_species){
      return(node)
    }
    if(!node$isLeaf){
      for(child in node$children){
        n <- findFocalSpecies(child)
        if(!is.null(n)){
          return(n)
        }
      }
    }
    return(NULL)
  }

  setAncestor <- function(node){
    # You should be at an leaf node only at the first level, and the leaf must
    # be the focal species
    if(node$isLeaf){
      if(node$name != config$focal_species){
        warning('You should start setAncestor from the focal species')
      }
    } else {
      for(child in node$children){
        if(child$visited == 0){
          stopifnot(names(node$cls) == names(classify(child)))
          node$cls <- child$cls
          node$gen <- sum(node$cls == 'gen')
          node$non <- sum(node$cls == 'non')
          node$unk <- sum(node$cls == 'unk')
        }
      }
    }
    node$visited <- node$visited + 1
    if(!node$isRoot){
      setAncestor(node$parent)
    }
  }
  fs <- findFocalSpecies(root)
  root$Set(visited=0)
  setAncestor(fs)

  d <- fs$Get('cls', traversal='ancestor')
  d[[1]] <- NULL
  backbone <- do.call(cbind, d) %>%
    as.data.frame %>%
    set_names(paste0('ps', 1:(root$height-1)))

  list(
    root=root,
    backbone=backbone
  )
}
