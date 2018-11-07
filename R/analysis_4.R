default_class_rule <- function(a, b){
  cls <- ifelse(              a == "ORFic"   | b == "ORFic"    , "ORFic"     , NA )
  cls <- ifelse(is.na(cls) & (a == "Unknown" | b == "Unknown") , "Unknown"   , cls)
  cls <- ifelse(is.na(cls)                                     , "Non-ORFic" , cls)
  cls
}

by_primary <- function(node_labels){
  node_labels$primary
}

by_secondary <- function(node_labels){
  node_labels$secondary
}

classify <- function(
  node,
  labels,
  get_classifier = by_primary,
  class_rule     = default_class_rule
){

    "For each node in the phylogenetic tree, count the number of members that
     are ORFic, non-ORFic, and unknown"

    if(node$isLeaf){
      if(node$name == con@input@focal_species){
        node$cls <- NULL
      } else {
        node_labels <- labels[[node$name]] %>% dplyr::arrange(.data$seqid)
        node$cls <- get_classifier(node_labels)
        names(node$cls) <- node_labels$seqid
      }
    } else {
      child_cls <- lapply(node$children, classify, labels, get_classifier, class_rule)
      if(length(child_cls) != 2){
        warning('Species tree must be bifurcating')
      }
      a <- child_cls[[1]]
      b <- child_cls[[2]]
      if(!is.null(a) && !is.null(b)){
        stopifnot(names(node$cls) == names(a))
        stopifnot(names(node$cls) == names(b))
        node$cls <- class_rule(a, b)
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
  if(node$name == con@input@focal_species){
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

setAncestor <- function(node, ...){
  # You should be at an leaf node only at the first level, and the leaf must
  # be the focal species
  if(node$isLeaf){
    if(node$name != con@input@focal_species){
      stop('You should start setAncestor from the focal species')
    }
  } else {
    for(child in node$children){
      if(child$visited == 0){
        stopifnot(names(node$cls) == names(classify(child, ...)))
        node$cls <- child$cls
        node$gen <- sum(child$cls == 'ORFic')
        node$non <- sum(child$cls == 'Non-ORFic')
        node$unk <- sum(child$cls == 'Unknown')
      }
    }
  }

  # the `visited` field should never be greater than 1
  node$visited <- node$visited + 1

  if(!node$isRoot){
    setAncestor(node$parent, ...)
  }
}


#' Given a species tree and a set of labels, determine orphan origin
#'
#' @param labels A list of data that includes a \code{labels} field
#' @param con A Config object
#' @param ... Arguments that will be passed to \code{classify}
#' @export
determine_origins <- function(m, con){

  message("Determining origins")

  m %>%
    rmonad::view('query_labels')   %>>%
      .determine_origins(con) %>% rmonad::tag('query_origins') %>%
    rmonad::view('control_labels') %>>%
      .determine_origins(con) %>% rmonad::tag('control_origins')
}
.determine_origins <- function(labels, con, ...){

  labels <- labels$labels

  if((sapply(labels, nrow) %>% unique %>% length) != 1)
    stop("All label tables must have equal numbers of rows")

  root <- ape::read.tree(con@input@tree) %>% data.tree::as.Node(replaceUnderscores=FALSE)

  fs <- findFocalSpecies(root)
  root$Set(visited=0)
  setAncestor(fs, labels=labels, ...)

  d <- fs$Get('cls', traversal='ancestor')
  d[[1]] <- NULL

  backbone <- do.call(cbind, d) %>%
    as.data.frame %>%
    set_names(paste0('ps', 1:(root$height-1))) %>% {
      l <- levels(.[[1]])
      l <- sub("Non-ORFic", "N", l)
      l <- sub("ORFic",     "O", l)
      l <- sub("Unknown",   "U", l)
      lapply(., function(s) {levels(s) <- l; s})
    } %>%
    as.data.frame

  classStr <- apply(backbone, 1, paste, collapse="")

  classSum <- classStr %>% factor %>% summary(maxsum=Inf) %>% {
    data.frame(
      group = names(.),
      count = .
    )} %>% { rownames(.) <- NULL; . }

  list(
    root     = root,
    backbone = backbone,
    classStr = classStr,
    classSum = classSum
  )

}
