#' Given a species tree and a set of labels, determine orphan origin
#' @export
determine_origins <- function(labels, con){

  labels <- labels$labels

  if((sapply(labels, nrow) %>% unique %>% length) != 1)
    stop("All label tables must have equal numbers of rows")

  root <- ape::read.tree(con@input@tree) %>% data.tree::as.Node(replaceUnderscores=FALSE)

  apply_class_rule <- function(a, b){
    cls <- ifelse(              a == "ORFic"   | b == "ORFic"    , "ORFic"     , NA )
    cls <- ifelse(is.na(cls) & (a == "Unknown" | b == "Unknown") , "Unknown"   , cls)
    cls <- ifelse(is.na(cls)                                     , "Non-ORFic" , cls)
    cls
    # cls <- ifelse(a == "Unknown" & b == "Non-ORFic", "Non-ORFic", "Unknown")
    # cls <- ifelse(a == "ORFic" | b == "ORFic", "ORFic", cls)
  }

  classify <- function(node){

      "For each node in the phylogenetic tree, count the number of members that
       are ORFic, non-ORFic, and unknown"

      if(node$isLeaf){
        if(node$name == con@input@focal_species){
          node$cls <- NULL
        } else {
          node_labels <- labels[[node$name]] %>% dplyr::arrange(.data$seqid)
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
          node$cls <- apply_class_rule(a, b)
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

  setAncestor <- function(node){
    # You should be at an leaf node only at the first level, and the leaf must
    # be the focal species
    if(node$isLeaf){
      if(node$name != con@input@focal_species){
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
