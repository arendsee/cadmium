# A local (for now) utility function
.extract <- function(m, tags, simplify_names=TRUE){
  xs <- rmonad::get_value(m, tag=tags)
  for(tag in tags){
    names(xs) <- sub(pattern=tag, replacement="", names(xs))
  }
  if(simplify_names){
    names(xs) <- sub(pattern="//", replacement="/", names(xs))
    names(xs) <- sub(pattern="^/", replacement="", names(xs))
  }
  xs
}

# extract a tag from a list of things
.select_tag <- function(xs, pattern){
  xs <- xs[grepl(pattern=pattern, x=names(xs), perl=TRUE)]
  names(xs) <- sub(pattern=pattern, replacement="", x=names(xs))
  xs
}

# Create a table from a named list of named vectors
#
# @param xs named list of named vectors
# @param fill ANY - value to replace missing cells with
# @param sort logical - sort the columns
# @param remove_empty should columns filled with only this value be removed?
.rbindlist <- function(xs, fill=0, sort=TRUE, remove_empty=FALSE){
  chars <- lapply(xs, names) %>% unlist %>% unique
  if(sort)
    chars <- sort(chars)
  tbl <- do.call(rbind, lapply(xs, function(x) x[match(chars, names(x))]))
  colnames(tbl) <- chars
  tbl[is.na(tbl)] <- fill
  if(remove_empty){
    tbl <- tbl[, apply(tbl, 2, function(x) !all(x == fill))]
  }
  tbl
}

# Calculate GC content from a named vector of counts 
.get_gc <- function(x){
  (x['G'] + x['C']) / (x['G'] + x['C'] + x['A'] + x['T'])
}

# Summarize a list numeric vectors
make_numeric_summary_table <- function(nd){
    data.frame(
        min    = sapply(nd, function(x) x@min)
      , q25    = sapply(nd, function(x) x@q25)
      , median = sapply(nd, function(x) x@median)
      , q75    = sapply(nd, function(x) x@q75)
      , max    = sapply(nd, function(x) x@max)
      , mean   = sapply(nd, function(x) x@mean)
      , sd     = sapply(nd, function(x) x@sd)
      , n      = sapply(nd, function(x) x@n)
    )
}

# Given a data.frame with species as rownames: 1) create a new 'species'
# column, 2) remove the rownames, 3) sort by the vector 'levels' (this will
# normally be the phylogenetic order).
rownames_as_species <- function(x, levels=unique(rownames(x))){
  x$species <- rownames(x) %>% factor(levels=levels)
  rownames(x) <- NULL
  dplyr::arrange(x, species) %>%
    { .[c(ncol(.), 1:(ncol(.)-1))] } # put species first
}


#' Transpose the species data summaries
#'
#' @param m Rmonad
#' @export
get_summaries <- function(m){
  list(
    aa          = .extract(m, "summary_aa"),
    genome      = .extract(m, "summary_genome"),
    gff         = .extract(m, "summary_gff"),
    nstring     = .extract(m, "summary_nstring"),
    orffaa      = .extract(m, "summary_orffaa"),
    orfgff      = .extract(m, "summary_orfgff"),
    phase       = .extract(m, "summary_phase"),
    transfna    = .extract(m, "summary_transfna"),
    transorfaa  = .extract(m, "summary_transorfaa"),
    transorfgff = .extract(m, "summary_transorfgff")
  )
}

#' Get a table listing the scaffold count for each species
#'
#' @param m Rmonad
#' @export
number_of_scaffolds <- function(m){
  .extract(m, "summary_genome") %>% lapply(function(s) nrow(s@table)) %>%
    {data.frame(nscaffold = unlist(.))}
}

#' Get a table listing the first protein residue counts for each species
#'
#' @param m Rmonad
#' @export
initial_protein_residue_counts <- function(m){
  .extract(m, "summary_aa") %>%
    lapply(function(s) s@initial_residue) %>%
    lapply(function(x) {names(x)[which(names(x) == '*')] <- "STOP"; x}) %>%
    .rbindlist
}

#' Get a table listing the last protein residue proportions for each species
#'
#' @param m Rmonad
#' @export
final_protein_residue_counts <- function(m){
  .extract(m, "summary_aa") %>%
    lapply(function(s) s@final_residue) %>%
    lapply(function(x) {names(x)[which(names(x) == '*')] <- "STOP"; x}) %>%
    .rbindlist
}

#' Count the internal stop codons in each species protein models
#'
#' @param m Rmonad
#' @export
protein_has_internal_stop <- function(m){
  .extract(m, "summary_aa") %>%
    lapply(function(s) s@has_internal_stop) %>%
    lapply(sum)
}

#' Summarize the genomic composition of each species
#'
#' @param m Rmonad
#' @export
genomic_composition <- function(m){
  .extract(m, "summary_genome") %>%
    lapply(function(s) s@comp %>% colSums) %>%
    .rbindlist(sort=FALSE, remove_empty=TRUE)
}

#' Get the phylogenetic order of species (mostly useful for plotting)
#'
#' @param con Config object
#' @export
get_species_phylogenetic_order <- function(con){
  get_species(con)
}

#' Tabulate states for each genome
#'
#' @param m Rmonad
#' @param species_order all species in tree order
#' @export
make_genome_table <- function(m, species_order){
  nscaf <- number_of_scaffolds(m)
  initialRes <- initial_protein_residue_counts(m)
  finalRes <- final_protein_residue_counts(m)
  hasStop <- protein_has_internal_stop(m)
  genComp <- as.data.frame(genomic_composition(m))
  genComp$GC <- (genComp$G + genComp$C) / (genComp$G + genComp$C + genComp$A + genComp$T)

  species  <- factor(species_order, levels=species_order)
  scafs    <- as.matrix(nscaf)[species_order, ]
  bases    <- sapply(.extract(m, "summary_genome"), function(x) x@table$length %>% sum)
  prots    <- sapply(.extract(m, "summary_aa"), function(x) sum(x@initial_residue))
  GC       <- as.matrix(genComp)[species_order, 'GC']
  oddstart <- apply(initialRes, 1, function(x) { sum(x) - x['M'] })
  oddstop  <- apply(finalRes, 1, function(x) { sum(x) - x['STOP'] })
  p_N      <- as.matrix(genComp)[, "N"]
  n_stop   <- hasStop %>% unlist

  data.frame(
      species  = species
    , scafs    = scafs[species_order]
    , bases    = bases[species_order]
    , prots    = prots[species_order]
    , GC       = GC[species_order]
    , oddstart = oddstart[species_order]
    , oddstop  = oddstop[species_order]
    , p_N      = p_N[species_order]
    , n_stop   = n_stop[species_order]
  ) %>% dplyr::arrange(species)
}

#' Tabulate states for each genome
#'
#' @param m Rmonad
#' @export
make_synmap_table <- function(m){
  .extract(m, 'synmap_summary') %>%
    lapply(function(x) {
       w <- x@width
       values <- c(w@min, w@q25, w@median, w@q75, w@max, w@mean, w@sd, w@n)
       names(values) <- c("min", "q25", "median", "q75", "max", "mean", "sd", "N")
       values
    }) %>%
    .rbindlist(sort=FALSE)
}

#' Tabulate states for each genome
#'
#' @param m Rmonad object
#' @param genes character list of genes to subset the table by
#' @param group string - either "query" or "control"
#' @export
make_synder_table <- function(m, genes=NULL){
  .extract(m, 'synder_out') %>%
    # target and query are actually the same, so no need to do both
    .select_tag(paste0("/query$")) %>%
    lapply(function(x){
      if(!is.null(genes)){
        x <- x[S4Vectors::mcols(x)$attr %in% genes]
      }
      CNEr::second(x) %>%
      GenomicRanges::width() %>%
      fagin::summarize_numeric()
    }) %>%
      make_numeric_summary_table
}

#' Get tables of number of search intervals
#'
#' @param m Rmonad object
#' @export
get_synder_nsi <- function(m){
  xs <- .extract(m, 'synder_out') %>%
    # target and query are actually the same, so no need to do both
    .select_tag(paste0("/query$"))

  tabs <- lapply(xs, function(x){
    S4Vectors::mcols(x)$attr %>% as.factor %>% table %>%
      {data.frame(
        seqid = names(.),
        count = as.integer(.)
      )}
  })

  tabs <- lapply(names(tabs), function(n){
    d <- tabs[[n]]
    names(d)[2] <- n 
    d
  })

  .merge <- function(x, y){
    merge(x, y, by='seqid', all=TRUE)
  }

  tab <- Reduce(f=.merge, x=tabs[-1], init=tabs[[1]])
  tab[is.na(tab)] <- 0

  tab
}

#' Get summary table of synder search intervals per gene counts 
#'
#' @param m Rmonad object
#' @export
get_synder_count_table <- function(m, genes=NULL){
  x <- get_synder_nsi(m)
  if(! is.null(genes)){
    x <- x[x$seqid %in% genes, ]
  }
  x[-1] %>%
    lapply(summarize_numeric) %>%
    make_numeric_summary_table
}

make_query_target_table <- function(m){
  apply_names <- function(x, f){
      xnames <- names(x)
      x <- lapply(names(x), function(name) f(x[[name]], name))
      names(x) <- xnames
      x
  }

  # Gather binary features
  feats <- rmonad::get_value(m, tag='feature_table/query')
  names(feats) <- sub("feature_table/query/", "", names(feats))
  feats <- apply_names(feats, function(d, n) {d$target_species <- n; d}) 
  feats_table <- do.call(rbind, feats)
  rownames(feats_table) <- NULL
  # Gather class (e.g. UNO), primary, and secondary labels
  seclabels <- rmonad::get_value(m, tag="query_labels")[[1]]
  seclabels <- apply_names(seclabels$labels, function(d, n) {d$target_species <- n; d})
  sec_table <- do.call(rbind, seclabels)
  rownames(sec_table) <- NULL
  # Merge it
  big_table <- merge(feats_table, sec_table, by=c("seqid", "target_species"))
  # Add final class
  classStr <- rmonad::get_value(m, tag='query_origins')[[1]]$classStr
  big_table$class <- classStr[big_table$seqid]
  # Check it
  stopifnot(nrow(big_table) == nrow(feats_table))
  big_table
}

#' Summarize the synder flags for each species
#'
#' @param m Romand object
#' @param group string - either "query" or "control"
#' @export
make_synder_flag_table <- function(m, group="query"){
  .extract(m, 'summary_synder_flags') %>%
    .select_tag(paste0("/", group, "$")) %>%
    lapply(function(d){
        dplyr::select(d, -attr) %>% colSums
      }) %>%
    .rbindlist
}

label_desc <- data.frame(
   desc = c(
       "AA match to known protein"          , #O1
       "AA match to genic ORF"              , #O2
       "AA match to inter-genic ORF"        , #O3
       "maps off the scaffold"              , #U1
       "possible indel"                     , #U2
       "possible N-string"                  , #U3
       "syntenically scrambled"             , #U5
       "seriously unknown"                  , #U6
       "skipped for technical reasons"      , #U7
       "DNA match to CDS"                   , #N1
       "DNA match to exon"                  , #N2
       "DNA match to gene"                  , #N3
       "DNA match to inter-genic region"      #N4
   ),
   secondary = c(
       "O1" ,
       "O2" ,
       "O3" ,
       "U1" ,
       "U2" ,
       "U3" ,
       "U5" ,
       "U6" ,
       "U7" ,
       "N1" ,
       "N2" ,
       "N3" ,
       "N4"  
   )
)

#' Make a table of the secondary labels for each query gene
#'
#' @param m Rmonad object
#' @export
make_secondary_genewise_table <- function(m){
  labels <- .extract(m, "query_labels")[[1]]$labels
  species_order <- names(labels)
  lapply(
    names(labels),
    function(x) {
      dplyr::select(labels[[x]], seqid, secondary) %>%
        { names(.)[2] <- x; . }
    }
  ) %>% Reduce(f=merge) %>% { .[, c('seqid', species_order[-1])] } 
}

#' Summarize secondary labels as a table
#'
#' @param m Rmonad
#' @param species_order all target species in tree order
#' @export
make_secondary_labels_table <- function(m, species_order){
  parse_labels <- function(labels, group){
    labels %>%
      lapply(
        function(x) {
            dplyr::group_by(x[-1], primary, secondary) %>% dplyr::count()
        }
      ) %>%
      dplyr::bind_rows(.id='species') %>%
      {
        .$species <- factor(.$species, levels=species_order)
        . <- merge(., label_desc)
        .$desc <- paste(.$secondary, .$desc, sep=': ')
        .$group <- group
        .
      }
  }

  rbind(
    parse_labels(.extract(m, "query_labels")[[1]]$labels, "query"),
    parse_labels(.extract(m, "control_labels")[[1]]$labels, "control")
  )
}

#' Plot the secondary labels
#'
#' @param m Rmonad
#' @param species_order all target species in tree order
#' @param fill character - either 'secondary' or not
#' @export
plot_secondary_labels <- function(m, species_order, fill='secondary'){ 

  dat <- make_secondary_labels_table(m, species_order)

  if(fill == 'secondary'){
    ggplot2::ggplot(dat) +
      ggplot2::geom_bar(ggplot2::aes(x=species, y=n, fill=desc), position="dodge", stat="identity") +
      ggplot2::scale_fill_brewer(palette="Paired") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=325, hjust=0, vjust=1)) +
      ggplot2::labs(
        fill="Classification",
        x="Target species",
        y="# of focal genes"
      ) +
      ggplot2::facet_grid(group ~ .)
  } else {
    ggplot2::ggplot(dat) +
      ggplot2::scale_fill_brewer(palette="Paired") +
      ggplot2::geom_bar(ggplot2::aes(x=desc, y=n, fill=species), position="dodge", stat="identity")+
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=325, hjust=0, vjust=1)) +
      ggplot2::labs(
        fill="Target species",
        x="Classification",
        y="# of focal genes"
      ) +
      ggplot2::facet_grid(group ~ .)
  }
}

#' Make a table of the origin summaries for each gene
#'
#' @param m Rmonad object
#' @export
make_origin_table <- function(m){
  query <- get_value(m, tag='query_origins')[[1]]$classSum
  cntrl <- get_value(m, tag='control_origins')[[1]]$classSum
  tbl <- merge(query, cntrl, by='group', all=TRUE)
  tbl[is.na(tbl)] <- 0
  names(tbl) <- c('group', 'query', 'control')
  tbl
}


#' Organize output data into a small archive
#'
#' Running \code{run_fagin} will produce a large archive with cached data and
#' many intermediate files. This function digests that archive into a clean,
#' succinct, well-organized folder.
#'
#' @param con Config object
#' @return filename of the output archive
make_result_archive <- function(con){
  stop("This function is a stub")
}

#' Create an excel spreadsheet of fagin results
#'
#' Contains all the tabular data created by the \code{make_result_archive} function.
#'
#' @param con Config object
make_excel_spreadsheet <- function(m, filename="fagin-result.xlsx"){
  # TODO: if I can save the image as an EMF, it can be edited in Excel (at
  # least on windows), but currently there seems to be a bug in XLConnect (or
  # some dependency) that prevents this. See
  # https://github.com/miraisolutions/xlconnect issue #22

  ss <- get_summaries(m)
  species_order <- get_species_phylogenetic_order(con)

  wb <- XLConnect::loadWorkbook(filename, create=TRUE)

  gentab <- make_genome_table(ss, species_order)
  XLConnect::createSheet(wb, "Genome Summaries")
  XLConnect::writeWorksheet(wb, data=gentab, sheet="Genome Summaries")

  maptab <- make_synmap_table(con)
  XLConnect::createSheet(wb, "Synmap Summaries")
  XLConnect::writeWorksheet(wb, data=maptab, sheet="Synmap Summaries")

  syntab <- make_synder_table(con)
  XLConnect::createSheet(wb, "Synder Summaries")
  XLConnect::writeWorksheet(wb, data=syntab, sheet="Synder Summaries")

  key <- label_desc[, c(2,1)]
  names(key) <- c("label", "description")
  XLConnect::createSheet(wb, "Key")
  XLConnect::writeWorksheet(wb, data=key, sheet="Key")

  labtab <- make_secondary_table(con, species_order)
  XLConnect::createSheet(wb, "Labels")
  XLConnect::writeWorksheet(wb, data=labtab, sheet="Labels")

  dir <- dirname(filename)
  dir.create(file.path(dir, 'figures'), recursive=TRUE)

  fig1_path <- file.path(dir, 'figures', 'fig1.png')
  fig2_path <- file.path(dir, 'figures', 'fig2.png')

  png(filename=fig1_path)
    plot_secondary_labels(con, species_order)
  dev.off()
  XLConnect::createSheet(wb, "Fig1")
  XLConnect::createName(
    wb,
    name    = 'Fig1',
    formula = paste('Fig1', XLConnect::idx2cref(c(1, 1)), sep="!")
  )
  XLConnect::addImage(
    wb,
    filename     = fig1_path,
    name         = 'Fig1',
    originalSize = TRUE
  )

  png(filename=fig2_path)
    plot_secondary_labels(con, species_order, fill='species')
  dev.off()
  XLConnect::createSheet(wb, "Fig2")
  XLConnect::createName(
    wb,
    name    = 'Fig2',
    formula = paste('Fig2', XLConnect::idx2cref(c(1, 1)), sep="!")
  )
  XLConnect::addImage(
    wb,
    filename     = fig2_path,
    name         = 'Fig2',
    originalSize = TRUE
  )

  XLConnect::saveWorkbook(wb)
}
