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
.select_tag <- function(xs, tag){
  xs <- xs[grepl(pattern=pattern, x=names(xs), perl=TRUE)]
  names(xs) <- sub(pattern=pattern, replacement="", x=names(xs))
  xs
}

# Create a table from a named list of named vectors
#
# @param xs named list of named vectors
# @param fill ANY - value to replace missing cells with
# @param sort logical - sort the columns
# @param remove_empty logical - should columns with no data be removed?
.rbindlist <- function(xs, fill=0, sort=TRUE, remove_empty=TRUE){
  chars <- lapply(xs, names) %>% unlist %>% unique
  if(sort)
    chars <- sort(chars)
  tbl <- do.call(rbind, lapply(xs, function(x) x[match(chars, names(x))]))
  colnames(tbl) <- chars
  tbl[is.na(tbl)] <- fill
  tbl
}

# Calculate GC content from a named vector of counts 
.get_gc <- function(x){
  (x['G'] + x['C']) / (x['G'] + x['C'] + x['A'] + x['T'])
}


#' Transpose the species data summaries
#'
#' @export
getSummaries <- function(m){
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
#' @export
number_of_scaffolds <- function(m){
  .extract(m, "summary_genome") %>% lapply(function(s) nrow(s@table)) %>%
    {data.frame(nscaffold = unlist(.))}
}

#' Get a table listing the first protein residue counts for each species
#'
#' @export
initial_protein_residue_counts <- function(x){
  .extract(m, "summary_aa") %>%
    lapply(function(s) s@initial_residue) %>%
    lapply(function(x) {names(x)[which(names(x) == '*')] <- "STOP"; x}) %>%
    .rbindlist
}

#' Get a table listing the last protein residue proportions for each species
#'
#' @export
final_protein_residue_counts <- function(m){
  .extract(m, "summary_aa") %>%
    lapply(function(s) s@final_residue) %>%
    lapply(function(x) {names(x)[which(names(x) == '*')] <- "STOP"; x}) %>%
    .rbindlist
}

#' Count the internal stop codons in each species protein models
#'
#' @export
protein_has_internal_stop <- function(x){
  .extract(m, "summary_aa") %>%
    lapply(function(s) s@has_internal_stop) %>%
    lapply(sum)
}

#' Summarize the genomic composition of each species
#'
#' @export
genomic_composition <- function(x){
  .extract(m, "summary_genome") %>%
    lapply(function(s) s@comp %>% colSums) %>%
    .rbindlist(sort=FALSE)
}

#' Get the phylogenetic order of species (mostly useful for plotting)
#'
#' @param con Config object
#' @export
getSpeciesPhylogeneticOrder <- function(con){
  get_species(con)
}

#' Tabulate states for each genome
#'
#' @param list  Summarized data as output by \code{invertSummaries}
#' @export
makeGenomeTable <- function(m, speciesOrder){
  nscaf <- number_of_scaffolds(m)
  initialRes <- initial_protein_residue_counts(m)
  finalRes <- final_protein_residue_counts(m)
  hasStop <- protein_has_internal_stop(m)
  genComp <- as.data.frame(genomic_composition(m))
  genComp$GC <- (genComp$G + genComp$C) / (genComp$G + genComp$C + genComp$A + genComp$T)

  species  <- factor(speciesOrder, levels=speciesOrder)
  scafs    <- as.matrix(nscaf)[speciesOrder, ]
  bases    <- sapply(.extract(m, "summary_genome"), function(x) x@table$length %>% sum)
  prots    <- sapply(.extract(m, "summary_aa"), function(x) sum(x@initial_residue))
  GC       <- as.matrix(genComp)[speciesOrder, 'GC']
  oddstart <- apply(initialRes, 1, function(x) { sum(x) - x['M'] })
  oddstop  <- apply(finalRes, 1, function(x) { sum(x) - x['STOP'] })
  p_N      <- as.matrix(genComp)[, "N"]
  n_stop   <- hasStop %>% unlist

  data.frame(
      species  = species
    , scafs    = scafs[speciesOrder]
    , bases    = bases[speciesOrder]
    , prots    = prots[speciesOrder]
    , GC       = GC[speciesOrder]
    , oddstart = oddstart[speciesOrder]
    , oddstop  = oddstop[speciesOrder]
    , p_N      = p_N[speciesOrder]
    , n_stop   = n_stop[speciesOrder]
  ) %>% dplyr::arrange(species)
}

# Summarize a list numeric vectors
makeNumericSummaryTable <- function(nd){
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
rownamesAsSpecies <- function(x, levels=unique(rownames(x))){
  x$species <- rownames(x) %>% factor(levels=levels)
  rownames(x) <- NULL
  dplyr::arrange(x, species) %>%
    { .[c(ncol(.), 1:(ncol(.)-1))] } # put species first
}

#' Tabulate states for each genome
#'
#' @param con Config object
#' @export
makeSynmapTable <- function(m){
  .extract(m, 'synmap_summary') %>%
    lapply(function(x) {
       x <- x@width
       values <- c(x@min, x@q25, x@median, x@q75, x@max, x@mean, x@sd)
       names(values) <- c("min", "q25", "median", "q75", "max", "mean", "sd")
       values
    }) %>%
    .rbindlist(sort=FALSE)
}

#' Tabulate states for each genome
#'
#' @param m Rmonad object
#' @param group string - either "query" or "control"
#' @export
makeSynderTable <- function(m, group="query"){
  .extract(m, 'synder_out') %>%
    .select_tag(paste0("/", group, "$")) %>%
    lapply(function(x){
       CNEr::second(x) %>%
       GenomicRanges::width() %>%
       fagin::summarize_numeric()
     }) %>%
       makeNumericSummaryTable
}

#' Summarize the synder flags for each species
#'
#' @param m Romand object
#' @param group string - either "query" or "control"
makeSynderFlagTable <- function(m, group="query"){
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

#' Plot the secondary labels
#'
#' @export
plotSecondaryLabels <- function(m, speciesOrder, fill='secondary'){ 
    parse_labels <- function(labels, group){
      labels %>%
        lapply(
          function(x) {
              dplyr::group_by(x[-1], primary, secondary) %>% dplyr::count()
          }
        ) %>%
        dplyr::bind_rows(.id='species') %>%
        {
          .$species <- factor(.$species, levels=speciesOrder)
          . <- merge(., label_desc)
          .$desc <- paste(.$secondary, .$desc, sep=': ')
          .$group <- group
          .
        }
    }

    dat <- rbind(
      parse_labels(.extract(m, "query_labels")[[1]]$labels, "query"),
      parse_labels(.extract(m, "control_labels")[[1]]$labels, "control")
    )

    if(fill == 'secondary'){
      ggplot2::ggplot(dat) +
        ggplot2::geom_bar(ggplot2::aes(x=species, y=n, fill=desc), position="dodge", stat="identity") +
        ggplot2::scale_fill_brewer(palette="Paired") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle=270, hjust=0, vjust=1)) +
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

#' Organize output data into a small archive
#'
#' Running \code{run_fagin} will produce a large archive with cached data and
#' many intermediate files. This function digests that archive into a clean,
#' succinct, well-organized folder.
#'
#' @param con Config object
#' @return filename of the output archive
makeResultArchive <- function(con){
  stop("This function is a stub")
}


#' Make a table of the secondary labels for each query gene
#'
#' @param m Rmonad object
#' @return
makeSecondaryTable <- function(m){
  labels <- .extract(m, "query_labels")[[1]]$labels
  speciesOrder <- names(labels)
  lapply(
    names(labels),
    function(x) {
      dplyr::select(labels[[x]], seqid, secondary) %>%
        { names(.)[2] <- x; . }
    }
  ) %>% Reduce(f=merge) %>% { .[, c('seqid', speciesOrder[-1])] } 
}

#' Create an excel spreadsheet of fagin results
#'
#' Contains all the tabular data created by the \code{makeResultArchive} function.
#'
#' @param con Config object
makeExcelSpreadsheet <- function(m, filename="fagin-result.xlsx"){
  # TODO: if I can save the image as an EMF, it can be edited in Excel (at
  # least on windows), but currently there seems to be a bug in XLConnect (or
  # some dependency) that prevents this. See
  # https://github.com/miraisolutions/xlconnect issue #22

  ss <- getSummaries(m)
  speciesOrder <- getSpeciesPhylogeneticOrder(con)

  wb <- XLConnect::loadWorkbook(filename, create=TRUE)

  gentab <- makeGenomeTable(ss, speciesOrder)
  XLConnect::createSheet(wb, "Genome Summaries")
  XLConnect::writeWorksheet(wb, data=gentab, sheet="Genome Summaries")

  maptab <- makeSynmapTable(con, speciesOrder)
  XLConnect::createSheet(wb, "Synmap Summaries")
  XLConnect::writeWorksheet(wb, data=maptab, sheet="Synmap Summaries")

  syntab <- makeSynderTable(con, speciesOrder)
  XLConnect::createSheet(wb, "Synder Summaries")
  XLConnect::writeWorksheet(wb, data=syntab, sheet="Synder Summaries")

  key <- label_desc[, c(2,1)]
  names(key) <- c("label", "description")
  XLConnect::createSheet(wb, "Key")
  XLConnect::writeWorksheet(wb, data=key, sheet="Key")

  labtab <- makeSecondaryTable(con, speciesOrder)
  XLConnect::createSheet(wb, "Labels")
  XLConnect::writeWorksheet(wb, data=labtab, sheet="Labels")

  dir <- dirname(filename)
  dir.create(file.path(dir, 'figures'), recursive=TRUE)

  fig1_path <- file.path(dir, 'figures', 'fig1.png')
  fig2_path <- file.path(dir, 'figures', 'fig2.png')

  png(filename=fig1_path)
    plotSecondaryLabels(con, speciesOrder)
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
    plotSecondaryLabels(con, speciesOrder, fill='species')
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
