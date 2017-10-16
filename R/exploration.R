#' Transpose the species data summaries
#'
#' @export
invertSummaries <- function(x){
  xs <- list(
    gff         = list(),
    dna         = list(),
    aa          = list(),
    trans       = list(),
    orfgff      = list(),
    orffaa      = list(),
    transorfgff = list(),
    transorffaa = list(),
    nstring     = list()
  )
  for(s in names(x@species)){
    xs$gff         = append(xs$gff         , from_cache(x@species[[s]]@summaries)@gff.summary         )
    xs$dna         = append(xs$dna         , from_cache(x@species[[s]]@summaries)@dna.summary         )
    xs$aa          = append(xs$aa          , from_cache(x@species[[s]]@summaries)@aa.summary          )
    xs$trans       = append(xs$trans       , from_cache(x@species[[s]]@summaries)@trans.summary       )
    xs$orfgff      = append(xs$orfgff      , from_cache(x@species[[s]]@summaries)@orfgff.summary      )
    xs$orffaa      = append(xs$orffaa      , from_cache(x@species[[s]]@summaries)@orffaa.summary      )
    xs$transorfgff = append(xs$transorfgff , from_cache(x@species[[s]]@summaries)@transorfgff.summary )
    xs$transorffaa = append(xs$transorffaa , from_cache(x@species[[s]]@summaries)@transorffaa.summary )
    xs$nstring     = append(xs$nstring     , from_cache(x@species[[s]]@summaries)@nstring.summary     )
  }
  lapply(xs, function(xx) { names(xx) <- names(x@species); xx } )
}

#' Get a table listing the scaffold count for each species
#'
#' @export
number_of_chromosomes <- function(x){
    lapply(x$dna, function(s) nrow(s@table))
}

#' Get a table listing the first protein residue proportions for each species
#'
#' @export
initial_protein_residue_counts <- function(x){
    lapply(x$aa, function(s) s@initial_residue)
}

#' Get a table listing the last protein residue proportions for each species
#'
#' @export
final_protein_residue_counts <- function(x){
    lapply(x$aa, function(s) s@final_residue)
}

#' Count the internal stop codons in each species protein models
#'
#' @export
protein_has_internal_stop <- function(x){
    lapply(x$aa, function(s) s@has_internal_stop %>% sum)
}

#' Summarize the genomic composition of each species
#'
#' @export
genomic_composition <- function(x){
    x$dna %>% lapply(function(s) s@comp %>% colSums %>% { . / sum(.) })
}

#' Get the phylogenetic order of species (mostly useful for plotting)
#'
#' @param con Config object
#' @export
getSpeciesPhylogeneticOrder <- function(con){
  # FIXME: loading everything like this is an ugly hack
  load(file.path(con@archive, 'd1.Rda'))
  d1@tree$tip.label
}

#' Tabulate states for each genome
#'
#' @param list  Summarized data as output by \code{invertSummaries}
#' @export
makeGenomeTable <- function(ss, speciesOrder){
  nscaf <- number_of_chromosomes(ss)
  initialRes <- initial_protein_residue_counts(ss)
  finalRes <- final_protein_residue_counts(ss)
  hasStop <- protein_has_internal_stop(ss)
  genComp <- genomic_composition(ss)
  data.frame(
      species = names(nscaf) %>% factor(levels=speciesOrder)
    , scafs = nscaf %>% unlist %>% unname
    , bases = sapply(ss$dna, function(x) x@table$length %>% sum)
    , prots = sapply(initialRes, sum) %>% unname
    , GC = sapply(genComp, function(x) ((x['G'] + x['C']) * 100) %>% round)
    , oddstart = sapply(initialRes, function(x) { sum(x) - x['M'] }) %>% unname
    , oddstop = sapply(finalRes, function(x) { sum(x) - x['*'] }) %>% unname
    , p_N = sapply(genComp, function(x) x['N'] %>% signif(3))
    , n_stop = hasStop %>% unlist %>% unname
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
rownamesAsSpecies <- function(x, levels){
  x$species <- rownames(x) %>% factor(levels=levels)
  rownames(x) <- NULL
  dplyr::arrange(x, species) %>%
    { .[c(ncol(.), 1:(ncol(.)-1))] } # put species first
}

#' Tabulate states for each genome
#'
#' @param con Config object
#' @export
makeSynmapTable <- function(con, speciesOrder){
  load(file.path(con@archive, 'd1.Rda'))
  d1@synmaps %>%
    lapply(function(x) x@synmap.summary@width) %>%
    makeNumericSummaryTable %>%
    rownamesAsSpecies(speciesOrder)
}

#' Tabulate states for each genome
#'
#' @param con Config object
#' @export
makeSynderTable <- function(con, speciesOrder){
  load('ARCHIVE/d2.Rda')
  d2$query %>%
    lapply(function(x){
       x$si %>%
           CNEr::second() %>%
           GenomicRanges::width() %>%
           fagin::summarize_numeric()
     }) %>%
       makeNumericSummaryTable %>%
       rownamesAsSpecies(speciesOrder)
}

#' Plot the secondary labels
#'
#' @export
plotSecondaryLabels <- function(con, fill='secondary'){
     desc <- data.frame(
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
    load(file.path(con@archive, 'd4.Rda'))

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
          . <- merge(., desc)
          .$desc <- paste(.$secondary, .$desc, sep=': ')
          .$group <- group
          .
        }
    }

    dat <- rbind(
      parse_labels(d4$query$labels,   group="orphan"),
      parse_labels(d4$control$labels, group="control")
    )

    if(fill == 'secondary'){
      ggplot(dat) +
        geom_bar(aes(x=species, y=n, fill=desc), position="dodge", stat="identity") +
        scale_fill_brewer(palette="Paired") +
        theme(axis.text.x = element_text(angle=270, hjust=0, vjust=1)) +
        labs(
          fill="Classification",
          x="Target species",
          y="# of focal genes"
        ) +
        facet_grid(group ~ .)
    } else {
      ggplot(dat) +
        scale_fill_brewer(palette="Paired") +
        geom_bar(aes(x=desc, y=n, fill=species), position="dodge", stat="identity")+
        theme(axis.text.x = element_text(angle=325, hjust=0, vjust=1)) +
        labs(
          fill="Target species",
          x="Classification",
          y="# of focal genes"
        ) +
        facet_grid(group ~ .)
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

#' Create an excel spreadsheet of fagin results
#'
#' Contains all the tabular data created by the \code{makeResultArchive} function.
#'
#' @export
#' @param con Config object
makeExcelSpreadsheet <- function(con, filename="fagin-result.xlsx"){
  # TODO: if I can save the image as an EMF, it can be edited in Excel (at
  # least on windows), but currently there seems to be a bug in XLConnect (or
  # some dependency) that prevents this. See
  # https://github.com/miraisolutions/xlconnect issue #22

  load(file.path(con@archive, "d5.Rda"))
  load(file.path(con@archive, "d1.Rda"))
  ss <- invertSummaries(d1)
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

  png(filename='fig1.png')
    plotSecondaryLabels(con)
  dev.off()
  XLConnect::createSheet(wb, "Fig1")
  XLConnect::createName(
    wb,
    name    = 'Fig1',
    formula = paste('Fig1', XLConnect::idx2cref(c(1, 1)), sep="!")
  )
  XLConnect::addImage(
    wb,
    filename     = 'fig1.png',
    name         = 'Fig1',
    originalSize = TRUE
  )

  png(filename='fig2.png')
    plotSecondaryLabels(con, fill='species')
  dev.off()
  XLConnect::createSheet(wb, "Fig2")
  XLConnect::createName(
    wb,
    name    = 'Fig2',
    formula = paste('Fig2', XLConnect::idx2cref(c(1, 1)), sep="!")
  )
  XLConnect::addImage(
    wb,
    filename     = 'fig2.png',
    name         = 'Fig2',
    originalSize = TRUE
  )

  XLConnect::saveWorkbook(wb)
}
