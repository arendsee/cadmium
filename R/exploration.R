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
        theme(
          axis.text.x = element_text(angle=270, hjust=0, vjust=1),
          legend.position = 'bottom'
        ) +
        guides(fill = guide_legend(ncol = 2)) +
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
        theme(
          axis.text.x = element_text(angle=325, hjust=0, vjust=1)
        ) +
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
  load(file.path(con@archive, "d5.Rda"))
  wb <- XLConnect::loadWorkbook(filename, create=TRUE)

  s1 <- "query_origins"
  t1 <- d5$query
  test_img <- "test_img"

  # TODO: if I can save the image as an EMF, it can be edited in Excel (at
  # least on windows), but currently there seems to be a bug in XLConnect (or
  # some dependency) that prevents this. See
  # https://github.com/miraisolutions/xlconnect issue #22

  # test_img_filename <- paste0(test_img, ".emf")
  test_img_filename <- paste0(test_img, ".png")

  XLConnect::createSheet(wb, s1)
  XLConnect::writeWorksheet(wb, data=t1, sheet=s1) 

  # devEMF::emf(file=test_img_filename)
  png(file=test_img_filename)
  qplot(rnorm(100))
  dev.off()

  XLConnect::createName(
    wb,
    name    = test_img,
    formula = paste(s1, XLConnect::idx2cref(c(5, ncol(t1)) + 2), sep="!")
  )
  XLConnect::addImage(
    wb,
    filename     = test_img_filename,
    name         = test_img,
    originalSize = TRUE
  )

  XLConnect::saveWorkbook(wb)
}
