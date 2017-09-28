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
    xs$gff         = append(xs$gff         , x@species[[s]]@summaries@gff.summary         )
    xs$dna         = append(xs$dna         , x@species[[s]]@summaries@dna.summary         )
    xs$aa          = append(xs$aa          , x@species[[s]]@summaries@aa.summary          )
    xs$trans       = append(xs$trans       , x@species[[s]]@summaries@trans.summary       )
    xs$orfgff      = append(xs$orfgff      , x@species[[s]]@summaries@orfgff.summary      )
    xs$orffaa      = append(xs$orffaa      , x@species[[s]]@summaries@orffaa.summary      )
    xs$transorfgff = append(xs$transorfgff , x@species[[s]]@summaries@transorfgff.summary )
    xs$transorffaa = append(xs$transorffaa , x@species[[s]]@summaries@transorffaa.summary )
    xs$nstring     = append(xs$nstring     , x@species[[s]]@summaries@nstring.summary     )
  }
  lapply(xs, function(xx) { names(xx) <- names(x@species); xx } )
}

number_of_chromosomes <- function(x){
    lapply(x$dna, function(s) nrow(s@table))
}

initial_protein_residue_counts <- function(x){
    lapply(x$aa, function(s) s@initial_residue)
}

final_protein_residue_counts <- function(x){
    lapply(x$aa, function(s) s@final_residue)
}

protein_has_internal_stop <- function(x){
    lapply(x$aa, function(s) s@has_internal_stop %>% sum)
}

genomic_composition <- function(x){
    x$dna %>% lapply(function(s) s@comp %>% colSums %>% { . / sum(.) })
}

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
    dat <- d4$labels %>%
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
                .
            }
    if(fill == 'secondary'){
        ggplot(dat) +
            geom_bar(aes(x=species, y=n, fill=desc), position="dodge", stat="identity") +
            scale_fill_brewer(palette="Paired") +
            theme(
                axis.text.x = element_text(angle=270, hjust=0, vjust=1)
              , legend.position = 'bottom'
            ) +
            guides(fill = guide_legend(ncol = 2)) +
            labs(
                fill="Classification",
                x="Target species",
                y="# of focal genes"
            )
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
            )
    }
}
