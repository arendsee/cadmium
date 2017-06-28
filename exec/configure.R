#!/usr/bin/Rscript --vanilla

deps <- c(
    "argparse",
    "knitr",
    "ape",
    "data.tree",
    "devtools",
    "dplyr",
    "fitdistrplus",
    "ggplot2",
    "intervals",
    "magrittr",
    "reshape2",
    "robustreg",
    "scales",
    "tidyr",
    "gridExtra",
    "xtable",
    "yaml"
)
bio.deps <- c('Biostrings', 'GenomicRanges')
current <- installed.packages()

chooseCRANmirror(ind=1)
for (x in c(deps, bio.deps)) {
    cat(paste("Checking R package", x, "... "))
    if(! x %in% current){
        cat("installing\n")
        if(x %in% bio.deps){
            require(devtools, quietly=TRUE, warn.conflicts=FALSE)
            source("https://bioconductor.org/biocLite.R")
            biocLite(x)
        } else {
            install.packages(x)
        }
    } else {
        cat("OK\n")
    }
}
