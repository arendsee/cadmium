# Currently all of these simply save the respective input as Rdata. Eventually
# I may specialize them.

archive_1 <- function(d1, archive_dir){

  "Archive data from the first analysis - gathering and checking raw data."

  save(d1, file=file.path(archive_dir, "d1.Rda"))

}

archive_2 <- function(d2, archive_dir){

  "Archive data from the second analysis - comparisons of target species to the
  focal species."

  save(d2, file=file.path(archive_dir, "d2.Rda"))

}

archive_3 <- function(d3, archive_dir){

  "Archive data from the third analysis - building a feature table from the
  target/focal comparisons."

  save(d3, file=file.path(archive_dir, "d3.Rda"))

}

archive_4 <- function(d4, archive_dir){

  "Archive data from the fourth analysis - labels and trees built from the
  feature table."

  save(d4, file=file.path(archive_dir, "d4.Rda"))

}

archive_5 <- function(d5, archive_dir){

  "Archive data from the fifth analysis - superimposing labels on the
  phylogenetic tree."

  save(d5, file=file.path(archive_dir, "d5.Rda"))

}

archive_rmonad <- function(m, archive_dir){

  "Archive the final rmonad object"

  save(m, file=file.path(archive_dir, "final_obj.Rda"))

}
