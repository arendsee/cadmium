#' Build a query to target map for all intervals in the target GFF
analyzeTargetFeature <- function(query, target){
  require(GenomicRanges)

  si <- target$si$target

  # I assume the order of elements in the query and search interval GRange
  # objects is the same. The  *id* column was set at load time to be the same
  # in the two objects. Here I assert that they have not been scrambled.
  stopifnot(mcols(si)$id == mcols(target$si$query)$id)

  # A Hits object
  # from(o) accesses the si indices
  # to(o) acceses the ft indices
  o <- findOverlaps(si, target$gff)

  # NOTE: There is a teeny-tiny difference between the results I previously got
  # with Genome_intervals and what I now get. Using GRanges I get 103306 rows
  # in query2target for Arabidopsis thaliana versus lyrata data, with
  # Genomic_intervals (and the old code), I got 103274. Maybe a small
  # difference in how the algorithms are implemented? But I'm not sure.
  data.frame(
    query=mcols(target$si$query)$seqid[from(o)],
    target=mcols(target$gff)$seqid[to(o)],
    type=mcols(target$gff)$type[to(o)],
    stringsAsFactors=FALSE
  ) %>%
    unique %>%
    {split(., .$type)} %>%
    lapply(function(x) x[1:2])
}

findQueryGaps <- function(nstring, target){

  sitar <- target$si$target
  sique <- target$si$query

  stopifnot(mcols(sitar)$id == mcols(sique)$id)

  # Find overlaps between search intervals and N-strings
  over <- findOverlaps(sitar, nstring)

  stopifnot(max(from(over)) <= length(sitar))
  stopifnot(max(to(over)) <= length(nstring))

  data.frame(
    query = sique$seqid[from(over)],
    length = nstring[to(over)] %>% width,
    stringsAsFactors=FALSE
  )
}
