# #' Assert that x is what it ought to be
# #'
# #' @param x a thing to check
# #' @param type the general type is ought to have
# #' @param ... Additional arguments passed to specific type checkers
# #' @return The thing, raising warnings or errors as appropriate
# #' @export
# typecheck <- function(x, type, ...){
#
# }
#
# .typecheck_SpeciesName
# .typecheck_FolderPath
# .typecheck_GenomeDb
# .typecheck_GenomeSeq
# .typecheck_TranscriptomeDb
# .typecheck_TranscriptomeSeq
# .typecheck_CdsDna
# .typecheck_GenomeSeqinfo
# .typecheck_Length
# .typecheck_Count
# .typecheck_NString
# .typecheck_OrfRanges
# .typecheck_GeneModelsDB
# .typecheck_GeneModelList
# .typecheck_NStringSummary
# .typecheck_Phase
# .typecheck_DNASummary
# .typecheck_AASummary
# .typecheck_GFFSummary
# .typecheck_NumericSummary
# .typecheck_GRangesSummary

# internal function, mustly for testing
top_class <- function(x){
  class(x)[1]
}
