#' Collate the warnings of a Biostrings::translate function
#'
#' The Biostrings::translate function emits one warning for each CDS that is
#' not a multiple of 3. This function collates all these warnings into one.
#' Also it specifies the name of each sequence and the species from whch they
#' derive.
#'
#' The warning is assumed to have the form:
#' "in 'x[[123]]': last 2 bases were ignored"
#' "in 'x[[456]]': last base was ignored"
#'
#' @param node Rmonad object wrapping the result of a Biostrings::translate function
#' @param species The species name
#' @return Rmonad An Rmonad object with a modified warning field
collate_translation_warnings <- function(node, species){
  if(m_OK(node)){
    node_ids <- which(grepl("base.*ignored", rmonad::m_warnings(node)))
    model_ids <- as.integer(sub(
      "in 'x\\[\\[(\\d+)\\]\\]': .*base.*ignored",
      "\\1",
      rmonad::m_warnings(node)[node_ids],
      perl=TRUE
    ))
    if(length(model_ids) > 0){ 
      msg <- "%s of %s gene models in %s are truncated (CDS length is not a multiple of 3): [%s]"
      msg <- sprintf(
        msg,
        length(model_ids),
        length(m_value(node)),
        species,
        paste0(names(m_value(node)[model_ids]), collapse=', ')
      )
      rmonad::m_warnings(node) <- c(msg, rmonad::m_warnings(node)[-node_ids])
    }
  }
  node
}
