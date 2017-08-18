run_fagin <- function(con){
  primary_data(con) %>>%
    secondary_data(con) %>>%
    tertiary_data(con) %>>%
    determine_labels(con) %>>%
    determine_origins(con)
}
