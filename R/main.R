run_fagin <- function(con){
  primary_data(con) %>>% secondary_data(con)
}
