config_get_tree_file <- function(config) {
  if(class(config) == "config_input") {
    return(config@tree)
  } else if(class(config) == "fagin_config"){
    return(config@input@tree)
  } else {
    stop("Cannot get tree file from object of class '%s'", class(config))
  }
}
