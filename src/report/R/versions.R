get_synder_version <- function(){
  synder_path <- file.path(config$d_home, 'bin', 'synder')
  system2(synder_path, args=list('-v'), stdout=TRUE)
}

get_fagin_version <- function(config){
  # TODO: Make this portable
  system2('cat', args=list(file.path(config$d_home, 'VERSION')), stdout=TRUE) 
}
