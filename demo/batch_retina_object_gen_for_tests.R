
load_function_and_run <- function(string_identifier){
  path_to_file <- paste0('retina/extdata/test_retinas/',
                     string_identifier,
                     '/',
                     'main.r')
  path_to_retina_dir <- paste0('retina/extdata/test_retinas/',
                       string_identifier,
                       '/diagram_retina/')
  #source the file from the installed packages directory
  path_to_test_retina_folder <- file.path(.libPaths(), path_to_retina_dir)
  browser()
  source(file.path(.libPaths(), path_to_file))
  run_command <- paste0(string_identifier, "(", path_to_retina_dir, ")")
  eval(parse(run_command))
}

message('Batch computing retina objects')
list_of_retina_identifiers <- c('3hgbqg',
'bjii4t',
'jhvbwt',
'40oik5',
'8pm223',
'gv3igs',
'raxp91')

lapply(list_of_retina_identifiers, load_function_and_run)
