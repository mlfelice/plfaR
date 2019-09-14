# function for checking correct file_type at start of processing funs
check_format <- function(x){
  if (class(x) == 'data.frame') return('data.frame')
  else if (class(x) == 'list') return('list')
  else stop('Input must be a dataframe or list of datframes')
}
