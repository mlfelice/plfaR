# function for checking correct file_type at start of processing funs
check_format <- function(x){
  if ('data.frame' %in% class(x)) return('data.frame')
  else if (class(x) == 'list') return('list')
  else stop('Input must be a dataframe or list of datframes')
}

Check_duplicate_lipids <- function(df){
  if(nrow(unique(df[, c('DataFileName', 'Name')])) != nrow(df)){
    stop('Warning: Duplicate lipids detected in at least on sample\n',
         'Correct peak identification before proceeding')
  }
}

check_sample_match <- function(df1, df2){

}
