library(data.table)
options(datatable.fread.datatable = F)

read_cols = function(cols) {
  if(!is.null(cols)) {
    strsplit(cols, ',')[[1]]
  }
  else {
    NULL
  }
}

myread = function(filename) {
  message(tools::file_ext(filename))
  if(tools::file_ext(filename) == 'gz') {
    cmd = paste0('zcat < ', filename)
  } else {
    cmd = paste0(filename)
  }
  fread(cmd, header = T)
}
