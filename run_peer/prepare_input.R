library(optparse)

option_list <- list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="Input data matrix (example: sample x gene matrix)",
                metavar="character"),
    make_option(c("-t", "--if_transpose"), type="character", default='No',
                help="'No'/'Yes'. If we should transpose the data matrix before applying PEER (example 'No'; default 'No')",
                metavar="character"),
    make_option(c("-s", "--skip_cols"), type="character", default=NULL,
                help="Skip columns (separated by ',')",
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="output CSV (X.csv in peer)",
                metavar="character"),
    make_option(c("-n", "--output_indiv_list"), type="character", default=NULL,
                help="output individual list (the rownames)",
                metavar="character"),
    make_option(c("-d", "--indiv_col"), type="character", default=NULL,
                help="column name of individual id column. It will be used only when if_transpose = No",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# This script prepare expression matrix for PEER factor workflow suggested here: https://www.nature.com/articles/nprot.2011.457

library(data.table)

myread = function(fname) {
  if(tools::file_ext(fname) == 'gz') {
    exec = 'zcat'
  } else {
    exec = 'cat'
  }
  df = fread(paste(exec, '<', fname), header = T, sep = '\t', data.table = F)
  return(df)
}

skip_cols = function(df, col_str) {
  if(is.null(col_str)) {
    return(df)
  }
  cols = strsplit(col_str, ',')[[1]]
  df = df[, ! colnames(df) %in% cols]
  return(df)
}

take_care_of_transpose = function(df, if_transpose) {
  if(if_transpose == 'No') {
    message('No transformation so assuming data matrix is N x G')
  } else if(if_transpose == 'Yes') {
    message('No transformation so assuming data matrix is G x N')
    df = t(df)
  } else {
    message('Not supported if_transpose = ', if_transpose)
  }
  return(df)
}

mat = myread(opt$input)

if(opt$if_transpose == 'No') {
  indiv_list = as.character(mat[, opt$indiv_col])
}

mat = skip_cols(mat, opt$skip_cols)

if(opt$if_transpose == 'Yes') {
  indiv_list = colnames(mat)
}

mat = take_care_of_transpose(mat, opt$if_transpose)
message('Dimension of data matrix: ', 'N = ', dim(mat)[1], '; G = ', dim(mat)[2])

## Impute missing cells as suggested in paper
mat_imp = impute::impute.knn(as.matrix(t(mat)))
write.table(t(as.matrix(mat_imp$data)), opt$output, row = F, col = F, quo = F, sep = ',')

write.table(data.frame(x = indiv_list), opt$output_indiv_list, col = F, row = F, quo = F)
