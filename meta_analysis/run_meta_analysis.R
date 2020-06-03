library(optparse)

option_list <- list(
    make_option(c("-G", "--gwas1"), type="character", default=NULL,
                help="input GWAS",
                metavar="character"),
    make_option(c("-a", "--gwas2"), type="character", default=NULL,
                help="column name of p-value",
                metavar="character"),
    make_option(c("-n", "--sample_size_cols"), type="character", default=NULL,
                help="list sample size column names of gwas1 and gwas2 separated by ,",
                metavar="character"),
    make_option(c("-z", "--zscore_cols"), type="character", default=NULL,
                help="list z-score column names of gwas1 and gwas2 separated by ,",
                metavar="character"),
    make_option(c("-b", "--beta_cols"), type="character", default=NULL,
                help="list effect size column names of gwas1 and gwas2 separated by ,",
                metavar="character"),
    make_option(c("-s", "--se_cols"), type="character", default=NULL,
                help="list standard deviation column names of gwas1 and gwas2 separated by ,",
                metavar="character"),
    make_option(c("-t", "--other_cols_in_1"), type="character", default=NULL,
                help="list other column names of gwas1 to keep",
                metavar="character"),
    make_option(c("-e", "--other_cols_in_2"), type="character", default=NULL,
                help="list other column names of gwas2 to keep",
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="output with new additional column indicating effect size and standard error converted from OR and p-value",
                metavar="character"),
    make_option(c("-u", "--source_file_dir"), type="character", default='.',
                help="specify the directory of source file",
                metavar="character"),
    make_option(c("-j", "--joinby_col"), type="character", default=NULL,
                help="specify column names of gwas1 and gwas2 to join by",
                metavar="character"),
    make_option(c("-f", "--af_col_in_1"), type="character", default=NULL,
                help="specify allele-frequency column names of gwas1 to use for convert z-score to effect size",
                metavar="character"),
    make_option(c("-x", "--af_col_in_2"), type="character", default=NULL,
                help="specify allele-frequency column names of gwas2 to use for convert z-score to effect size",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

source(paste0(opt$source_file_dir, '/meta_analysis.R'))
source(paste0(opt$source_file_dir, '/util.R'))

library(dplyr)

joinbys = read_cols(opt$joinby_col)

if(length(joinbys) > 0) {
  df1 = myread(opt$gwas1)
  df2 = myread(opt$gwas2)
  df1$joinbycol = df1[, joinbys[1]]
  df2$joinbycol = df2[, joinbys[2]]
}



if (!is.null(opt$sample_size_cols) & !is.null(opt$zscore_cols)) {
  message('Detected both sample size and zscore columns, run sample size based meta-analysis ... ')
  zcols = read_cols(opt$zscore_cols)
  ncols = read_cols(opt$sample_size_cols)
  df1$myzscore1 = df1[, zcols[1]]
  df2$myzscore2 = df2[, zcols[2]]
  df1$myN1 = df1[, ncols[1]]
  df2$myN2 = df2[, ncols[2]]
  keepcol1 = c('joinbycol', 'myzscore1', 'myN1', read_cols(opt$other_cols_in_1))
  keepcol2 = c('joinbycol', 'myzscore2', 'myN2', read_cols(opt$other_cols_in_2))
  if (!is.null(args$af_col_in_1)) {
    df1$myaf = df1[, args$af_col_in_1]
    keepcol1 = c(keepcol1, 'myaf')
  } else if (!is.null(args$af_col_in_2)) {
    df2$myaf = df2[, args$af_col_in_2]
    keepcol2 = c(keepcol2, 'myaf')
  } else {
    message('No AF is specified so that z-score to effect size conversion will NOT be performed')
  }
  df1 = df1[, keepcol1]
  df2 = df2[, keepcol2]
  df = inner_join(df1, df2, by = 'joinbycol', suffix = c('.gwas1', '.gwas2'))
  if (nrow(df) == 0) {
    message('No rows after inner join. Exit ..')
    quit()
  }
  result = sample_size_based_meta(df$myzscore1, df$myzscore2, df$myN1, df$myN2)
  result = result %>% mutate(p_meta = zval2pval(z_meta))
  if('myaf' %in% colnames(df)) {
    result$effect_size_meta = zval2beta(df$myaf, df$myN1 + df$myN2, result$z_meta)
    result = result %>% mutate(se_meta = effect_size_meta / z_meta)
  }
} else if (!is.null(opt$beta_cols) & !is.null(opt$se_cols)) {
  message('Detected both effect size and SE columns, run inverse variance based meta-analysis ... ')
  bcols = read_cols(opt$beta_cols)
  scols = read_cols(opt$se_cols)
  df1$myeffectsize1 = df1[, bcols[1]]
  df2$myeffectsize2 = df2[, bcols[2]]
  df1$myse1 = df1[, scols[1]]
  df2$myse2 = df2[, scols[2]]
  keepcol1 = c('joinbycol', 'myeffectsize1', 'myse1', read_cols(opt$other_cols_in_1))
  keepcol2 = c('joinbycol', 'myeffectsize2', 'myse2', read_cols(opt$other_cols_in_2))
  df1 = df1[, keepcol1]
  df2 = df2[, keepcol2]
  df = inner_join(df1, df2, by = 'joinbycol', suffix = c('.gwas1', '.gwas2'))
  if (nrow(df) == 0) {
    message('No rows after inner join. Exit ..')
    quit()
  }
  result = inverse_variance_based_meta(df$myeffectsize1, df$myeffectsize2, df$myse1, df$myse2)
  result = result %>% mutate(z_meta = effect_size_meta / se_meta)
  result = result %>% mutate(p_meta = zval2pval(z_meta))
} else {
  message('Cannot determine which method to use. Exit ..')
  quit()
}

gz1 <- gzfile(opt$output, "w")
write.table(cbind(df, result), gz1, col = T, row = F, quo = F, sep = '\t')
close(gz1)

