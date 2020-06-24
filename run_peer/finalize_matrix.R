library(optparse)

option_list <- list(
    make_option(c("-i", "--input_x"), type="character", default=NULL,
                help="X.csv from peertool run",
                metavar="character"),
    make_option(c("-n", "--indiv_list"), type="character", default=NULL,
                help="individual id (row name of X.csv)",
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="PEER factor x individual matrix",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# This script prepare expression matrix for PEER factor workflow suggested here: https://www.nature.com/articles/nprot.2011.457

library(data.table)
options(datatable.fread.datatable = FALSE)

peer = fread(opt$input_x, sep = ',', header = FALSE)
indiv_list = as.character(read.table(opt$indiv_list, header = FALSE)$V1)

peer = cbind(indiv_list, peer)
colnames(peer) = c('individual_id', paste0('PEER', 1 : (ncol(peer) - 1)))

gz1 <- gzfile(opt$output, "w")
write.table(peer, gz1, row = FALSE, col = TRUE, quo = F, sep = '\t')
close(gz1)
