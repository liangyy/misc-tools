library(dplyr)
baseline_lf_rds = '~/Desktop/tmp/playground/ld4m/baselineLF.snp.rds'

# load UK10K SNPs
if(!file.exists(baseline_lf_rds)) {
  ss = list()
  
  for(i in 1 : 22) {
    # dd[[length(dd) + 1]] = data.table::fread(paste0('~/Downloads/1000G_EUR_Phase3_plink/1000G.EUR.QC.', i, '.bim'), data.table = F)
    # ff[[length(ff) + 1]] = data.table::fread(paste0('~/Downloads/1000G_Phase3_frq/1000G.EUR.QC.', i, '.frq'), data.table = F)
    # nn[[length(nn) + 1]] = data.table::fread(paste0('~/Downloads/baselineLD_v1.1/baselineLD.', i, '.l2.M_5_50'), data.table = F)
    message(i)
    ss[[length(ss) + 1]] = data.table::fread(paste0('~/Downloads/baselineLF/baselineLF.', i, '.annot.gz'), data.table = F) %>% select(MAFbin_lowfrq_1, MAFbin_lowfrq_2, MAFbin_lowfrq_3, MAFbin_lowfrq_4, MAFbin_lowfrq_5, MAFbin_common_1, MAFbin_common_2, MAFbin_common_3, MAFbin_common_4, MAFbin_common_5, MAFbin_common_6, MAFbin_common_7, MAFbin_common_8, MAFbin_common_9, MAFbin_common_10, SNP) %>% 
      mutate(
        common = MAFbin_common_1 + MAFbin_common_2 + MAFbin_common_3 + MAFbin_common_4 + MAFbin_common_5 + MAFbin_common_6 + MAFbin_common_7 + MAFbin_common_8 + MAFbin_common_9 + MAFbin_common_10,
        lowfrq = MAFbin_lowfrq_1 + MAFbin_lowfrq_2 + MAFbin_lowfrq_3 + MAFbin_lowfrq_4 + MAFbin_lowfrq_5) %>% 
      select(SNP, common, lowfrq)
    # mm[[length(mm) + 1]] = data.table::fread(paste0('~/Downloads/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.', i, '.l2.ldscore.gz'), data.table = F)
  }
  # dd = do.call(rbind, dd)
  # ff = do.call(rbind, ff)
  df_uk10k = do.call(rbind, ss)
  saveRDS(df_uk10k, baseline_lf_rds)
} else {
  df_uk10k = readRDS(baseline_lf_rds)
}

# load l2l4 pre-computed 
l2l4_mat_file = '~/Desktop/tmp/playground/ld4m/base_rsid_mafbin.csv'
# in matlab, do the following
# load('baselineLD.1kg.l2l4.mat', 'annot', 'RSIDs');
# tmp = table(annot(:, 1), RSIDs, sum(annot(:, 60:69), 2));
# writetable(tmp, 'base_rsid_mafbin.csv', 'Delimiter', ',');
df_l2l4 = data.table::fread(l2l4_mat_file, data.table = F, sep = ',') %>% mutate(SNP = paste0('rs', RSIDs)) %>% 
  rename(base = Var1, maf_bin = Var3)

# load 1000G SNPs
onekg_file = '~/Desktop/tmp/playground/ld4m/1000G_EUR.snp.rds'
if(!file.exists(onekg_file)) {
  ff = list()
  for(i in 1 : 22) {
    message(i)
    ff[[length(ff) + 1]] = data.table::fread(paste0('~/Downloads/1000G_Phase3_frq/1000G.EUR.QC.', i, '.frq'), data.table = F)
  }
  df_1kg = do.call(rbind, ff)
  saveRDS(df_1kg, onekg_file)
} else {
  df_1kg = readRDS(onekg_file)
}

# annotate l2l4 snps
df_l2l4 = df_l2l4 %>% left_join(df_1kg %>% select(SNP, MAF) %>% rename(MAF.1kg = MAF), by = 'SNP')
df_l2l4 = df_l2l4 %>% left_join(df_uk10k %>% rename(common.uk10k = common, lowfrq.uk10k = lowfrq), by = 'SNP')

common_snp_pool0 = df_l2l4$SNP
common_snp_pool1 = union(df_uk10k$SNP[df_uk10k$common == 1], common_snp_pool0)
common_snp_pool = union(df_1kg$SNP[df_1kg$MAF > 0.05], common_snp_pool1)
# common_snp_pool = union(common_snp_pool1, common_snp_pool2)
df_common_snp = data.frame(SNP = common_snp_pool)
df_common_snp$UK10K = df_common_snp$SNP %in% df_uk10k$SNP[df_uk10k$common == 1]
df_common_snp$l2l4_maf_bin = df_common_snp$SNP %in% df_l2l4$SNP[df_l2l4$maf_bin == 1]
df_common_snp$l2l4_base = df_common_snp$SNP %in% df_l2l4$SNP[df_l2l4$base == 1]
df_common_snp$onekg = df_common_snp$SNP %in% df_1kg$SNP[df_1kg$MAF > 0.05]
df_common_snp$l2l4 = df_common_snp$SNP %in% df_l2l4$SNP
colSums(df_common_snp %>% select(-SNP))

message(
  'Count of common SNPs \n',
  'UK10K = ', sum(df_uk10k$common), '\n',
  'l2l4 (maf_bin) = ', sum(df_l2l4$maf_bin), '\n',
  'l2l4 (base) = ', sum(df_l2l4$base), '\n',
  '1000G = ', sum(df_1kg$MAF > 0.05), '\n',
  'l2l4 = ', nrow(df_l2l4), '\n'
)

check_2_is_subset_of_1 = function(l1, l2) {
  sum(l2 > l1) == 0
}


message('l2l4 base - l2l4 maf_bin')
table(paste(df_common_snp$l2l4_base, df_common_snp$l2l4_maf_bin))
message(
  'l2l4 maf_bin is subset of l2l4 base: ', 
  check_2_is_subset_of_1(df_common_snp$l2l4_base, df_common_snp$l2l4_maf_bin)
)

message('UK10K - l2l4 maf_bin')
table(paste(df_common_snp$UK10K, df_common_snp$l2l4_maf_bin))
message(
  'l2l4 maf_bin is subset of UK10K: ', 
  check_2_is_subset_of_1(df_common_snp$UK10K, df_common_snp$l2l4_maf_bin)
)

message('UK10K - l2l4 base')
table(paste(df_common_snp$UK10K, df_common_snp$l2l4_base))
message(
  'l2l4 base is subset of UK10K: ', 
  check_2_is_subset_of_1(df_common_snp$UK10K, df_common_snp$l2l4_base)
)

message('1000G - l2l4')
table(paste(df_common_snp$onekg, df_common_snp$l2l4))
message(
  'l2l4 is subset of 1000G: ', 
  check_2_is_subset_of_1(df_common_snp$onekg, df_common_snp$l2l4)
)

message('UK10K - l2l4')
table(paste(df_common_snp$UK10K, df_common_snp$l2l4))
message(
  'l2l4 is subset of UK10K: ', 
  check_2_is_subset_of_1(df_common_snp$UK10K, df_common_snp$l2l4)
)

message('Conclusion: l2l4 maf_bin is a subset of UK10K common SNPs.')

not_l2l4_maf_bin = df_l2l4$SNP[df_l2l4$maf_bin != 1]
message(length(not_l2l4_maf_bin), ' SNPs in l2l4 non-maf_bin.')
message(sum(not_l2l4_maf_bin %in% df_uk10k$SNP[df_uk10k$common == 1]), ' SNPs in l2l4 non-maf_bin are common in UK10K.')
not_l2l4_maf_bin_common_uk10k = not_l2l4_maf_bin[not_l2l4_maf_bin %in% df_uk10k$SNP[df_uk10k$common == 1]]
message('Among these, ', sum(not_l2l4_maf_bin_common_uk10k %in% df_1kg$SNP[df_1kg$MAF > 0.05]), ' SNPs are common in 1000G.')
message(sum(not_l2l4_maf_bin %in% df_1kg$SNP[df_1kg$MAF > 0.05]), ' SNPs in l2l4 non-maf_bin are common in 1000G.')



message('Final conclusion: ')
message('1. l2l4 maf_bin is a subset of UK10K common SNPs. ')
message('2. l2l4 non-maf_bin is almost not common in UK10K. ')
message('3. |l2l4 maf_bin| = ', sum(df_common_snp$l2l4_maf_bin), '; |UK10K common| = ', sum(df_common_snp$UK10K))
message('4. Aggregate the results from the non-overlapping l2l4 maf_bin bins. It may be slightly smaller than the one from UK10K common since it contains fewer SNPs.')
