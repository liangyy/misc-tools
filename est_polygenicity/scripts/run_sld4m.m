function [] = run_sld4m(l2l4_file, gwas_file, out_file, mode_)
% INPUT:
% l2l4_file: MAT file
% gwas_file: MAT file
% out_file: csv format
% mode: annot = base (mode = base)
%       annot = all annotations (mode = all)
%       annot = all annotations with one extra that aggregates common MAF bins 
%               (mode = all_w_agg)

load(l2l4_file);
load(gwas_file, 'chisq', 'sumstat_RSIDs_as_ints');

Annot_names = string(header)';
no_blocks=100;
[~, idx_ref, idx_sumstat] = intersect(RSIDs, sumstat_RSIDs_as_ints, 'stable');
fprintf('number of SNPs to use = %d \n', size(idx_ref, 1));
fprintf('running mode = %s \n', mode_);

if strcmp(mode_, 'base')
  ell2 = ell2(:, 1);
  ell4 = ell4(:, 1);
  annot = annot(:, 1);
  Annot_names = Annot_names(1);
  report_annot_mat = [ 1 ];
  ref_idx = 1;
elseif strcmp(mode_, 'all_w_agg')
  % add aggregation as the first category
  tmp = zeros(size(ell4, 2), 1);
  tmp(60:69) = 1;
  disp(['Aggregate:', header(60:69)])
  report_annot_mat = [ tmp, eye(size(ell4, 2)) ];
  ref_idx = 2;
elseif strcmp(mode_, 'all')
  report_annot_mat = eye(size(ell4, 2));
  ref_idx = 1;
else
  fprintf('Wrong mode = %s', mode_); 
end
fprintf('Mode = %s', mode_);

[ Ma_est, Ma_err, ...
  h2_est, h2_err, ...
  Maenrich_est, Maenrich_err, ...
  h2enrich_est, h2enrich_err, ...
  ma_jk, h2_jk, ...
  kurtexp_est, kurtexp_err, ...
  propkurt_est, propkurt_err ] = SLD4M( chisq(idx_sumstat), ell2, ell4, idx_ref, annot, EW2(idx_ref), EW4(idx_ref), no_blocks, ref_idx, report_annot_mat);

if strcmp(mode_, 'base') | strcmp(mode_, 'all')
  outnames = Annot_names;
elseif strcmp(mode_, 'all_w_agg')
  outnames = [ "Manual_aggregated", Annot_names' ]';
end 

out_table = table(outnames, ...
  Ma_est, Ma_err, ...
  h2_est, h2_err, ...
  Maenrich_est, Maenrich_err, ...
  h2enrich_est, h2enrich_err);
writetable(out_table, out_file, 'Delimiter', ',');

exit

