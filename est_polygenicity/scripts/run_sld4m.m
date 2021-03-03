function [] = run_sld4m(l2l4_file, gwas_file, out_file, mode_)
% INPUT:
% l2l4_file: MAT file
% gwas_file: MAT file
% out_file: csv format
% mode: annot = 1 (mode = simple) or annot = all (mode = all)

load(l2l4_file);
load(gwas_file, 'chisq', 'sumstat_RSIDs_as_ints');

Annot_names = string(header)';
no_blocks=100;
[~, idx_ref, idx_sumstat] = intersect(RSIDs, sumstat_RSIDs_as_ints, 'stable');
fprintf('number of SNPs to use = %d \n', size(idx_ref, 1));

if strcmp(mode_, 'simple')
  ell2 = ell2(:, 1);
  ell4 = ell4(:, 1);
  annot = annot(:, 1);
  Annot_names = Annot_names(1);
  report_annot_mat = [ 1 1 ];
elseif strcmp(mode_, 'all')
  % looks good
  tmp = zeros(size(ell4, 2), 1);
  tmp(60:69) = 1;
  disp(['Aggregate:', header(60:69)])
  report_annot_mat = [ tmp, eye(size(ell4, 2)) ];
else
  fprintf('Wrong mode = %s', mode_); 
end
fprintf('Mode = %s', mode_);

[Ma_est, Ma_err] = SLD4M( chisq(idx_sumstat), ell2, ell4, idx_ref, annot, EW2(idx_ref), EW4(idx_ref), no_blocks, 1, report_annot_mat);

out_table = table([ "Manual_aggregated", Annot_names' ]', Ma_est, Ma_err);
writetable(out_table, out_file, 'Delimiter', ',');

exit

