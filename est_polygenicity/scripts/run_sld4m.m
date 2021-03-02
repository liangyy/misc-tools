function [] = run_sld4m(l2l4_file, gwas_file, out_file, mode)
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


if mode == 'simple'
  ell2 = ell2(:, 1)
  ell4 = ell4(:, 1)
  annot = annot(:, 1)
  Annot_names = Annot_names(1)
end

[Ma_est, Ma_err] = SLD4M( chisq(idx_sumstat), ell2, ell4, idx_ref, annot, EW2(idx_ref), EW4(idx_ref), no_blocks);

out_table = table(Annot_names, Ma_est, Ma_err);
writetable(out_table, out_file, 'Delimiter', ',');
