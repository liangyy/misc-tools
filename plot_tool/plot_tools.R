qqplot_by_group <- function(pval, group, pval_cutoff = 1, ...) {
  n <- length(pval)
  pexp <- rank(pval) / n
  df <- data.frame(p.val = pval, grp = group) %>% group_by(grp) %>% mutate(p.exp = pval_cutoff * rank(p.val) / (n() + 1)) %>% ungroup()
  p <- ggplot(df) + 
    geom_point(aes(x = -log10(p.exp), y = -log10(p.val), color = grp), ...) + 
    geom_hline(yintercept = -log10(0.05 / n)) + 
    geom_abline(slope = 1, intercept = 0, linetype = 2)
  return(p)
}

qqplot_quick_by_group <- function(pval, group, pval_cutoff = 1, reference = NULL, ...) {
  if(is.null(reference)) {
    reference = c(
      seq(0, 0.005, by = 0.01 / 10000),
      seq(0.005, 0.01, by = 0.0001),
      seq(0.01, 1, by = 0.001)
    )
    reference = reference[ reference >= 1 / length(pval) ]
  }
  n <- length(pval)
  observe = as.numeric(quantile(pval, probs = reference, type = 1))
  dup_ind = duplicated(observe)
  observe = observe[!dup_ind]
  reference = reference[!dup_ind]
  df <- data.frame(p.val = observe, grp = group) %>% group_by(grp) %>% mutate(p.exp = reference * pval_cutoff) %>% ungroup()
  p <- ggplot(df) + 
    geom_point(aes(x = -log10(p.exp), y = -log10(p.val), color = grp), ...) + 
    geom_hline(yintercept = -log10(0.05 / n)) + 
    geom_abline(slope = 1, intercept = 0, linetype = 2)
  return(p)
}