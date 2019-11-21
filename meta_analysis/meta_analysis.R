sample_size_based_meta = function(z1, z2, n1, n2) {
  w1 = sqrt(n1)
  w2 = sqrt(n2)
  z_meta = (z1 * w1 + z2 * w2) / sqrt(w1 ^ 2 + w2 ^ 2)
  data.frame(z_meta = z_meta)
}

inverse_variance_based_meta = function(b1, b2, se1, se2) {
  w1 = 1 / (se1 ^ 2)
  w2 = 1 / (se2 ^ 2)
  se_meta = sqrt(1 / (w1 + w2))
  b_meta = (b1 * w1 + b2 * w2) / (w1 + w2)
  o = data.frame(effect_size_meta = b_meta, se_meta = se_meta)
  if_na = is.na(b1) | is.na(se1)
  o$effect_size_meta[if_na] = b2[if_na]
  o$se_meta[if_na] = se2[if_na]
  if_na = is.na(b2) | is.na(se2)
  o$effect_size_meta[if_na] = b1[if_na]
  o$se_meta[if_na] = se1[if_na]
  o
}

zval2pval = function(z) {
   exp(pnorm(abs(z), log.p=T, low=F)) * 2
}

zval2beta = function(af, n, z) {
  z / sqrt(2 * n * af * (1 - af))
}