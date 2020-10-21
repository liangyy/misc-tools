set.seed(1)

n = 1000
k = 1500
Sigma = matrix(rnorm(n * k), nrow = n)
Sigma = apply(Sigma, 2, function(x) { (x - mean(x)) / sd(x) })
Sigma = Sigma %*% t(Sigma) / k
# image(Sigma)
Vr = 5
Ve = 8
evd = eigen(Sigma)
y = as.numeric(mvtnorm::rmvnorm(1, sigma = Vr * Sigma + Ve * diag(n)))

ytil = t(evd$vectors) %*% y

obj = function(params, ytil, lambda) {
  u = params[1]
  v = params[2]
  ve = exp(u)
  vr = exp(v)
  sigma = ve + vr * lambda
  - 0.5 * sum(log(sigma)) - 0.5 * sum(ytil ^ 2 / sigma)
}

grad = function(params, ytil, lambda) {
  u = params[1]
  v = params[2]
  ve = exp(u)
  vr = exp(v)
  denom = ve + vr * lambda
  grad_ve = - 0.5 * sum(1 / denom) + 0.5 * sum(ytil ^ 2 / denom ^ 2)
  grad_vr = - 0.5 * sum(lambda / denom) + 0.5 * sum(ytil ^ 2 * lambda / denom ^ 2)
  c(grad_ve * ve, grad_vr * vr)
}

jacob = function(params, ytil, lambda) {
  u = params[1]
  v = params[2]
  ve = exp(u)
  vr = exp(v)
  denom = ve + vr * lambda
  grad2_ve = - 0.5 * sum(1 / denom ^ 2) - sum(ytil ^ 2 / denom ^ 3)
  grad2_vr = - 0.5 * sum(lambda ^ 2 / denom ^ 2) - sum(ytil ^ 2 * lambda ^ 2 / denom ^ 3)
  grad2_vevr = - 0.5 * sum(lambda / denom ^ 2) - sum(ytil ^ 2 * lambda / denom ^ 3)
  grad_ve = - 0.5 * sum(1 / denom) + 0.5 * sum(ytil ^ 2 / denom ^ 2)
  grad_vr = - 0.5 * sum(lambda / denom) + 0.5 * sum(ytil ^ 2 * lambda / denom ^ 2)
  a11 = (grad2_ve * ve + grad_ve) * ve
  a12 = grad2_vevr * ve * vr
  a22 = (grad2_vr * vr + grad_vr) * vr
  matrix(c(a11, a12, a12, a22), byrow = T, nrow = 2)
}

newton_raphson = function(obj, jacob, x_init, tol = 1e-10, ...) {
  x_curr = x_init
  diff = Inf
  while(diff > tol) {
    x = x_curr - solve(jacob(x_curr, ...), obj(x_curr, ...))
    diff = sum((x - x_curr) ^ 2)
    x_curr = x
  }
  x_curr
}

mixed_effect_solver = function(y, Sigma) {
  evd = eigen(Sigma)
  ytil = t(evd$vectors) %*% y
  soln_curr = optim(c(0, 0), fn = obj, gr = grad, control = list(fnscale = -1), lambda = evd$values, ytil = ytil)
  soln_refined = exp(newton_raphson(grad, jacob, soln_curr$par, ytil = ytil, lambda = evd$values))
  ll = obj(log(soln_refined), ytil, evd$values)
  sigma02 = sum(y ^ 2) / (length(y) - 1)
  ll0 = -0.5 * length(y) * log(sigma02) - 0.5 * sum(y ^ 2) / sigma02
  return(list(lr = ll - ll0, ve = soln_refined[1], vg = soln_refined[2]))
}

oo = mixed_effect_solver(y, Sigma)
oo

oo2 = emma.MLE.noX(y, Sigma)
oo2

