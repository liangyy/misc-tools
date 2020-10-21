import numpy as np
import pyemma
from imp import reload
pyemma = reload(pyemma)
from rpy2.robjects.packages import importr, data
from rpy2.robjects import numpy2ri
numpy2ri.activate()
emma = importr('emma')

n = 1000
k = 1500
Sigma = np.random.normal(size=(n, k))
Sigma = np.apply_along_axis(pyemma.standardize, axis=0, arr=Sigma)
Sigma = Sigma @ Sigma.T / k
Vr = 5
Ve = 8
val, vec = pyemma.pyemma_mle_mat_fac(Sigma)
val_r, vec_r = pyemma.pyemma_reml_mat_fac(np.ones((n, 1)), Sigma)

y = np.random.normal(size=(n))
tmp_val, tmp_vec = np.linalg.eigh(Vr * Sigma + Ve * np.eye(n))
y = tmp_vec @ np.diag(np.sqrt(tmp_val)) @ y
res = emma.emma_MLE_noX(y, Sigma)
print(res.rx2('ve'), res.rx2('vg'), res.rx2('ML'))
res = emma.emma_MLE(y, np.ones((n, 1)), Sigma)
print(res.rx2('ve'), res.rx2('vg'), res.rx2('ML'))
res = emma.emma_REMLE(y, np.ones((n, 1)), Sigma)
print(res.rx2('ve'), res.rx2('vg'), res.rx2('REML'))
print('pyemma_no_X')
print(pyemma.pyemma_no_X(y, vec, val))
print('pyemma_w_X')
print(pyemma.pyemma_w_X(y, np.ones((n, 1)), vec, val))
print('pyemma_reml')
print(pyemma.pyemma_reml(y, vec_r, val_r, vec, val, np.ones((n, 1))))