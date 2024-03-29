import numpy as np
import scipy.optimize
import pandas as pd

import gzip

'''
Note that here the obj and corresponding grad and hessian 
are taking a negative sign (relative to emma_solver.R)
so that we do minimization.
'''

def load_grm_id(grm_id):
    o = []
    with open(grm_id, 'r') as f:
        for i in f:
            i = i.strip()
            o.append(i.split('\t')[1])
    return o

def load_grm(grm_txt_gz, grm_id):
    grm_indiv = load_grm_id(grm_id)
    nindiv = len(grm_indiv)
    grm_mat = np.zeros((nindiv, nindiv))
    with gzip.open(grm_txt_gz, 'rt') as f:
        for i in f:
            i = i.strip()
            x, y, _, val = i.split('\t')
            x, y, val = int(x), int(y), float(val)
            grm_mat[x - 1, y - 1] = val
            grm_mat[y - 1, x - 1] = val
    return grm_mat, grm_indiv

def standardize(x):
    return (x - x.mean()) / x.std()

def get_ll0(y, X=None):
    if X is not None:
        q, _ = np.linalg.qr(X)
        r = y - q @ (q.T @ y)
    else:
        r = y 
    sigma02 = (r ** 2).sum() / (r.shape[0] - 1)
    ll0 = -0.5 * r.shape[0] * np.log(sigma02) - 0.5 * (r ** 2).sum() / sigma02
    return ll0

def _calc_grad(ve, vr, denom, ytil, lambda_):
    grad_ve = - 0.5 * (1 / denom).sum() + 0.5 * (ytil ** 2 / denom ** 2).sum()
    grad_vr = - 0.5 * (lambda_ / denom).sum() + 0.5 * (ytil ** 2 * lambda_ / denom ** 2).sum()
    return grad_ve, grad_vr

def _calc_grad2(ve, vr, denom, ytil, lambda_):
    grad2_ve = + 0.5 * (1 / (denom ** 2)).sum() - (ytil ** 2 / denom ** 3).sum()
    grad2_vr = + 0.5 * ((lambda_ ** 2) / (denom ** 2)).sum() - ((ytil ** 2) * (lambda_ ** 2) / denom ** 3).sum()
    grad2_vevr = + 0.5 * (lambda_ / (denom ** 2)).sum() - ((ytil ** 2) * lambda_ / (denom ** 3)).sum()
    return grad2_ve, grad2_vr, grad2_vevr

def obj(params, ytil, lambda_):
    ve, vr = np.exp(params)
    sigma = ve + vr * lambda_
    res = - 0.5 * np.log(sigma).sum() - 0.5 * (ytil ** 2 / sigma).sum()
    return - res

def grad(params, ytil, lambda_):
    ve, vr = np.exp(params)
    denom = ve + vr * lambda_
    grad_ve, grad_vr = _calc_grad(ve, vr, denom, ytil, lambda_)
    return - np.array([grad_ve * ve, grad_vr * vr])

def jacob(params, ytil, lambda_):
    ve, vr = np.exp(params)
    denom = ve + vr * lambda_
    grad_ve, grad_vr = _calc_grad(ve, vr, denom, ytil, lambda_)
    grad2_ve, grad2_vr, grad2_vevr = _calc_grad2(ve, vr, denom, ytil, lambda_)
    a11 = (grad2_ve * ve + grad_ve) * ve
    a12 = grad2_vevr * ve * vr
    a22 = (grad2_vr * vr + grad_vr) * vr
    return - np.array([[a11, a12], [a12, a22]])

def grad2_beta_vevr(W, X, resid, lambda_):
    g2_ve = X.T @ (- resid * (W ** 2))
    g2_vr = X.T @ (- resid * lambda_ * (W ** 2))
    return g2_ve, g2_vr

def fisher_information_full(params, resid, W, Xtil, lambda_):
    nbeta = Xtil.shape[1]
    hessian = np.zeros((nbeta + 2, nbeta + 2))
    hessian[:2, :][:, :2] = - fisher_information(params, resid, lambda_)
    hessian[2:, :][:, 2:] = calc_XtWX(Xtil, W)
    grad2_b_ve, grad2_b_vr = grad2_beta_vevr(W, Xtil, resid, lambda_)
    hessian[0, :][2:] = hessian[:, 0][2:] = grad2_b_ve
    hessian[1, :][2:] = hessian[:, 1][2:] = grad2_b_vr
    return - hessian

def calc_XtWX(X, W):
    return np.einsum('ij,jk,j->ik', X.T, X, W)

def fisher_information(params, ytil, lambda_):
    '''
    in the scale of np.exp(params)
    '''
    ve, vr = np.exp(params)
    denom = ve + vr * lambda_
    grad2_ve, grad2_vr, grad2_vevr = _calc_grad2(ve, vr, denom, ytil, lambda_)
    return - np.array([[grad2_ve, grad2_vevr], [grad2_vevr, grad2_vr]])

def newton_raphson(obj_, jacob_, x_init, args, tol=1e-10):
    x_curr = x_init
    diff = float('inf')
    try:
        while diff > tol:
            x = x_curr - np.linalg.solve(
                jacob_(x_curr, *args), 
                obj_(x_curr, *args)
            ) 
            diff = ((x - x_curr) ** 2).sum()
            x_curr = x
        return x_curr
    except np.linalg.LinAlgError:
        return x_curr

def delta_method(mu, varcov, h_fun, h_grad):
    '''
    Return the estimate as SE of h_fun(x) given
    x ~ N(mu, varcov)
    '''
    tmp = h_grad(mu)
    h_mu = h_fun(mu)
    h_var = tmp.T @ varcov @ tmp
    return h_mu, h_var

def vp_fun(est):
    '''
    est = (ve, vr)
    '''
    return est[0] + est[1]

def vp_grad(est):
    return np.ones((2))

def ratio_fun(est):
    '''
    est = (ve, vr)
    '''
    return est[1] / (est[0] + est[1])

def ratio_grad(est):
    return np.array([- est[1] / (est[0] + est[1]) ** 2, est[0] / (est[0] + est[1]) ** 2])

def calc_vp(est, varcov):
    '''
    est = (ve, vr)
    '''
    return delta_method(est, varcov, vp_fun, vp_grad)

def calc_h2(est, varcov):
    return delta_method(est, varcov, ratio_fun, ratio_grad)

def null_df():
    return pd.DataFrame({'LR': np.nan, 'Vg': np.nan, 'Ve': np.nan, 'Vp': np.nan, 'Vg_SE': np.nan, 'Ve_SE': np.nan, 'Vp_SE': np.nan, 'h2': np.nan, 'h2_SE': np.nan, 'L0': np.nan, 'L1': np.nan}, index=[0])

def pyemma_reml_mat_fac(X, grm, min_max_ratio=None):
    Q, _ = np.linalg.qr(X)
    Pgrm = grm - Q @ (Q.T @ grm)
    PgrmP = Pgrm.T - Q @ (Q.T @ Pgrm.T)
    q = X.shape[1]
    val, vec = pyemma_mle_mat_fac(PgrmP, min_max_ratio=min_max_ratio)
    return val[q:], vec[:, q:]

def pyemma_mle_mat_fac(grm, min_max_ratio=None):
    if min_max_ratio is None:
        return np.linalg.eigh(grm)
    else:
        val, vec = np.linalg.eigh(grm)
        to_remove_idx = np.where((val / val[-1]) <= min_max_ratio)[0]
        val[to_remove_idx] = 0
        return val, vec

def pyemma_reml(y, grm_eig_vec, grm_eig_val, ori_grm_eig_vec, ori_grm_eig_val, X):
    res = pyemma_no_X(y, grm_eig_vec, grm_eig_val)
    # print(res)
    ve, vr = res['Ve'][0], res['Vg'][0]
    ytil = ori_grm_eig_vec.T @ y
    Xtil = ori_grm_eig_vec.T @ X
    
    W = 1 / (ve + vr * ori_grm_eig_val)
    beta = np.linalg.solve(calc_XtWX(Xtil, W), Xtil.T @ (W * ytil))
    resid = ytil - Xtil @ beta
    est = np.array([ve, vr])
    soln = np.log(est)
    ll = - obj(soln, resid, ori_grm_eig_val)
    ll0 = get_ll0(y, X=X)
    varcov = fisher_information_full(soln, resid, W, Xtil, ori_grm_eig_val) 
    varcov = np.linalg.inv(varcov[:2, :2])
    # print(varcov)
    var_ve = varcov[0, 0] 
    var_vr = varcov[1, 1] 
    
    vp, vp_var = calc_vp(est, varcov[:2, :2])
    h2, h2_var = calc_h2(est, varcov[:2, :2])
    res = pd.DataFrame(
        {
            'LR': ll - ll0,
            'Vg': est[1],
            'Ve': est[0],
            'Vp': vp,
            'Vg_SE': np.sqrt(var_vr),
            'Ve_SE': np.sqrt(var_ve),
            'Vp_SE': np.sqrt(vp_var),
            'h2': h2,
            'h2_SE': np.sqrt(h2_var),
            'L0': ll0,
            'L1': ll # - 0.5 * ytil.shape[0] * np.log(2 * np.pi)
        }, 
        index=[0]
    )
    return res

def initialize_with_grid(y, ytil, lambda_, Xtil=None, ngrid=20):
    scale = y.var()
    scale_seq = np.array([scale * 0.5, scale, scale * 2])
    ratio_seq = np.arange(ngrid) / ngrid
    ratio_seq = ratio_seq[1:] # remove the leading zero
    lld = np.zeros((3, ratio_seq.shape[0]))
    for si, scale in enumerate(scale_seq):
        for ri, ratio in enumerate(ratio_seq):
            ve = ratio * scale
            vr = (1 - ratio) * scale
            soln = np.log([ve, vr])
            if Xtil is not None:
                W = 1 / (ve + vr * lambda_)
                beta = np.linalg.solve(calc_XtWX(Xtil, W), Xtil.T @ (W * ytil))
                resid = ytil - Xtil @ beta
                lld[si][ri] = - obj(soln, resid, lambda_)
            else:
                lld[si][ri] = - obj(soln, ytil, lambda_)
    opt_idx = np.unravel_index(np.nanargmax(lld, axis=None), lld.shape)
    ratio_opt = ratio_seq[opt_idx[1]]
    scale_opt = scale_seq[opt_idx[0]]
    ve = ratio_opt * scale_opt
    vr = (1 - ratio_opt) * scale_opt
    return np.log([ve, vr])    

def pyemma_w_X(y, X, grm_eig_vec, grm_eig_val, tol=1e-10):
    ytil = grm_eig_vec.T @ y
    Xtil = grm_eig_vec.T @ X
    beta = np.zeros((X.shape[1]))
    Q, _ = np.linalg.qr(X)
    soln_v = initialize_with_grid(y - Q @ (Q.T @ y), ytil, grm_eig_val, Xtil=Xtil)
    diff = tol + 1
    while diff > tol:
        resid = ytil - Xtil @ beta
        soln_curr = scipy.optimize.minimize(
            fun=obj, 
            x0=soln_v, 
            method='Nelder-Mead',
            args=(resid, grm_eig_val),
            # jac=grad
        )
        if soln_curr.success is False:
            return null_df() 
        soln_refined = newton_raphson(
            grad, 
            jacob, 
            soln_curr.x, 
            args=(resid, grm_eig_val)
        )
        
        ve, vr = np.exp(soln_refined)
        W = 1 / (ve + vr * grm_eig_val)
        beta_new = np.linalg.solve(calc_XtWX(Xtil, W), Xtil.T @ (W * ytil))
        
        diff = ((beta - beta_new) ** 2).sum() + ((soln_v - soln_refined) ** 2).sum()
        
        soln_v = soln_refined
        beta = beta_new
    
    resid = ytil - Xtil @ beta
    ve, vr = np.exp(soln_v)
    W = 1 / (ve + vr * grm_eig_val)
    ll = - obj(soln_refined, resid, grm_eig_val)
    ll0 = get_ll0(y, X=X)
    
    est = np.exp(soln_v)
    varcov = fisher_information_full(soln_v, resid, W, Xtil, grm_eig_val) 
    varcov = np.linalg.inv(varcov)
    
    var_ve = varcov[0, 0] 
    var_vr = varcov[1, 1] 
    
    vp, vp_var = calc_vp(est, varcov[:2, :2])
    h2, h2_var = calc_h2(est, varcov[:2, :2])
    res = pd.DataFrame(
        {
            'LR': ll - ll0,
            'Vg': est[1],
            'Ve': est[0],
            'Vp': vp,
            'Vg_SE': np.sqrt(var_vr),
            'Ve_SE': np.sqrt(var_ve),
            'Vp_SE': np.sqrt(vp_var),
            'h2': h2,
            'h2_SE': np.sqrt(h2_var),
            'L0': ll0,
            'L1': ll # - 0.5 * ytil.shape[0] * np.log(2 * np.pi)
        }, 
        index=[0]
    )
    return res

def pyemma_no_X(y, grm_eig_vec, grm_eig_val):
    ytil = grm_eig_vec.T @ y
    soln_init = initialize_with_grid(y, ytil, grm_eig_val)
    soln_curr = scipy.optimize.minimize(
        fun=obj, 
        method='Nelder-Mead',
        x0=soln_init,   # np.array([-0.5, -0.5]), # np.zeros(2), 
        args=(ytil, grm_eig_val),
        # jac=grad
    )
    if soln_curr.success is False:
        return null_df()
    soln_refined = newton_raphson(
        grad, 
        jacob, 
        soln_curr.x, 
        args=(ytil, grm_eig_val)
    )
    
    ll = - obj(soln_refined, ytil, grm_eig_val)
    ll0 = get_ll0(y)
    
    est = np.exp(soln_refined)
    varcov = fisher_information(soln_refined, ytil, grm_eig_val) 
    varcov = np.linalg.inv(varcov)
    # print(varcov)
    var_ve = varcov[0, 0] 
    var_vr = varcov[1, 1] 
    
    vp, vp_var = calc_vp(est, varcov)
    h2, h2_var = calc_h2(est, varcov)
    res = pd.DataFrame(
        {
            'LR': ll - ll0,
            'Vg': est[1],
            'Ve': est[0],
            'Vp': vp,
            'Vg_SE': np.sqrt(var_vr),
            'Ve_SE': np.sqrt(var_ve),
            'Vp_SE': np.sqrt(vp_var),
            'h2': h2,
            'h2_SE': np.sqrt(h2_var),
            'L0': ll0,
            'L1': ll # - 0.5 * ytil.shape[0] * np.log(2 * np.pi)
        }, 
        index=[0]
    )
    return res
    
