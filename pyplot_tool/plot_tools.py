import matplotlib.pyplot as plt
import numpy as np

def rank_(vec):
    '''
    input: np.array 1-dimensional
    '''
    temp = vec.argsort()
    ranks = np.empty_like(temp)
    ranks[temp] = np.arange(len(vec))
    return ranks

def get_expected_p(obs_p):
    return (rank_(obs_p) + 1) / (obs_p.shape[0] + 1) 

def _qqplot(axe, exp, obs, label):
    axe.scatter(-np.log10(exp), -np.log10(obs), label=label)
    return -np.log10(exp.min())

def qqplot(lpval, labels, axe, pval_cutoff=1):
    '''
    input: list of np.array 1-dimensional
    '''
    max_range = 0
    for pval, ll in zip(lpval, labels):
        pexp = (rank_(pval) + 1) / (pval.shape[0] + 2)   
        xmax = _qqplot(axe, pexp * pval_cutoff, pval, ll)
        # axe.scatter(-np.log10(pexp), -np.log10(pval), label=ll)
        max_range = max(max_range, xmax)
    axe.plot([0, max_range], [0, max_range], c='black')
    axe.legend()

def qqplot_quick(lpval, labels, axe, pval_cutoff=1, reference=None):
    if reference is None:
        reference = np.concatenate([
            np.arange(0, 0.005, step=0.01/10000), 
            np.arange(0.005, 0.01, step=0.0001), 
            np.arange(0.01, 1, step=0.001)
        ])
        reference = reference[1:]
    max_range = 0
    for pval, ll in zip(lpval, labels):    
        observe = np.quantile(pval, reference)
        xmax = _qqplot(axe, reference * pval_cutoff, observe, ll)
        # axe.scatter(-np.log10(reference * pval_cutoff), -np.log10(observe), label=ll)
        max_range = max(max_range, xmax) 
    axe.plot([0, max_range], [0, max_range], c='black')
    axe.legend()

def scatter(lx, ly, labels, axe, kwargs, diag_line=True):
    '''
    input lx and ly: list of np.array 1-dimensional
    '''
    max_ = -np.float('inf')
    min_ = np.float('inf')
    for x, y, l in zip(lx, ly, labels):
        axe.scatter(x, y, label=l, **kwargs)
        if diag_line is True:
                max_ = max(np.max(x), np.max(y), max_)
                min_ = min(np.min(x), np.min(y), min_)
#     print(max_, min_)    
    if diag_line is True:
        axe.plot([min_, max_], [min_, max_], color='black', linestyle='dashed')
    axe.legend()
