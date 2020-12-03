import numpy as np
import scipy.stats

def return_constant_ind(mat):
    std = np.std(mat, axis=0)
    return std == 0

def inv_norm_by_col(mat):
    res_mat, good_ind = _prep_mat(mat)
    res_mat[:, good_ind] = np.apply_along_axis(inv_norm_vec, 0, mat[:, good_ind])
    return res_mat
    
def inv_norm_vec(vec, offset = 1):
    rank = myrank(vec)
    return scipy.stats.norm.ppf(rank / (len(rank) + offset), loc=0, scale=1)

def standardize_by_col(mat):
    res_mat, good_ind = _prep_mat(mat)
    res_mat[:, good_ind] = np.apply_along_axis(standardize_vec, 0, mat[:, good_ind])
    return res_mat

def _prep_mat(mat):
    good_ind = np.logical_not(return_constant_ind(mat))
    res_mat = np.zeros(mat.shape)
    return res_mat, good_ind

def standardize_vec(vec):
    return _divide(vec - np.mean(vec), np.std(vec))

def _divide(a, b):
    return np.divide(a, b, out=np.zeros_like(a), where=(b != 0))

def myrank(vec):
    argsort = np.argsort(vec)
    ranks = np.empty_like(argsort)
    ranks[argsort] = np.arange(len(vec))
    return ranks + 1  # rank starts from 1
    
def z2p(zscore):
    '''
    Input 1d np.array zscore and return the corresponding two-sided p-value.
    '''
    return scipy.stats.norm.sf(np.abs(zscore)) * 2
    
