import gc
import sqlite3
import pandas as pd
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

class UKBReader:
    
    def __init__(self, bgen, bgi, no_var=True):
        self.rbgen = importr("rbgen")
        self.bgen_path = bgen
        self.bgi_path = bgi
        if no_var is False:
            self._init_rsid()
    
    def _init_rsid(self):
        with sqlite3.connect(self.bgi_path) as conn:
            variants = conn.execute('select * from Variant').fetchall()
        self.rsids = set([ v[2] for v in variants ])
        # self.a0 = [ v[4] for v in variants ]
        # self.a1 = [ v[5] for v in variants ]
    
    def check_rsid(self, rsid_list):
        '''
        Check if each element of rsid_list is in the bgen file.
        Return the list containing the elements that are in the bgen only. 
        '''
        if hasattr(self, 'rsids'):
            o = []
            for rsid in rsid_list:
                if rsid in self.rsids:
                    o.append(rsid)
            return o
        else:
            self._init_rsid()
            self.check_rsid(rsid_list)
        
    
    def get_dosage_by_id(self, snpid):
        '''
        Extract genotype dosage prob
        and convert to expected dosage of second allele by
        dosage[:, 1] + 2 * dosage[:, 2].
        The samples labelled as missing will be imputed to population mean.
        '''
        cached_data = self.rbgen.bgen_load(
            self.bgen_path,
            index_filename=self.bgi_path,
            rsids=snpid, 
            max_entries_per_sample=4
        )
        dosage = pandas2ri.ri2py(cached_data[4])
        df_var = pandas2ri.ri2py(cached_data[0])
        # dosage: nvariant (we have 1 here) x nsample x num_of_allele_combination
        with np.errstate(invalid='ignore'):
            dosage = dosage[0, :, 1] + 2 * dosage[0, :, 2]
        missing = np.isnan(dosage)
        dosage[missing] = np.nanmean(dosage)
        return dosage, df_var.allele0[0], df_var.allele1[0]
        
    def get_dosage_by_chunk(self, snpid_list):
        '''
        Extract genotype dosage prob
        and convert to expected dosage of second allele by
        dosage[:, 1] + 2 * dosage[:, 2].
        The samples labelled as missing will be imputed to population mean.
        '''
        cached_data = self.rbgen.bgen_load(
            self.bgen_path,
            index_filename=self.bgi_path,
            rsids=StrVector(snpid_list), 
            max_entries_per_sample=4
        )
        with localconverter(ro.default_converter + pandas2ri.converter):
            df_var = ro.conversion.rpy2py(cached_data[0])
            dosage = ro.conversion.rpy2py(cached_data[4])
        # dosage: nvariant (we have several here) x nsample x num_of_allele_combination
        with np.errstate(invalid='ignore'):
            dosage = dosage[:, :, 1] + 2 * dosage[:, :, 2]
        missing_ind = np.isnan(dosage)
        if missing_ind.sum() > 0:
            missing = np.where(missing_ind)
            dosage[missing[0], missing[1]] = np.nanmean(dosage, axis=1)[missing[0]]
        # need to handle potential duplicated snpid (due to multialleleic)
        _, dup_idx = np.unique(df_var.rsid.tolist(), return_index=True)
        dosage = dosage[dup_idx, :]
        df_var = df_var.iloc[ dup_idx, : ]
        # need to reorder so that it is the same as input query
        match_idx = self.match_y_to_x(df_var.rsid.values, np.array(snpid_list))
        df_var = df_var.iloc[ match_idx, : ]
        return dosage[match_idx, :], df_var.allele0.tolist(), df_var.allele1.tolist(), df_var.rsid.tolist()
    
    @staticmethod
    def match_y_to_x(x, y):
        '''
        x, y are 1d np.array 
        y is a subset of x.
        return idx such that x[idx] = y
        '''
        return np.where(y.reshape(y.size, 1) == x)[1]

    def dosage_generator(self, id_list):
        for snpid in id_list:
            yield self.get_dosage_by_id(snpid)
    
    @staticmethod
    def _split_list_into_chunks(mylist, chunk_size):
        ntotal = len(mylist)
        nchunk = ntotal // chunk_size 
        if ntotal > chunk_size * nchunk:
            nchunk += 1
        curr_chunk = []
        curr_idx = []
        chunk_list = []
        idx_list = []
        for idx, i in enumerate(mylist):
            curr_chunk.append(i)
            curr_idx.append(idx)
            if len(curr_chunk) == chunk_size:
                chunk_list.append(curr_chunk)
                idx_list.append(curr_idx)
                curr_chunk = []
                curr_idx = []
        if len(curr_chunk) != 0:
            chunk_list.append(curr_chunk)
            idx_list.append(curr_idx)
        return chunk_list, idx_list
    
    def get_nchunk(self, id_list, chunk_size=20):
        chunk_size = min(chunk_size, len(id_list)) 
        chunk_list, idx_list = self._split_list_into_chunks(id_list, chunk_size)
        return len(chunk_list)
    
    def dosage_generator_by_chunk(self, id_list, chunk_size=20):
        chunk_size = min(chunk_size, len(id_list)) 
        chunk_list, idx_list = self._split_list_into_chunks(id_list, chunk_size)
        
        for snp_chunk, idx_chunk in zip(chunk_list, idx_list):
            yield self.get_dosage_by_chunk(snp_chunk), snp_chunk, idx_chunk


