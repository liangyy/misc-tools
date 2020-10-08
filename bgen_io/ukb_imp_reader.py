import gc
import sqlite3
import pandas as pd
import numpy as np
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
pandas2ri.activate()

class UKBReader:
    
    def __init__(self, bgen, bgi, sample):
        self.rbgen = importr("rbgen")
        self.bgen_path = bgen
        self.bgi_path = bgi
        self.sample_path = sample
        self._init_rsid()
        # self._index_variant()
    
    def _init_rsid(self):
        with sqlite3.connect(self.bgi_path) as conn:
            variants = conn.execute('select * from Variant').fetchall()
        self.rsids = [ v[2] for v in variants ]


    # def _index_variant(self):
    #     if hasattr(self, 'variant_index'):
    #         return
    #     with sqlite3.connect(self.bgi_path) as conn:
    #         variants = conn.execute('select * from Variant').fetchall()
    #     self.variant_index = {}
    #     for i_ in range(len(variants)):
    #         i = variants[i_]
    #         varid = i[2]
    #         my_var_id = '{}:{}:{}:{}'.format(i[0], i[1], i[4], i[5])
    #         self.variant_index[my_var_id] = varid
    
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
        # dosage: nvariant (we have 1 here) x nsample x num_of_allele_combination
        dosage = dosage[0, :, 1] + 2 * dosage[0, :, 2]
        missing = np.isnan(dosage)
        dosage[missing] = np.nanmean(dosage)
        return dosage 
        
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
        dosage = pandas2ri.ri2py(cached_data[4])
        # dosage: nvariant (we have 1 here) x nsample x num_of_allele_combination
        dosage = dosage[:, :, 1] + 2 * dosage[:, :, 2]
        missing = np.where(np.isnan(dosage))[0]
        dosage[missing[0], missing[1]] = np.nanmean(dosage, axis=1)[missing[0]]
        return dosage 

    def dosage_generator(self, id_list):
        for snpid in id_list:
            yield self.get_dosage_by_id(snpid)
