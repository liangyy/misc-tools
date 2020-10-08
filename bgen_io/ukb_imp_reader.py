import bgen_reader
import numpy as np
class UKBReader:
    
    def __init__(self, bgen, metafile, sample, copy_from=None, byid='rsid'):
        '''
        Initialize self.id_dict as a dictionary with 
          key being byid = rsid or id; 
          value being the index in genotype file.
        '''
        self.bgen = bgen_reader.open_bgen(bgen, metafile_filepath=metafile, samples_filepath=sample)
        self.byid = byid

        if copy_from is None:
            self._init_id_dict()
        else:
            self.id_dict = copy_from.id_dict.copy()
        
    def _init_id_dict(self):
        if self.byid == 'rsid':
            self.id_dict = { self.rsids[i]: i for i in range(len(self.rsids)) }
        elif self.byid == 'id':
            self.id_dict = { self.ids[i]: i for i in range(len(self.ids)) }
        else:
            raise ValueError('Only support rsid and id.')
    
    def get_dosage_by_id(self, snpid):
        '''
        Extract genotype dosage prob
        and convert to expected dosage of second allele by
        dosage[:, 1] + 2 * dosage[:, 2].
        The samples labelled as missing will be imputed to population mean.
        '''
        if snpid not in self.id_dict:
            # the id is not valid
            return None
        dosage = self.bgen.read([self.id_dict[snpid]])
        # dosage: nsample x nvariant (=1) x n_allele_combinations
        dosage = dosage[:, 0, 1] + 2 * dosage[:, 0, 2]
        missing = np.isnan(dosage, axis=1)
        dosage[missing] = np.nanmean(dosage)
        return dosage 

    def dosage_generator(self, id_list):
        for snpid in id_list:
            yield self.get_dosage_by_id(snpid)
