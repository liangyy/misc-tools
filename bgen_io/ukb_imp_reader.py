import bgen_reader
class UKBReader:
    
    def __init__(self, bgen, metafile, sample, byid='rsid'):
        '''
        Initialize self.id_dict as a dictionary with 
          key being byid = rsid or id; 
          value being the index in genotype file.
        '''
        self.bgen = bgen_reader.read_bgen(bgen, metafile_filepath=metafile, samples_filepath=sample)
        self.variant_df = self.bgen['variants'].compute() 
        self.samples = self.bgen["samples"].tolist()
        self.byid = byid
        self._init_id_dict()
        
    def _init_id_dict(self):
        if self.byid == 'rsid':
            self.id_dict = { self.variant_df.rsid[i]: i for i in range(self.variant_df.shape[0]) }
        elif self.byid == 'id':
            self.id_dict = { self.variant_df.id[i]: i for i in range(self.variant_df.shape[0]) }
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
        res = self.bgen['genotype'][self.id_dict[snpid]].compute()
        dosage = res['probs']
        dosage = dosage[:, 1] + 2 * dosage[:, 2]
        dosage[res['missing']] = dosage[~res['missing']].mean()
        return dosage 

    def dosage_generator(self, id_list):
        for snpid in id_list:
            yield self.get_dosage_by_id(snpid)
