from pyliftover import LiftOver
import pandas as pd
import numpy as np



def liftover(chr, pos, chainfile):
    # chr: number or chrN
    # pos: 1-base position
    
    lo = LiftOver(chainfile)
    
    # formatting chromosome
    if (not isinstance(chr[0], str)) or ('chr' not in chr[0]):
        chr = [ 'chr' + str(i) for i in chr ]
        
    pos = pos - 1  # pyliftover uses base-0
    lo_out = [ _tidy_liftover(i, j, lo) for i, j in zip(chr, pos) ]
    
    
    out = pd.DataFrame(lo_out, columns = ['liftover_chr', 'liftover_pos'])
    out.iloc[:,1] = out.iloc[:,1] + 1  # convert back to base-1
    return out
    

def _tidy_liftover(chr, pos, lo):
    # print(chr, pos)
    tmp = lo.convert_coordinate(chr, pos)
    if tmp is not None and len(tmp) > 0:
        return tmp[0][:2]
    else:
        return ('chr0', -1)

def ref_alt_to_code(ref, alt):
    if ref + alt not in coding_dic:
        return None
    else:
        return coding_dic[ref + alt]        
