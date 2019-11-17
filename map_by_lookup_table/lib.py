import numpy as np
import sys, os

global coding_dic
# ref/alt coding rule:
m = '''
- A T G C
A - 1 -2 -3
T 1 - -3 -2
G 2 3 - 4
C 3 2 4 -'''
m = m.split('\n')[1:]
m = np.array([ i.split(' ') for i in m ])
x = m[0][1:]
y = m.transpose()[0][1:]
coding_dic = {}
for i in range(len(x)):
    for j in range(len(y)):
        # print(i, j)
        val = m[i + 1][j + 1]
        if x[i] == y[j]:
            continue
        coding_dic[ x[i] + y[j] ] = int(val)


def ref_alt_to_code(ref, alt):
    if ref + alt not in coding_dic:
        return None
    else:
        return coding_dic[ref + alt]        

def my_read(filename):
    filep, filee = os.path.splitext(filename)
    if filee == '.gz':
        return gzip.open(filename, 'rt')
    else:
       return open(filename, 'r')

def read_in_file_name_allowing_by_chromosome(filename):
    out = []
    if '@' in filename:
        file_split = filename.split('@')
        if len(file_split) != 2:
            print('Wrong filename:', filename, ', it contains more than one @!')
            sys.exit()
        for i in range(1, 23):
            out.append(file_split[0] + str(i), file_split[1])
    else:
        out.append(filename)
    return out
    