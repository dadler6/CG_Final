'''
Dan Adler, Farhan Damani
Computational Genomics Final

Initialization

    - randomly partition reads into contigs and to positions in each contig
'''
from random import randint
import math
#CONTIG_LENGTH = 10000
CONTIG_LENGTH = 1000
NUM_CONTIGS = 3
DEFAULT_PROB = 1*10**(-100)
K_MER = 100

def _process(f):
    '''
        reads dictionary maps read -> [ [s_1,o_1,p_1],[s_2,o_1,p_2],...[s_n,o_n,p_n] ]
            where n is the number of mappings of different reads from different viruses but are the same sequence
        @param f: file of reads
        @return: reads_dict
    '''

    reads = _read_fa_file(f)
    reads_dict = {} # map read -> [s,o]
    for r in reads:
        if r not in reads_dict:
            reads_dict[r] = []
        l = []
        s = randint(0,NUM_CONTIGS - 1) # randomly assign read to one of the 7 contigs
        o = randint(0,CONTIG_LENGTH-len(r)+1)
        l.append(s),l.append(o),l.append(-math.log(DEFAULT_PROB))
        reads_dict[r].append(l)
    return reads_dict
def _read_fa_file(f):
    '''
        Read fa file and return list of reads
    '''
    reads = []

    for i,l in enumerate(f):
        w = l.strip().split()
        if len(w) > 0:
            if i%2 == 1: reads.append(''.join(w))
    return reads
