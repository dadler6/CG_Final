'''
Dan Adler, Farhan Damani
Computational Genomics Final

Initialization

    - randomly partition reads into contigs and to positions in each contig
'''
from random import randint

CONTIG_LENGTH = 10000

def _process(f):
    reads = _read_fa_file(f)
    reads_dict = {} # map read -> [s,o]
    for r in reads:
        reads_dict[r] = []
        s = randint(1,7) # randomly assign read to one of the 7 contigs
        o = randint(0,CONTIG_LENGTH-len(r)+1)
        reads_dict[r].append(s)
        reads_dict[r].append(o)
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

_process()