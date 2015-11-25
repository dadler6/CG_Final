'''
Dan Adler, Farhan Damani
Computational Genomics Final

Library for read mapping

Maps each read to a contig, offset in contig, and the optimal alignment of read to the contig.
'''

import math

def _compute_p_s_o_y(x_i,s_i,o_i,y_i):
    '''
        Compute p(s,o,y | . )
        We simplify this computation to
        p(s,o|y*,.) = p(s_i = s | s_-i) * p(o_i = o | s_i = s) * p(x_i, y*_so | s_i = s, o_i = o, b)

    '''
    return _compute_p_s(s_i)*_compute_p_o(o_i)*_compute_p_x_y(x_i,y_i)


def _compute_p_s(s_i):
    '''
        Compute p(s_i = s | s_-i)
    '''
    return 1


def _compute_p_o(o_i):
    '''
        Compute p(o_i = o | s_i = s)
    '''
    return 0.5 # uniform distribution

def _compute_p_x_y(x,y):
    '''
        Compute p(x,y*_so | s_i = s, o_i = 0, b)
        x_i, y_i ~ A(s_i, o_i, b, p_mis)
        SCORE(read_i) = log(0.5) + n_hit * log(1-p_mis) + n_mis * log(p_mis / (|B|-1))

        @param x: read i
        @param y: optimal alignment of x to contigs
    '''
    n_hit = 0
    p_mis = 0.1
    for i in range(len(x)):
        if x[i] == y[i]: n_hit = n_hit + 1
    return math.log(0.5) + n_hit * math.log(1-p_mis) + (len(x)-n_hit)*math.log(p_mis / (3))

def _compute_y_star(x, contigs):
    '''
        Compute y* = arg min (hamming distance) over all possible (s,o)
        @param x: read i
        @param contigs: list of contigs
        @return y* sequence

    '''
    min_hamming = len(x) + 1
    for c in contigs:
        for i in xrange(0,len(c)-len(x)+1): # from 0 to n-k+1
            hamming_dist = _compute_hamming(x,c[i:i+len(x)]) # x and substring starting at index i up to length of x
            if hamming_dist < min_hamming: min_hamming = hamming_dist
    return min_hamming

def _compute_hamming(x,y):
    '''
        Compute hamming distance of x and y
        @param x,y : two sequences
        @return hamming score
    '''
    score = 0
    for i in range(len(x)):
        if x[i] != y[i]: score = score + 1
    return score

def map_reads(reads_dict,contigs):
    reads_dict_new = {}
    return reads_dict