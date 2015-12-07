'''
Dan Adler, Farhan Damani
Computational Genomics Final

Library for read mapping

Maps each read to a contig, offset in contig, and the optimal alignment of read to the contig.
'''

import math
import random

# constant
geometric_prob = 0.7

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
        P(s) is modeled as a geometric distribution
        P(s) = (1-p)^(k-1) * p -> this is the probability that the first occurrence of success
        requires k number of independent trials, each with success probability p.

        @param s_i : contig for read i
        @return: probability of read i belonging to contig s_i

    '''
    return (1-geometric_prob)**(s_i-1)*geometric_prob


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
        @return: SCORE_i(read)
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
    y_star = ''
    for c in contigs:
        for i in xrange(0,len(c)-len(x)+1): # from 0 to n-k+1
            hamming_dist = _compute_hamming(x,c[i:i+len(x)]) # x and substring starting at index i up to length of x
            if hamming_dist < min_hamming:
                min_hamming = hamming_dist
                y_star = c[i:i+len(x)]
    return y_star # ERROR: RETURN Y* SEQUENCE

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



def run(reads_dict, contigs):
    percent_sampling = .8
    key_list = reads_dict.keys()
    for i in xrange(0,int(percent_sampling*len(key_list))):
        candidate_probabilities = {}

        # generate a random read
        key_index = random.randint(0,int(percent_sampling*len(key_list))-1)
        x = key_list[key_index]
        if len(reads_dict.get(x)) > 1: # multiple same reads mapping to different contigs
            read_index = random.randint(0,len(reads_dict.get(x))-1)
        else:
            read_index = 0
        read_mapping_list = reads_dict.get(x)[read_index]
        s,o,p_prev = read_mapping_list[0],read_mapping_list[1],read_mapping_list[2]
        y_star = _compute_y_star(x,contigs)

        for s,contig in enumerate(contigs):
            for o,b_so in enumerate(contig):
                p_new = _compute_p_s_o_y(x,s,o,y_star)
                candidate_probabilities[p_new] = (s,o)
        # take max of candidate probabilities
        max_prob = max(map(float,candidate_probabilities))
        max_s,max_o = candidate_probabilities.get(max_prob)

        u = random.random()
        alpha = min(1,float(max_prob/p_prev))
        if u < alpha: # accept
            reads_dict[x][read_index][0],reads_dict[x][read_index][1],reads_dict[x][read_index][2] = max_s,max_o,max_prob

