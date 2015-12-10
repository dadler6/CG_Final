__author__ = 'farhan_damani'

import matplotlib.pyplot as plt
import numpy as np


def _create_likelihood_figure(ll_list):
    iters = []
    for i in range(len(ll_list)): iters.append(i)
    plt.plot(np.array(iters),np.array(ll_list), label ='Likelhood')
    plt.xlabel('Iteration')
    plt.ylabel('-log(likelihood)')
    plt.show()


def _process_likelihood_output_file(f):
    ll = []
    for i,l in enumerate(f):
        if i == 0 or i == 1: continue # skip ll = 0
        w = l.split('\t')
        if len(w) > 0:
            ll.append(float(w[1]))
    return ll



def _percent_change_of_contigs():
    '''
        Evaluate percent chance of each contig over time steps.
        Use edit distance to show this.
    '''

    f = open('../SRC_OUTPUT/trial_three/contig.txt', 'r')
    c_prev = ''
    contig_one_list = []
    contig_two_list = []
    contig_three_list = []
    for l in f:
        w = l.split('\t')
        if len(w) > 0:
           if str(w[1].strip())=='merge': continue
           if str(w[1].strip()) =='end': continue
           contig_one_list.append(str(w[2].strip())), contig_two_list.append(str(w[3].strip())), contig_three_list.append(str(w[4].strip()))

    for i in xrange(0,len(contig_one_list)):
        if i == 0: continue
        print edDistDp(contig_one_list[i], contig_one_list[i-1]), edDistDp(contig_two_list[i], contig_two_list[i-1]), \
            edDistDp(contig_three_list[i], contig_three_list[i-1])

def edDistDp(x, y):
    """ Calculate edit distance between sequences x and y using
        matrix dynamic programming.  Return distance. """
    D = np.zeros((len(x)+1, len(y)+1), dtype=int)
    D[0, 1:] = range(1, len(y)+1)
    D[1:, 0] = range(1, len(x)+1)
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            delt = 1 if x[i-1] != y[j-1] else 0
            D[i, j] = min(D[i-1, j-1]+delt, D[i-1, j]+1, D[i, j-1]+1)
    return D[len(x), len(y)]

def _driver():
    #f = open('../SRC_OUTPUT/trial_one/likelihood.txt', 'r')
    #_create_likelihood_figure(_process_likelihood_output_file(f))
    _percent_change_of_contigs()

_driver()
