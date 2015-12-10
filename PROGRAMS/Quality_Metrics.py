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



def _driver():
    f = open('../SRC_OUTPUT/trial_one/likelihood.txt', 'r')
    _create_likelihood_figure(_process_likelihood_output_file(f))

