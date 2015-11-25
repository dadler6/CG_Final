'''
Dan Adler, Farhan Damani
Computational Genomics Final

Driver Program
'''
import Initialization as init
import Consensus_Sequence as cs
import Read_Mapping as rm
import Merge_Contigs as mc
import Likelihood as ll

def _main():

    f = open('../data/reads.fa', 'r')
    reads_dict = init._process(f)
    contigs = cs.consensus_sequence(reads_dict)
    likelihood = 0
    likelihood_new = 100
    while abs(abs(likelihood) - abs(likelihood_new)) > 1: # until likelihood converges
        likelihood = likelihood_new
        reads_dict = rm.map_reads(reads_dict, contigs)
        contigs = cs.consensus_sequence(reads_dict)
        contigs = mc.merge_check_global(contigs,reads_dict) # how do we know if a merge has happened..do we need to know?
        likelihood_new = ll._likelihood(reads_dict,contigs)


        assert likelihood_new-likelihood > 0 # likelihood is always increasing after a consensus sequence step