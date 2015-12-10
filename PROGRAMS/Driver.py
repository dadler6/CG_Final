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

    f = open('../DATA/0_wuhan_input.txt', 'r')
    reads_dict = init._process(f)
    contigs = cs.run_consensus(reads_dict)
    likelihood = 0
    likelihood_new = 100
    contig_file = open('../SRC_OUTPUT/trial_one/contig.txt', 'w+')
    ll_file = open('../SRC_OUTPUT/trial_one/likelihood.txt', 'w+')
    #while abs(abs(likelihood) - abs(likelihood_new)) > 1: # until likelihood converges
    for i in range(1000):
        for c in contigs:
            contig_file.write('%s\tstart%s\t' %(str(i),str(c)))
        contig_file.write('\n')
        contig_file.flush()
        ll_file.write('%s\t%s\t%s\n' %(str(i), str(likelihood), str(len(contigs))))
        ll_file.flush()
        likelihood = likelihood_new
        reads_dict = rm.run(reads_dict, contigs)
        contigs = cs.run_consensus(reads_dict)
        for c in contigs:
            contig_file.write('%s\tmerge%s\t' %(str(i),str(c)))

        contigs, reads_dict = mc.run_merge(contigs,reads_dict) # how do we know if a merge has happened..do we need to know?
        likelihood_new = ll._likelihood(reads_dict,contigs)
        #assert likelihood_new-likelihood > 0 # likelihood is always increasing after a consensus sequence step

    for c in contigs:
        contig_file.write('1000\tend%s\t' %(str(c)))
    ll_file.write('1000\t%s\t%s\n' %(str(likelihood), str(len(contigs))))
    ll_file.flush()

_main()