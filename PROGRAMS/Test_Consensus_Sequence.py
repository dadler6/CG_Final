'''
Dan Adler, Farhan Damani
Computational Genomics Final

Tests for Consensus Sequencing.
'''

import unittest
import Consensus_Sequence as cs 

class Test_Consensus_Seuqence(unittest.TestCase):

    # Test with one contig
    def test_single_contig_1(self):
        contigs_1 = ['ATCGCTGATT']
        reads_1 = {\
        'ATCG':[[0,0]],'TCGC':[[0,1]],'CCCT':[[0,2]],'GCTG':[[0,3]],\
        'TTCA':[[0,4]],'TGAT':[[0,5]],'GATT':[[0,6]]}
        self.assertEqual(contigs_1,cs.consensus_sequence(reads_1))

    # Test with single contig and equal probabilities of certain nts
    def test_single_contig_2(self):
        contigs_2 = ['ATCGCTGATT','ATCCCTCATT','ATCGCTCATT','ATCCCTGATT']
        reads_2 = {\
        'ATCC':[[0,0]],'TCGC':[[0,1]],'CCCT':[[0,2]],'GCTC':[[0,3]],\
        'TTCA':[[0,4]],'TGAT':[[0,5]],'GATT':[[0,6]]}
        self.assertTrue(cs.consensus_sequence(reads_2)[0] in contigs_2)

    # Test with multiple contigs
    def test_multiple_contigs_1(self):
        contigs_3 = ['ATCGCTGATT','TTTACGATGC']
        reads_3 = {\
        'ATCG':[[0,0]],'TCGC':[[0,1]],'CCCT':[[0,2]],'GCTG':[[0,3]],\
        'TTCA':[[0,4]],'TGAT':[[0,5]],'GATT':[[0,6]],\
        'TTAA':[[1,0]],'TTAC':[[1,1]],'TGCG':[[1,2]],'ACGA':[[1,3]],\
        'AGAT':[[1,4]],'GTTG':[[1,5]],'AAGC':[[1,6]],}
        self.assertEqual(contigs_3,cs.consensus_sequence(reads_3))

    # Test with multiple contigs and equal probabilities of certain nts
    def test_multiple_contigs_2(self):
        contigs_4 = ['ATCGCTGATT','ATCCCTCATT','ATCGCTCATT','ATCCCTGATT',\
            'TTTACGATGC','TTTAAGATGC']
        reads_4 = {\
        'ATCG':[[0,0]],'TCGC':[[0,1]],'CCCT':[[0,2]],'GCTG':[[0,3]],\
        'TTCA':[[0,4]],'TGAT':[[0,5]],'GATT':[[0,6]],\
        'TTAA':[[1,0],[1,1]],'TGCG':[[1,2]],'CCGA':[[1,3]],\
        'AGAT':[[1,4]],'GTTG':[[1,5]],'AAGC':[[1,6]]}
        self.assertTrue(cs.consensus_sequence(reads_4)[0] in contigs_4\
            and cs.consensus_sequence(reads_4)[1] in contigs_4)

if __name__ == '__main__':
    unittest.main()


