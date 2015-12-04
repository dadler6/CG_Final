'''
Dan Adler, Farhan Damani
Computational Genomics Final

Tests for the Merge_Contigs.py file.
'''

import unittest
import Merge_Contigs as mc

class Test_Merge_Contigs(unit.TestCase):

    # Test whether two contigs are merged
    def two_contig_test_1(self):
        pass


    # Test whether no contigs are merged
    def two_contig_test_2(self):
        pass

    # Test three contigs where two have overlap, but one has more overlap than the other
    def three_contig_test_1(self):
        pass


    # Test three contigs that should all be successively merged
    def three_contig_test_2(self):
        pass

    # Test three contigs that should not be merged at all
    def three_contig_test_3(self):
        pass

if __name__ == '__main__':
    unittest.main()