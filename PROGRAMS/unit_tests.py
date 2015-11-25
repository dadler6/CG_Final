'''
Dan Adler, Farhan Damani
Computational Genomics Final

Unit Tests
'''

import unittest
import Read_Mapping as rm
class Test_Read_Mapping(unittest.TestCase):

    '''
    def _test_compute_p_s_o_y(self):
        return 1
    def _test_compute_p_s(self):
        return 1
    def _test_compute_p_x_y(self):
        return 1
    def _test_compute_y_star(self):
        return 1
    '''
    def _test_compute_hamming(self):
        x = 'ACCT'
        y = 'ATCA'
        hamming_dist = rm._compute_hamming(x,y)
        self.assertEqual(2,hamming_dist)


if __name__ == '__main__':
    unittest.main()
