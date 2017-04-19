import unittest
from Utilities import BisPin_util

class TestAlignmentTokenizers(unittest.TestCase):
 
    def setUp(self):
        pass
 
    def test_tokenizeCigar01(self):
        cigar_list = [(2, 'M'), (1, 'I'), (27, 'M'), (1, 'D'), (11, 'M')]
        self.assertEqual( BisPin_util.tokenizeCigar('2M1I27M1D11M'), cigar_list)
        
    def test_tokenizeCigar02(self):
        cigar_list = [(99, 'M')]
        self.assertEqual( BisPin_util.tokenizeCigar('99M'), cigar_list)
 
    def test_tokenizeMDtag01(self):
        mdtag_list = [5, 'T', 1, 'T', 0, 'T', 21, 'G', 1, '^A', 2]
        self.assertEqual( BisPin_util.tokenizeMDtag('5T1T0T21G1^A2'), mdtag_list)
        
    def test_tokenizeMDtag02(self):
        mdtag_list = [51]
        self.assertEqual( BisPin_util.tokenizeMDtag('51'), mdtag_list)
    
    def test_tokenizeMDtag03(self):
        mdtag_list = [5, 'T', 100, 'T', 0, 'T']
        self.assertEqual( BisPin_util.tokenizeMDtag('5T100T0T'), mdtag_list)
        
    def test_tokenizeMDtag04(self):
        mdtag_list = [5, 'T', 1, 'T', 0, 'T', 21, 'G', 1, '^ACT', 2]
        value = BisPin_util.tokenizeMDtag('5T1T0T21G1^ACT2')
        self.assertEqual(value , mdtag_list)
        
if __name__ == '__main__':
    unittest.main()