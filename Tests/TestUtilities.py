import unittest
from Utilities import BisPin_util
from Utilities import Constants

class TestUtilities(unittest.TestCase):
 
    def setUp(self):
        pass
 
    def test_getReadStrFastaStr(self):
        (r_str, f_str) = BisPin_util.getReadStrFastaStr(Constants.CONV_CT_GA_CT)
        self.assertEqual(Constants.CONV_CT_GA, r_str)
        self.assertEqual(Constants.CONV_CT, f_str)
       
        
if __name__ == '__main__':
    unittest.main()