import unittest
from Tests import TestAlignmentTokenizers
from Tests import TestMismatchString

def mySuite():
    suite = unittest.TestSuite()
    #suite.addTest(unittest.defaultTestLoader.loadTestsFromName("BisFAST_TestAlignmentTokenizers.TestAlignmentTokenizers"))
    # TestAlignmentTokenizers
    suite.addTest(unittest.makeSuite(TestAlignmentTokenizers.TestAlignmentTokenizers))
    suite.addTest(unittest.makeSuite(TestMismatchString.TestMismatchString))
    # suite.addTest(BisFAST_TestAlignmentTokenizers.TestAlignmentTokenizers())
    #suite.addTest()
    return suite

#if __name__ == '__main__':
#unittest.main()
runner = unittest.TextTestRunner()
theSuite = mySuite()
runner.run(theSuite) 