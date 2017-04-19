import unittest
import BisPin_postprocess


class TestMismatchString(unittest.TestCase):
    
    def setUp(self):
        pass
    
    def test_mismatchWithZeroAtFront(self):
        print "WithZeroAtFront"
        reference_sequence = "TAA"
        bisulfite_sequence = "AAA"
        converted_sequence = "AAA"
        cigar_string = "3M"
        md_string = "0 T 2"
        new_md_string = BisPin_postprocess.makeMismatchString(bisulfite_sequence, converted_sequence, reference_sequence, cigar_string, md_string, printOut=True)
        expected_md_string = "0T2"
        self.assertEqual(expected_md_string , new_md_string[0])
        
    def test_mismatchWithZeroAtEnd(self):
        print "WithZeroAtEnd"
        reference_sequence = "AAT"
        bisulfite_sequence = "AAA"
        converted_sequence = "AAA"
        cigar_string = "3M"
        md_string = "2 T 0"
        new_md_string = BisPin_postprocess.makeMismatchString(bisulfite_sequence, converted_sequence, reference_sequence, cigar_string, md_string, printOut=True)
        expected_md_string = "2T0"
        self.assertEqual(expected_md_string , new_md_string[0])
    
    def test_mismatchBig(self):
        print "mismatchBig"
        bisulfite_sequence = "GTGGGTATAAATTTATTGTTTTATGAAGAGATTGTAAATTTTTTTAATTTTTTTTTTTTTTTTTTTGGAAATTTTTGTCCTTG"
        converted_sequence = "GTAAGTATAAATTTATTGTTTTATGAAGAGATTGTAAATTTTTTTAATTTTTTTTTTTTTTTTTTTGGAAATTTTTAATTTTG"
        cigar_string = "2M 1I 27M 1I 3M 1D 11M 2I 19M 5I12M"
        md_string = "5 T 1 T 0 T 21 G 1 ^A 2T5A26T0T5"
        new_md_string = BisPin_postprocess.makeMismatchString(bisulfite_sequence, converted_sequence, bisulfite_sequence, cigar_string, md_string, printOut=True)
        expected_md_string = "2 G 2 T 1 T 0 T 21 G 1 ^A 2T5A26T1C0C3"
        self.assertEqual(expected_md_string.replace(" ", "") , new_md_string[0])
        
    def test_mismatchAtBegin(self):
        print "At Begin"
        reference_sequence = "TTAAA"
        bisulfite_sequence = "GTGGA"
        converted_sequence = "ATAAA"
        cigar_string = "5M"
        md_string = "0T4"
        new_md_string = BisPin_postprocess.makeMismatchString(bisulfite_sequence, converted_sequence, reference_sequence, cigar_string, md_string, printOut=True)
        expected_md_string = "0T1A0A1"
        self.assertEqual(expected_md_string , new_md_string[0])
    
    def test_mismatchEmpty(self):
        print "Empty"
        bisulfite_sequence = ""
        converted_sequence = ""
        cigar_string = ""
        md_string = ""
        new_md_string = BisPin_postprocess.makeMismatchString(bisulfite_sequence, converted_sequence, bisulfite_sequence, cigar_string, md_string, printOut=True)
        expected_md_string = ""
        self.assertEqual(expected_md_string , new_md_string[0])
        
    def test_mismatchAtEnd(self):
        print "mismatchAtEnd"
        reference_sequence = "GTGGT"
        bisulfite_sequence = "GTGGA"
        converted_sequence = "GTGGA"
        cigar_string = "5M"
        md_string = "4T0"
        new_md_string = BisPin_postprocess.makeMismatchString(bisulfite_sequence, converted_sequence, reference_sequence, cigar_string, md_string, printOut=True)
        expected_md_string = "4T0"
        self.assertEqual(expected_md_string , new_md_string[0])
        
    def test_mismatchInsertBegin(self):
        print "Insert Begin"
        bisulfite_sequence = "CCAAG"
        converted_sequence = "CCAAG"
        cigar_string = "2I3M"
        md_string = "0T2"
        new_md_string = BisPin_postprocess.makeMismatchString(bisulfite_sequence, converted_sequence, bisulfite_sequence, cigar_string, md_string, printOut=True)
        expected_md_string = "0T2"
        self.assertEqual(expected_md_string , new_md_string[0])
        
    def test_mismatchManyDeletes(self):
        print "Many Deletes"
        reference_sequence = "CAATGCATTCTACTGAC"
        bisulfite_sequence = "CCAAGGTCTT"
                          #   IMMMDMDIMMMDM
        converted_sequence = "CCAAGGTTTT"
        cigar_string = "1I 3M 1D 1M 3D 1I 3M 5D 1M"
        md_string = "3 ^T 1 ^CAT 3 ^ACTGA 0 C 0"
        new_md_string = BisPin_postprocess.makeMismatchString(bisulfite_sequence, converted_sequence, reference_sequence, cigar_string, md_string, printOut=True)
        expected_md_string = "3^T1^CAT3^ACTGA0C0"
        self.assertEqual(expected_md_string , new_md_string[0])
        
    def test_mismatchDeleteWithZeros(self):
        print "DeleteWithZeros"
        reference_sequence = "CTTT"
        bisulfite_sequence = "CT"
        converted_sequence = "CT"
        cigar_string = "1M2D1M"
        md_string = "0G0^TT0A0"
        new_md_string = BisPin_postprocess.makeMismatchString(bisulfite_sequence, converted_sequence, reference_sequence, cigar_string, md_string, printOut=True)
        expected_md_string = "0G0^TT0A0"
        self.assertEqual(expected_md_string , new_md_string[0])
        
    def test_mismatchMismatchDiff(self):
        print "MismatchDiff"
        bisulfite_sequence = "CT"
        converted_sequence = "TT"
        cigar_string = "2M"
        md_string = "0G1"
        new_md_string = BisPin_postprocess.makeMismatchString(bisulfite_sequence, converted_sequence, bisulfite_sequence, cigar_string, md_string, printOut=True)
        expected_md_string = "0G1"
        self.assertEqual(expected_md_string , new_md_string[0])
        
    def test_mismatchMismatchDiff2(self):
        print "MismatchDiff2"
        bisulfite_sequence = "TC"
        converted_sequence = "TT"
        cigar_string = "2M"
        md_string = "1G0"
        new_md_string = BisPin_postprocess.makeMismatchString(bisulfite_sequence, converted_sequence, bisulfite_sequence, cigar_string, md_string, printOut=True)
        expected_md_string = "1G0"
        self.assertEqual(expected_md_string , new_md_string[0])
        
    def test_mismatchBisulfiteCarries(self):
        print "Bisulfite Carries"
        bisulfite_sequence = "GA"
        converted_sequence = "AA"
        cigar_string = "2M"
        md_string = "0G1"
        new_md_string = BisPin_postprocess.makeMismatchString(bisulfite_sequence, converted_sequence, bisulfite_sequence, cigar_string, md_string, printOut=True)
        expected_md_string = "2"
        self.assertEqual(expected_md_string , new_md_string[0])
        
    def test_mismatchBisulfiteCarries2(self):
        print "Bisulfite Carries2"
        bisulfite_sequence = "GA"
        converted_sequence = "AA"
        cigar_string = "2M"
        md_string = "2"
        new_md_string = BisPin_postprocess.makeMismatchString(bisulfite_sequence, converted_sequence, bisulfite_sequence, cigar_string, md_string, printOut=True)
        expected_md_string = "0G1"
        self.assertEqual(expected_md_string , new_md_string[0])
        
    def test_mismatchBisulfiteCarries3(self):
        print "Bisulfite Carries3"
        bisulfite_sequence = "AAG"
        converted_sequence = "AAA"
        cigar_string = "1M1D2M"
        md_string = "1^A2"
        new_md_string = BisPin_postprocess.makeMismatchString(bisulfite_sequence, converted_sequence, bisulfite_sequence, cigar_string, md_string, printOut=True)
        expected_md_string = "1^A1G0"
        self.assertEqual(expected_md_string , new_md_string[0])
        
    def test_mismatchBisulfiteCarries4(self):
        print "Bisulfite Carries4"
        bisulfite_sequence = "AG"
        converted_sequence = "AA"
        cigar_string = "2M"
        md_string = "1G0"
        new_md_string = BisPin_postprocess.makeMismatchString(bisulfite_sequence, converted_sequence, bisulfite_sequence, cigar_string, md_string, printOut=True)
        expected_md_string = "2"
        self.assertEqual(expected_md_string , new_md_string[0])
        
    def test_mismatchDelWithInserts(self):
        print "DelWithInserts"
        reference_sequence = "TAACTTA"
        bisulfite_sequence = "TACTGACTTA"
        converted_sequence = "TACTGACTTA"
        cigar_string = "1M 3I 2D 2I 4M"
        md_string = "1 ^AA 0 G 0 G 2"
        new_md_string = BisPin_postprocess.makeMismatchString(bisulfite_sequence, converted_sequence, reference_sequence, cigar_string, md_string, printOut=True)
        expected_md_string = "1^AA0G0G2"
        self.assertEqual(expected_md_string , new_md_string[0])
        


if __name__ == '__main__':
    unittest.main()