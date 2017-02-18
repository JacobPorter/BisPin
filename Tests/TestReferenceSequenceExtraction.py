import unittest
import BisPin_postprocess
from Utilities import SeqIterator
from Utilities import Constants
from Utilities import BisPin_util

class TestReferenceSequenceExtraction(unittest.TestCase):
    
    ref_location = "ref.fa"
    sam_location = "reads.sam"
    
    def setUp(self):
        self.ref_dictionary = BisPin_postprocess.makeReferenceGenomeDictionary(self.ref_location)
    
    def test_extractedSeq(self):
        sam_iterator = SeqIterator.SeqIterator(self.sam_location, file_type=Constants.SAM)
        for record in sam_iterator:
            key = record[Constants.SAM_KEY_RNAME]
            position = int(record[Constants.SAM_KEY_POS]) 
            seq = record[Constants.SAM_KEY_SEQ]
            seq_len = len(seq)
            token_cigar = BisPin_util.tokenizeCigar(record[Constants.SAM_KEY_CIGAR])
            reference = BisPin_postprocess.getReferenceSequence(self.ref_dictionary, key, position, seq_len, token_cigar)
            self.assertEqual(seq_len, len(reference))
            self.assertEqual(seq, reference)
    

if __name__ == '__main__':
    unittest.main()