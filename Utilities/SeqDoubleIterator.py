import SeqIterator
import Constants

"""
@author: Jacob Porter
@summary: An iterator class for iterating through two sequence record files simultaneously.
@requires: SeqIterator
"""

class SeqDoubleIterator:
    
    def __init__(self, file_name1, file_name2, file_type=Constants.FASTQ, gzip_switch = False):
        self.SeqIterator1 = SeqIterator.SeqIterator(file_name1, file_type=file_type, gzip_switch = gzip_switch)
        self.SeqIterator2 = SeqIterator.SeqIterator(file_name2, file_type=file_type, gzip_switch = gzip_switch)
    
    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        record1 = self.SeqIterator1.next()
        record2 = self.SeqIterator2.next()
        return (record1, record2)
