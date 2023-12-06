from line_profiler import LineProfiler
import unittest
import operator
from ProtParts.utils import hobohm1, read_seq
from ProtParts.Measure import Measure

class TestProfile(unittest.TestCase):
        
        def setUp(self):
            self.measure = Measure()
            self.blastp_exec = 'blastp'
            self.makeblastdb_exec = 'makeblastdb'
            self.tmp_dir = './tmp'
            self.threshold = 1e-9
            self.operator = operator.lt
            self.sequence_file = 'example.fa'

    
        def test_profile_hobohm1(self):


            sequences = read_seq(self.sequence_file)
            measurement = self.measure.blastp(sequences, self.makeblastdb_exec, self.blastp_exec, self.tmp_dir, evalue=10, num_threads=4)
    
            lp = LineProfiler()
            lp_wrapper = lp(hobohm1)
            lp_wrapper(sequences, measurement, self.threshold, self.operator)
            lp.print_stats()    


if __name__ == '__main__':
    unittest.main()