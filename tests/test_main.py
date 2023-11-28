import unittest
from ProtParts.main import clust_partition


class TestClustPartition(unittest.TestCase):
    
    def test_clust_partition(self):
        sequence_file = './example.fa'
        threshold_c = 1E-10
        threshold_r = None
        num_partition = 5
        output_file = 'test_example_partitions.fa'
        output_format = 'fa'
        makeblastdb_exec = None
        blastp_exec = None
        tmp_dir = None
        
        # We're just testing if the function runs without throwing an error
        try:
            clust_partition(sequence_file, threshold_c, threshold_r, num_partition, output_file, output_format, makeblastdb_exec, blastp_exec, tmp_dir)
        except Exception as e:
            self.fail(f'clust_partition raised an exception: {e}')

if __name__ == '__main__':
    unittest.main()