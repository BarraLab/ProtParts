import subprocess
import os
import shutil

class Measure:

    def __init__(self):
        pass

    def blastp(self, sequences, makeblastdb_exec, blastp_exec, tmp_dir=None, **kwargs):
        """
        Run blastp

        Parameters
        ----------
        sequences : dict
            Dict of sequences
        blastp_exec : str
            Path to blastp executable
        tmp_dir : str
            Path to temporary directory
        kwargs : dict
            Keyword arguments for blastp
        
        Returns
        -------
        measurement : tuple
            List of measurement (seq1, seq2, measurement)
        """
        if not tmp_dir:
            tmp_dir = os.getcwd()
            tmp_dir = os.path.join(tmp_dir, 'tmp')
        
        if not os.path.exists(tmp_dir):
            os.mkdir(tmp_dir)
        
        tmp_seq_file = os.path.join(tmp_dir, 'tmp.seq')
        tmp_blast_file = os.path.join(tmp_dir, 'tmp.blastp.tab')

        with open(tmp_seq_file, 'w') as f:
            for seqid, seq in sequences.items():
                f.write('>{}\n{}\n'.format(seq.id, seq.seq))

        # makeblastdb
        tmp_db_file = os.path.join(tmp_dir, 'tmp.db')
        cmd = [makeblastdb_exec, "-in", tmp_seq_file, "-dbtype", "prot", "-out", tmp_db_file]
        subprocess.run(cmd, shell=False, stdout=subprocess.DEVNULL)
        
        arglist = []
        for key, value in kwargs.items():
            arglist.append('-' + key)
            arglist.append(str(value))
        cmd = [blastp_exec, "-query", tmp_seq_file, "-db", tmp_db_file, "-out", tmp_blast_file, "-outfmt", "6"] + arglist
        # Run blastp without stdout
        subprocess.run(cmd, shell=False, stdout=subprocess.DEVNULL) 
        
        # read blastp output
        # blastp outfmt 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        measurement = []
        with open(tmp_blast_file, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                measurement.append((line[0], line[1], float(line[10])))
        
        # remove temporary directory
        # shutil.rmtree(tmp_dir)


        return measurement
    








