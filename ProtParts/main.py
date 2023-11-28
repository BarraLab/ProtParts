from .Clustering import Clustering
from .Measure import Measure
from .Partitioning import Partitioning
from .utils import read_seq, write_partition, write_cluster, reduce_redundancy
from .settings import MAKEBLASTDB_EXEC, BLASTP_EXEC, TMP_DIR
import logging

def clust_partition(sequence_file, threshold_c, threshold_r, num_partitions, output_file, output_format, makeblastdb_exec=None, blastp_exec=None, tmp_dir=None):
    """
    Clustering and partitioning

    Parameters
    ----------
    sequence_file : str
        Path to input sequence file
    threshold_c : float
        Threshold for clustering
    threshold_r : float
        Threshold for sequence redundancy reduction. None: skip redundancy reduction
    num_partition : int
        Number of partitions. 0: skip partitioning
    output_file : str
        Path to output file
    output_format : str
        Output format
    makeblastdb_exec : str
        Path to makeblastdb executable
    blastp_exec : str
        Path to blastp executable
    tmp_dir : str
        Path to temporary directory
    """

    # set config
    if makeblastdb_exec is None:
        makeblastdb_exec = MAKEBLASTDB_EXEC
    if blastp_exec is None:
        blastp_exec = BLASTP_EXEC
    if tmp_dir is None:
        tmp_dir = TMP_DIR
    
    # set logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # read sequences
    sequences = read_seq(sequence_file)
    logging.info(f"Number of sequences: {len(sequences)}")


    # sequence similarity measurement
    logging.info("Runing BLASTP...")
    measure = Measure()
    measurement = measure.blastp(sequences, makeblastdb_exec, blastp_exec, tmp_dir, evalue=10, num_threads=4)

    # redundancy reduction
    if threshold_r is not None:
        logging.info("Reducing redundancy...")
        sequences = reduce_redundancy(sequences, measurement, threshold_r)
        logging.info(f"Number of sequences after redundancy reduction: {len(sequences)}")

    # clustering
    logging.info(f"Clustering with graph {threshold_c}...")
    clust = Clustering(threshold=threshold_c, method='graph', measurement_type='distance')
    cluster = clust.clustering(sequences, measurement)
    logging.info(f"Number of clusters: {len(cluster)}")


    # partitioning
    if num_partitions and num_partitions > 0:
        logging.info(f"Partitioning with {num_partitions} partitions...")
        partitioner = Partitioning(num_partitions=num_partitions, method='random')
        partitions = partitioner.random_partitioning(cluster)
        logging.info("Writing partitions...")
        write_partition(partitions, output_file, output_format, sequences=sequences, method='graph', threshold=threshold_c)
    else:
        logging.info("Writing clusters...")
        write_cluster(cluster, output_file, output_format, sequences=sequences, method='graph', threshold=threshold_c)

