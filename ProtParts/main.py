from .Clustering import Clustering
from .Measure import Measure
from .Partitioning import Partitioning
from .utils import read_seq, write_partition, write_cluster, hobohm1, init_logging, remove_duplicate
from .settings import MAKEBLASTDB_EXEC, BLASTP_EXEC, TMP_DIR

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
    logger = init_logging(tmp_dir)
    
    # read sequences
    sequences = read_seq(sequence_file)
    logger.info(f"Number of sequences: {len(sequences)}")

    # remove duplicate sequences
    logger.debug("Removing duplicate sequences...")
    sequences = remove_duplicate(sequences)
    logger.info(f"Number of unique sequences: {len(sequences)}")


    # sequence similarity measurement
    logger.debug("Runing BLASTP...")
    measure = Measure()
    measurement = measure.blastp(sequences, makeblastdb_exec, blastp_exec, tmp_dir, evalue=10, num_threads=4)

    # redundancy reduction
    if threshold_r is not None:
        logger.debug("Reducing redundancy...")
        logger.info(f"Threshold for redundancy reduction: {threshold_r}")
        sequences = hobohm1(sequences, measurement, threshold_r, reduce_redundancy=True)
        logger.info(f"Number of sequences after redundancy reduction: {len(sequences)}")

    # clustering
    logger.debug("Clustering with graph...")
    logger.info(f"Threshold for clustering: {threshold_c}")
    clust = Clustering(threshold=threshold_c, method='graph', measurement_type='distance')
    cluster = clust.clustering(sequences, measurement)
    logger.info(f"Number of clusters: {len(cluster)}")


    # partitioning
    if num_partitions and num_partitions > 0:
        logger.debug("Partitioning...")
        logger.info(f"Number of Partitions: {num_partitions}")
        partitioner = Partitioning(num_partitions=num_partitions, method='random')
        partitions = partitioner.random_partitioning(cluster)
        logger.debug("Writing partitions...")
        write_partition(partitions, output_file, output_format, sequences=sequences, method='graph', threshold=threshold_c)
    else:
        logger.debug("Writing clusters...")
        write_cluster(cluster, output_file, output_format, sequences=sequences, method='graph', threshold=threshold_c)
    
    logger.debug("Done.")

