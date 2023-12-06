from Bio import SeqIO
import json
import warnings
import operator
import logging

def read_seq(seq_file):
    """
    Read sequences from file

    Parameters
    ----------
    seq_file : str
        Path to sequence file

    Returns
    -------
    sequences : list
        List of sequences
    """
    sequences = dict()
    with open(seq_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            if record.id in sequences:
                warnings.warn(f"Duplicate sequence ID: {record.id}, only the first one is used")
            else:
                sequences[record.id] = record
    return sequences


def read_blastp(blastp_file):
    """
    Read blastp output with outfmt 6

    Parameters
    ----------
    blastp_file : str
        Path to blastp output file

    Returns
    -------
    measurement : tuple
        List of measurement (seq1, seq2, measurement)
    """
    measurement = []
    first_hit = set()
    with open(blastp_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if (line[0], line[1]) in first_hit:
                continue
            else:
                first_hit.add((line[0], line[1]))
                measurement.append((line[0], line[1], float(line[10])))
    return measurement


def remove_duplicate(sequences):
    """
    Remove duplicate sequences

    Parameters
    ----------
    sequences : dict
        Dict of sequences
    
    Returns
    -------
    sequences : dict
        Dict of sequences
    """
    sequences_nodup = dict()
    sequences_uniq = set()
    for seq_id, seq in sequences.items():
        if seq.seq not in sequences_uniq:
            sequences_uniq.add(seq.seq)
            sequences_nodup[seq_id] = seq
    return sequences_nodup

def write_partition(partition, out_file, fmt='json', **kwargs):
    """
    Write partition to file

    Parameters
    ----------
    partition : dict
        Partitioned sequences
    out_file : str
        Path to output file
    fmt : str
        Output format
    kwargs : dict
        Keyword arguments for output format
    """
    
    if fmt == 'txt':
        with open(out_file, 'w') as f:
            f.write(f"# Clustering method: {kwargs['method']}\n")
            f.write(f"# Threshold: {kwargs['threshold']}\n")
            f.write(f"# Number of partitions: {len(partition)}\n")
            for pidx, par in partition.items():
                for cidx, c in par.items():
                    for name in c:
                        f.write(f"ClustID {cidx} PartID {pidx} {name}\n")    
    elif fmt == 'json':
        with open(out_file, 'w') as f:
            partition_named = {f"Partition_{pidx}":{f"Cluster_{cidx}":c for cidx, c in par.items()} for pidx, par in partition.items()}
            json.dump(partition_named, f, indent=4)
    elif fmt == 'csv':
        with open(out_file, 'w') as f:
            f.write('SequenceID,PartitionID,ClusterID\n')
            for pidx, par in partition.items():
                for cidx, c in par.items():
                    for name in c:
                        f.write(f"{name},{pidx},{cidx}\n")
    elif fmt in ('fasta', 'fa'):
        with open(out_file, 'w') as f:
            for pidx, par in partition.items():
                for cidx, c in par.items():
                    for name in c:
                        f.write(f">{name} Cluster_{cidx} Partition_{pidx}\n")
                        f.write(f"{kwargs['sequences'][name].seq}\n")
    else:
        raise ValueError(f"Unknown output format: {fmt}")


def write_cluster(cluster, out_file, fmt='json', **kwargs):
    """
    Write cluster to file

    Parameters
    ----------
    cluster : Cluster
        Cluster object
    out_file : str
        Path to output file
    fmt : str
        Output format
    kwargs : dict
        Keyword arguments for output format
    """

    if fmt == 'txt':
        with open(out_file, 'w') as f:
            f.write(f"# Clustering method: {kwargs['method']}\n")
            f.write(f"# Threshold: {kwargs['threshold']}\n")
            f.write(f"# Number of clusters: {cluster.num_data(by='sum')}\n")
            for cidx, c in cluster.items():
                for name in c:
                    f.write(f"ClustID {cidx} {name}\n")
    elif fmt == 'json':
        with open(out_file, 'w') as f:
            cluster_named = {f"Cluster_{cidx}":c for cidx, c in cluster.items()}
            json.dump(cluster_named, f, indent=4)
    elif fmt == 'csv':
        with open(out_file, 'w') as f:
            f.write('SequenceID,ClusterID\n')
            for cidx, c in cluster.items():
                for name in c:
                    f.write(f"{name},{cidx}\n")
    elif fmt in ('fasta', 'fa'):
        with open(out_file, 'w') as f:
            for cidx, c in cluster.items():
                for name in c:
                    f.write(f">{name} Cluster_{cidx}\n")
                    f.write(f"{kwargs['sequences'][name].seq}\n")
    else:
        raise ValueError(f"Unknown output format: {fmt}")


def hobohm1(sequences, measurement, threshold, op=operator.le, reduce_redundancy=True):
    """
    Redundancy reduction

    Parameters
    ----------
    sequences : dict
        Dict of sequences
    measurement : tuple
        List of measurement (seq1, seq2, measurement)
    threshold : float
        Threshold for redundancy reduction

    Returns
    -------
    sequences : dict
        Dict of sequences
    """
    # hobohm1
    sequences_id_s = sorted(sequences, key=lambda x:len(sequences[x].seq), reverse=True)
    measurement_dict = {(seq1, seq2):measure for seq1, seq2, measure in measurement}

    unique_seq = dict()
    for qseq_id in sequences_id_s:
        keep = True
        if len(unique_seq) == 0:
            unique_seq[qseq_id] = []
            continue
        for useq_id in unique_seq:
            if (qseq_id, useq_id) in measurement_dict:
                measure = measurement_dict[(qseq_id, useq_id)]
            # elif (useq_id, qseq_id) in measurement_dict:
            #     measure = measurement_dict[(useq_id, qseq_id)]
            else:
                measure = 11
            
            if op(measure, threshold):
                keep = False
                unique_seq[useq_id].append(qseq_id)
                break
            else:
                keep = True
        
        if keep:
            unique_seq[qseq_id] = []

    if reduce_redundancy:
        sequences_r = {seq_id:sequences[seq_id] for seq_id in unique_seq}
        return sequences_r
    else:
        return unique_seq

def init_logging(tmp_dir):
    """
    Initialize logging

    Parameters
    ----------
    tmp_dir : str
        Path to temporary directory
    
    Returns
    -------
    logger : logging.Logger
        Logger
    """
    logger = logging.getLogger('protparts')
    logger.setLevel(logging.DEBUG)

    # create stream handler and set level to debug
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    stream_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(stream_handler)

    # create file handler and set level to info
    file_handler = logging.FileHandler(f'{tmp_dir}/tmp.log')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(file_handler)

    return logger



