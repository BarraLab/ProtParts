from Bio import SeqIO
import json
import warnings
import operator


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


def reduce_redundancy(sequences, measurement, threshold, operator=operator.lt):
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
    measurement_dict = {(seq1, seq2):value for seq1, seq2, value in measurement}

    unqiue_seq = set()
    for qseq in sequences_id_s:
        keep = True
        for useq_id in unqiue_seq:
            if (qseq.id, useq_id) in measurement_dict:
                measure = measurement_dict[(qseq.id, useq_id)]
            elif (useq_id, qseq.id) in measurement_dict:
                measure = measurement_dict[(useq_id, qseq.id)]
            else:
                measure = 11
            
            if operator(measure, self.threshold):
                keep = False
                break
            else:
                keep = True
        
        if keep:
            unqiue_seq.add(qseq.id)
    
    sequences_r = {seq_id:seq for seq_id, seq in sequences.items() if seq_id in unqiue_seq}

    return sequences_r




