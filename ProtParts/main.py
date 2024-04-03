from .Clustering import Clustering
from .Measure import Measure
from .Partitioning import Partitioning
from .Report import Report
from .utils import read_seq, write_partition, write_cluster, hobohm1, init_logging, remove_duplicate, draw_figures, plot_sizebar #create_report
from .settings import MAKEBLASTDB_EXEC, BLASTP_EXEC, TMP_DIR
import zipfile
import os

#def clust_partition(sequence_file, threshold_c, threshold_r, num_partitions, output_file, output_format, makeblastdb_exec=None, blastp_exec=None, tmp_dir=None):
def clust_partition(args):
    """
    Clustering and partitioning

    Parameters
    ----------
    input_file : str
        Path to input sequence file
    threshold_c : float
        Threshold for clustering
    threshold_r : float
        Threshold for sequence redundancy reduction. None: skip redundancy reduction
    num_partition : int
        Number of partitions. 0: skip partitioning
    output_file : str
        Path to output file
    fmt : str
        Output format
    makeblastdb_exec : str
        Path to makeblastdb executable
    blastp_exec : str
        Path to blastp executable
    tmp_dir : str
        Path to temporary directory
    """

    # set config
    if args.makeblastdb_exec is None:
        args.makeblastdb_exec = MAKEBLASTDB_EXEC
    if args.blastp_exec is None:
        args.blastp_exec = BLASTP_EXEC
    if args.tmp_dir is None:
        args.tmp_dir = TMP_DIR
    
    # set logging
    logger = init_logging(args.tmp_dir)

    # get the absolute path of the input file
    input_file = os.path.abspath(args.input_file)
    input_name = os.path.basename(input_file).split('.')[0]
    
    # read sequences
    sequences = read_seq(input_file)
    num_seq = len(sequences)
    logger.info(f"Number of sequences: {num_seq}")

    # remove duplicate sequences
    logger.debug("Removing duplicate sequences...")
    sequences = remove_duplicate(sequences)
    num_seq_nodup = len(sequences)
    logger.info(f"Number of unique sequences: {num_seq_nodup}")


    # sequence similarity measurement
    logger.debug("Runing BLASTP...")
    measure = Measure()
    measurement = measure.blastp(sequences, args.makeblastdb_exec, args.blastp_exec, args.tmp_dir, evalue=10, num_threads=4)

    # redundancy reduction
    if args.threshold_r is not None:
        logger.debug("Reducing redundancy...")
        logger.info(f"Threshold for redundancy reduction: {args.threshold_r}")
        sequences = hobohm1(sequences, measurement, args.threshold_r, reduce_redundancy=True)
        logger.info(f"Number of sequences after redundancy reduction: {len(sequences)}")
    
    # get the absolute path of the output file
    output_dir = os.path.abspath(args.output_dir)
    
    # clustering and partitioning
    logger.debug("Clustering with graph...")
    file_results = []
    clustering_results = [['Threshold', '# sequences', '# unique sequences', '# remaining sequences', '# clusters', 'Silhouette score', 'Download']]
    only_partition = False
    size_thres_dict = {}

    if args.threshold_c:
        threshold_c = [float(i) for i in args.threshold_c.split(',')]
    elif args.exp_s and args.exp_e:
        exp_min = min(abs(int(args.exp_s)), abs(int(args.exp_e)))
        exp_max = max(abs(int(args.exp_s)), abs(int(args.exp_e)))
        threshold_c = [10 ** -(exp) for exp in range(exp_min, exp_max+1)]
    elif (args.threshold_c is None) and (args.exp_s is None) and (args.exp_e is None) and args.num_partitions:
        threshold_c = [10 ** (-exp) for exp in range(1, 21)]
        only_partition = True
    else:
        raise ValueError("Threshold for clustering is not specified. At least one of the E-value threshold, range of threshold, or number of partitions should be specified.")
    
    for t_c in threshold_c:
        have_partition = False

        logger.info(f"Threshold for clustering: {t_c}")
        clust = Clustering(threshold=t_c, method='graph', measurement_type='distance')
        cluster = clust.clustering(sequences, measurement)
        logger.info(f"Number of clusters: {len(cluster)}")
        
        output_file = os.path.join(output_dir, input_name + f"_{t_c}.{args.fmt.lower()}")

        if args.num_partitions is None:
            logger.debug("Writing clusters...")
            write_cluster(cluster, output_file, args.fmt, sequences=sequences, method='graph', threshold=t_c)
        else:
            logger.debug("Partitioning...")
            logger.info(f"Number of Partitions: {args.num_partitions}")
            
            partitioner = Partitioning(num_partitions=args.num_partitions, num_sequences=len(sequences), method='random')
            partition_size = partitioner.partition_size()
            #logger.debug(f"Partition size: {partition_size}")
            max_partition_size = partition_size[max(partition_size, key=partition_size.get)]

            max_cluster_size = cluster.num_data(by='max')
            #logger.debug(f"Max size of clusters: {max_cluster_size}")
            
            size_thres_dict[t_c] = max_cluster_size

            if max_cluster_size > max_partition_size:
                output_file = "NA"
            else:
                partitions = partitioner.random_partitioning(cluster)
                logger.debug("Writing partitions...")
                write_partition(partitions, output_file, args.fmt, sequences=sequences, method='graph', threshold=t_c)
                have_partition = True
    
        # evaluate silhouette score
        logger.debug("Evaluating silhouette score...")
        silhouette, _ = cluster.silhouette(measurement)
        logger.info(f"Silhouette score: {silhouette:.3f}")

        row = [t_c, num_seq, num_seq_nodup, len(sequences), len(cluster), silhouette.round(3), output_file]
        clustering_results.append(row)

        # draw figures
        logger.debug("Drawing figures...")
        hist_file, silhouette_file = draw_figures(cluster, measurement, output_dir, threshold=t_c)
        file_results.append([t_c, hist_file, silhouette_file, output_file])

        if only_partition and have_partition:
            logger.debug("Drawing size bar...")
            sizebar_file = plot_sizebar(size_thres_dict, max_partition_size, output_dir)
            file_results[-1].append(sizebar_file)
            break

    # create a zip file
    logger.debug("Creating a zip output file...")
    out_zip_file = os.path.join(output_dir, input_name + '_protparts.zip')
    with zipfile.ZipFile(out_zip_file, 'w') as zipout:        
        for files in file_results:
            if files[3] != "NA":
                zipout.write(files[3], arcname=os.path.basename(files[3]), compress_type=zipfile.ZIP_DEFLATED)

    
    # write clustering report
    logger.debug("Creating clustering report...")
    report = Report()
    report.write_params(args)
    report.write_results(clustering_results, out_zip_file)
    report.write_figures(file_results)
    report.save_html(os.path.join(output_dir, input_name + '_protparts_report.html'))


    # create report
    # logger.debug("Creating report...")
    # create_report(cluster, measurement, output_file, threshold_c=threshold_c, threshold_r=threshold_r, num_partitions=num_partitions, num_seq=num_seq, num_seq_nodup=num_seq_nodup)
    
    logger.debug("Done.")

    
"""
    # clustering based on the number of partitions
    if (args.threshold_c is None) and (args.exp_s is None) and (args.exp_e is None):
        logger.debug("Clustering based on the number of partitions...")

        if args.num_partitions is None:
            raise ValueError("Number of partitions is not specified.")
        size_thres_dict = {}
        
        partitioner = Partitioning(num_partitions=args.num_partitions, num_sequences=len(sequences), method='random')
        partition_size = partitioner.partition_size()
        logger.info(f"Partition size: {partition_size}")
        max_partition_size = partition_size[max(partition_size, key=partition_size.get)]
        
        for exp in range(1, 21):
            t_c = 10 ** (-exp)
            logger.debug(f"Threshold for clustering: {t_c}")
            clust = Clustering(threshold=t_c, method='graph', measurement_type='distance')
            cluster = clust.clustering(sequences, measurement)
            logger.debug(f"Number of clusters: {len(cluster)}")

            # max size of clusters
            max_cluster_size = cluster.num_data(by='max')
            logger.debug(f"Max size of clusters: {max_cluster_size}")
            size_thres_dict[t_c] = max_cluster_size

            if max_cluster_size <= max_partition_size:
                partitions = partitioner.random_partitioning(cluster)
                
                logger.debug("Writing partitions...")
                output_file = os.path.join(output_dir, input_name + f"_{t_c}.{args.fmt.lower()}")
                write_partition(partitions, output_file, args.fmt, sequences=sequences, method='graph', threshold=t_c)
                
                logger.debug("Evaluating silhouette score...")
                silhouette, _ = cluster.silhouette(measurement)
                logger.info(f"Silhouette score: {silhouette:.3f}")
                
                row = [t_c, num_seq, num_seq_nodup, len(sequences), len(cluster), silhouette.round(3), output_file]
                clustering_results.append(row)

                break
        
        # draw figure 1: threshold vs max cluster size
        logger.debug("Drawing figure 1...")
        sizebar_file = plot_sizebar(size_thres_dict, max_partition_size, output_dir)

        # draw figure 2: historgram vs silhouette score
        logger.debug("Drawing figures...")
        hist_file, silhouette_file = draw_figures(cluster, measurement, output_dir, threshold=t_c)
        file_results.append([t_c, hist_file, silhouette_file, output_file, sizebar_file])
        
    else:
        
        if args.threshold_c is None:
            exp_min = min(abs(int(args.exp_s)), abs(int(args.exp_e)))
            exp_max = max(abs(int(args.exp_s)), abs(int(args.exp_e)))
            threshold_c = [10 ** -(exp) for exp in range(exp_min, exp_max+1)]
        # clustering based on thresholds
        else:
            threshold_c = [float(i) for i in args.threshold_c.split(',')]
        
        for t_c in threshold_c:

            logger.info(f"Threshold for clustering: {t_c}")
            clust = Clustering(threshold=t_c, method='graph', measurement_type='distance')
            cluster = clust.clustering(sequences, measurement)
            logger.info(f"Number of clusters: {len(cluster)}")

            # partitioning
            #output_name = os.path.basename(input_file).split('.')[0] + f"_{t_c}." + args.fmt.lower()
            output_file = os.path.join(output_dir, input_name + f"_{t_c}.{args.fmt.lower()}")
            
            if args.num_partitions and args.num_partitions > 0:
                logger.debug("Partitioning...")
                logger.info(f"Number of Partitions: {args.num_partitions}")
                
                partitioner = Partitioning(num_partitions=args.num_partitions, num_sequences=len(sequences), method='random')
                partition_size = partitioner.partition_size()
                logger.debug(f"Partition size: {partition_size}")
                max_partition_size = partition_size[max(partition_size, key=partition_size.get)]

                max_cluster_size = cluster.num_data(by='max')
                logger.debug(f"Max size of clusters: {max_cluster_size}")

                if max_cluster_size > max_partition_size:
                    output_file = "NA"
                else:
                    partitions = partitioner.random_partitioning(cluster)
                    logger.debug("Writing partitions...")
                    write_partition(partitions, output_file, args.fmt, sequences=sequences, method='graph', threshold=t_c)
            else:
                logger.debug("Writing clusters...")
                write_cluster(cluster, output_file, args.fmt, sequences=sequences, method='graph', threshold=t_c)

            # evaluate silhouette score
            logger.debug("Evaluating silhouette score...")
            silhouette, _ = cluster.silhouette(measurement)
            logger.info(f"Silhouette score: {silhouette:.3f}")

            row = [t_c, num_seq, num_seq_nodup, len(sequences), len(cluster), silhouette.round(3), output_file]
            clustering_results.append(row)

            # draw figures
            logger.debug("Drawing figures...")
            hist_file, silhouette_file = draw_figures(cluster, measurement, output_dir, threshold=t_c)
            file_results.append([t_c, hist_file, silhouette_file, output_file])
    """


