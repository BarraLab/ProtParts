import argparse
from ProtParts.main import clust_partition

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description="Protein clustering and partitioning", formatter_class=argparse.RawTextHelpFormatter)
    argparser.add_argument('-i', action='store', dest='input_file', required=True, help="Input fasta file")
    argparser.add_argument('-c', action='store', dest='threshold_c', type=str, help="Threshold for clustering (use comma , to separate multiple thresholds)")
    argparser.add_argument('--exps', action='store', dest='exp_s', type=int, help="Starting exponent for threshold")
    argparser.add_argument('--expe', action='store', dest='exp_e', type=int, help="Ending exponent for threshold")
    argparser.add_argument('-r', action='store', dest='threshold_r', type=float, default=None, help="Threshold for sequence redundancy reduction.\nNone: skip redundancy reduction\n(Default: None)")
    argparser.add_argument('-p', action='store', dest='num_partitions', type=int, help="Number of partitions. 0: skip partitioning")
    argparser.add_argument('-f', action='store', dest='fmt', default='JSON', choices=['JSON', 'TXT', 'CSV', 'FASTA'], help="Output format\n(Default: JSON)")
    argparser.add_argument('-o', action='store', dest='output_dir', required=True, help="Output directory")
    argparser.add_argument('--prune', action='store_true', dest='prune', help="Pruning clusters to improve clustering performance")
    argparser.add_argument('--makeblastdb', action='store', dest='makeblastdb_exec', help="Path to makeblastdb executable\n(Default: config.MAKEBLASTDB_EXEC)")
    argparser.add_argument('--blastp', action='store', dest='blastp_exec', help="Path to blastp executable\n(Default: config.BLASTP_EXEC)")
    argparser.add_argument('--tmpdir', action='store', dest='tmp_dir', help="Path to temporary directory\n(Default: config.TMP_DIR)")

    args = argparser.parse_args()
    # input_file = args.input_file
    # threshold_c = args.threshold_c
    # threshold_r = args.threshold_r
    # num_partitions = args.num_partitions
    # fmt = args.fmt
    # output_file = args.output_file
    # makeblastdb_exec = args.makeblastdb_exec
    # blastp_exec = args.blastp_exec
    # tmp_dir = args.tmp_dir

    clust_partition(args)
    #clust_partition(input_file, threshold_c, threshold_r, num_partitions, output_file, fmt, makeblastdb_exec, blastp_exec, tmp_dir)
    

