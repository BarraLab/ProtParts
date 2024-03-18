# ProtParts

ProtParts clusters and partitions  biological protein sequences for machine learning.

## Installation

```bash
pip install -r requirments.txt
```

### Usage

```txt
$ python protparts.py -h
usage: protparts.py [-h] -i INPUT_FILE [-c THRESHOLD_C [THRESHOLD_C ...]] [-r THRESHOLD_R] [-p NUM_PARTITIONS] [-f {JSON,TXT,CSV,FASTA}] -o OUTPUT_DIR [-s] [--makeblastdb MAKEBLASTDB_EXEC] [--blastp BLASTP_EXEC]
                    [--tmpdir TMP_DIR]

Protein clustering and partitioning

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILE         Input fasta file
  -c THRESHOLD_C [THRESHOLD_C ...]
                        Threshold for clustering
  -r THRESHOLD_R        Threshold for sequence redundancy reduction.
                        None: skip redundancy reduction
                        (Default: None)
  -p NUM_PARTITIONS     Number of partitions. 0: skip partitioning
  -f {JSON,TXT,CSV,FASTA}
                        Output format
                        (Default: JSON)
  -o OUTPUT_DIR         Output directory
  -s                    Separate clusters for improvement
  --makeblastdb MAKEBLASTDB_EXEC
                        Path to makeblastdb executable
                        (Default: config.MAKEBLASTDB_EXEC)
  --blastp BLASTP_EXEC  Path to blastp executable
                        (Default: config.BLASTP_EXEC)
  --tmpdir TMP_DIR      Path to temporary directory
                        (Default: config.TMP_DIR)
```

Clustering with a threshold

```bash
python protparts.py -i example.fa -c 1e-9 -o example_clustered.json
```

Clustering with multiple thresholds

```bash
python protparts.py -i example.fa -c 1e-8 1e-9 1e-10 -o example_clustered.json
```

Performance sequence redundancy reduction with a threshold before clustering

```bash
python protparts.py -i example.fa -c 1e-9 -r 1e-100 -o example_clustered.json
```

Clustering and partitioning

```bash
python protparts.py -i example.fa -c 1e-9 -p 5 -o example_clustered.json
```

Clustering based on the number of partitions 
This option will use the highest E-value threshold (lowest sequence homology) which can fit into the partition capacity for the given number of partitions
```bash
python protparts.py -i example.fa -p 5 -o example_clustered.json
```

Output with specific output format

```bash
python protparts.py -i example.fa -c 1e-9 -f fasta -o example_clustered.fasta
```

Speicify BLAST programs and temporary directory

```bash
python protparts.py -i example.fa -c 1e-9 -o example_clustered.json --makeblastdb blast_program_dir/makeblastdb --blastp  blast_program_dir/blastp --tmpdir your_dir/tmp
```
