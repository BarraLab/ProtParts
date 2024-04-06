# ProtParts

ProtParts clusters and partitions  biological protein sequences for machine learning.

## Installation

```bash
pip install -r requirments.txt
```

### Usage

```txt
$ python protparts.py -h
usage: protparts.py [-h] -i INPUT_FILE [-c THRESHOLD_C] [--exps EXP_S]
                    [--expe EXP_E] [-r THRESHOLD_R] [-p NUM_PARTITIONS]
                    [-f {JSON,TXT,CSV,FASTA}] -o OUTPUT_DIR [--prune]
                    [--makeblastdb MAKEBLASTDB_EXEC] [--blastp BLASTP_EXEC]
                    [--tmpdir TMP_DIR]

Protein clustering and partitioning

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILE         Input fasta file
  -c THRESHOLD_C        Threshold for clustering (use comma , to separate multiple thresholds)
  --exps EXP_S          Starting exponent for threshold
  --expe EXP_E          Ending exponent for threshold
  -r THRESHOLD_R        Threshold for sequence redundancy reduction.
                        None: skip redundancy reduction
                        (Default: None)
  -p NUM_PARTITIONS     Number of partitions. 0: skip partitioning
  -f {JSON,TXT,CSV,FASTA}
                        Output format
                        (Default: JSON)
  -o OUTPUT_DIR         Output directory
  --prune               Pruning clusters to improve clustering performance
  --makeblastdb MAKEBLASTDB_EXEC
                        Path to makeblastdb executable
                        (Default: config.MAKEBLASTDB_EXEC)
  --blastp BLASTP_EXEC  Path to blastp executable
                        (Default: config.BLASTP_EXEC)
```

Clustering with a threshold

```bash
python protparts.py -i example.fa -c 1e-9 -o results/
```

Clustering with multiple thresholds

```bash
python protparts.py -i example.fa -c 1e-8,1e-9,1e-10 -o results/
```

Clustering with a range of thresholds in exponent, for exmaple, using threshold between 1e-5 and 1e-10.

```bash
python protparts.py -i example.fa --exps 5 --expe 10 -o results/
```

Performance sequence redundancy reduction with a threshold before clustering

```bash
python protparts.py -i example.fa -c 1e-9 -r 1e-100 -o results/
```

Clustering and partitioning

```bash
python protparts.py -i example.fa -c 1e-9 -p 5 -o results/
```

Clustering based on the number of partitions 
This option will use the highest E-value threshold (lowest sequence homology) which can fit into the partition capacity for the given number of partitions

```bash
python protparts.py -i example.fa -p 5 -o results/
```

Clustering with a threshold and prune the result clusters to improve clustering performance

```bash
python protparts.py -i example.fa -c 1e-9 --prune -o results/
```

Output with specific output format

```bash
python protparts.py -i example.fa -c 1e-9 -f FASTA -o results/
```

Speicify BLAST programs and temporary directory

```bash
python protparts.py -i example.fa -c 1e-9 -o results/ --makeblastdb blast_program_dir/makeblastdb --blastp  blast_program_dir/blastp --tmpdir your_dir/tmp
```

### Results

ProtParts will create a report of clustering result in html format under the result directory, which contains parameters for clustering and partitioning, stastical description of clusters, and graphical analysis of clusters.

#### Output format

##### JSON

The JSON has python dictionary-like format. The sequence ID can be accessed by partition or cluster index.

```json
{
    "Cluster_0": [
        "A0002",
        "A0003",
        "A0004",
        ...
```

or

```json
{
    "Partition_0": {
        "Cluster_4": [
            "A0010",
            "A0118",
            "A0119",
            ...
```

##### TXT

The TXT file contains basic information of the clustering results, starting with `#`. The `ClustID 0` indicates the numbering of clusters, which is followed by sequence ID. If the number of partitions is provides, extra information `PartID 0` will be appended after the `ClustID`, showing the partition numbering.

```txt
# Clustering method: graph
# Threshold: 1e-09
# Number of clusters: 2030
ClustID 0 A0002
ClustID 0 A0003
ClustID 0 A0004
...
```

or 

```txt
# Clustering method: graph
# Threshold: 1e-09
# Number of partitions: 5
ClustID 4 PartID 0 A0010
ClustID 4 PartID 0 A0118
ClustID 4 PartID 0 A0119
...
```

##### CSV

The CSV file consists of `SequenceID`, `ClusterID` or optional `PartitionID`.

```csv
SequenceID,ClusterID
A0002,0
A0003,0
A0004,0
...
```

or 

```csv
SequenceID,PartitionID,ClusterID
A0010,0,4
A0118,0,4
A0119,0,4
...
```

##### FASTA

The clustering and partitioning results are added to the description line after protein IDs in FASTA file.

```fasta
>A0002 Cluster_0
MAQLTLLLLSLFLTLISLPPPGASISSCNGPCRDLNDCDGQLICIKGKCNDDPEVGTHICGGTTPSPQPGSCNPSGTLTCQGKSYPTYDCSPPVTSSTPAKLTNNDFSEGGDGGGPSECDESYHSNNERIVALSTGWYNGGSRCGKMIRITASNGKSVSAKVVDECDSRHGCDKEHAGQPPCRNNIVDGSNAVWSALGLDKNVGVVDITWSMA
>A0003 Cluster_0
MAQLTLLLLSLFFTLISLPPPGASISSCNGPCRDLNDCNGQLICIKGKCNDDPEVGTHICGGTTPSPQPGSCKPSGTLTCQGKSYPTYDCSPPVTSSTPAKLTNNDFSEGGDGGGPSECDESYHSNNERIVALSTGWYNGGSRCGKMIRITASNGKSVSAKVVDECDSRHGCDKEHAGQPPCRNNIVDGSNAVWSALGLDKNVGVVDITWSMA
>A0004 Cluster_0
MAQLTLLLLSLFLTLISLPPPGASISSCNGPCRDLNDCDGQLICIKGKCNDDPEVGTHICGGTTPSPQPGGCNPSGTLTCQGKSYPTYDCSPPVTSSTPAKLTNNDFSEGGDGGGPSECDESYHSNNERIVALSTGWYNGGSRCGKMIRITASNGKSVSAKVVDECDSRHGCDKEHAGQPPCRNNIVDGSNAVWSALGLDKNVGVVDITWSMA
...
```

or 

```fasta
>A0010 Cluster_4 Partition_0
MARPSFLSLVSLSLLVLSHSSAANRQPSKYQQQQKGECQIQRLNAQEPQQRIQAEAGVTEFWDWTDDQFQCAGVAACRNMIQPRGLLLPSYTNAPTLIYILKGRGITGVMIPGCPETYQSSQQSREGDVSHRQFRDQHQKIRRFQQGDVIALPAGVAHWCYNDGDSDLVTVSVEDTGNRQNQLDNNPRRFFLAGNPQQQQKEMYAKRPQQQHSGNVFRGFDTEVLAETFGVDMEMARRLQGKDDYRGHIIQVERELKIVRPPRTREEQEQQERGERDNGMEETICTARLVENIDNPSRADIFNPRAGRLTSVNSFNLPILNYLRLSAEKGVLYRNALMPPHWKLNAHCVLYATRGEAQMQIVDQRGEAVFNDRIREGQLVVVPQNFVVMKQAGNQGFEWVAIKTNENAMFNTLAGRTSALRAMPVDVLANAYQISQSEARRLKMGREEAVLFEPRSEGRDVD
>A0118 Cluster_4 Partition_0
PPTKFSFSLFLVSVLVLCLGFALAKIDPELKQCKHQCKVQRQYDEQQKEQCVKECEKYYKEKKGREREHEEEEEEWGTGGVDEPSTHEPAEKHLSQCMRQCERQEGGQQKQLCRFRCQERYKKERGQHNYKREDDEDEDEDEAEEEDENPYVFEDEDFTTKVKTEQGKVVLLPKFTQKSKLLHALEKYRLAVLVANPQAFVVPSHMDADSIFFVSWGRGTITKILENKRESINVRQGDIVSISSGTPFYIANNDENEKLYLVQFLRPVNLPGHFEVFHGPGGENPESFYRAFSWEILEAALKTSKDTLEKLFEKQDQGTIMKASKEQVRAMSRRGEGPKIWPFTEESTGSFKLFKKDPSQSNKYGQLFEAERIDYPPLEKLDMVVSYANITKGGMSVPFYNSRATKIAIVVSGEGCVEIACPHLSSSKSSHPSYKKLRARIRKDTVFIVPAGHPFATVASGNENLEIVCFEVNAEGNIRYTLAGKKNIIKVMEKEAKELAFKMEGEEVDKVFGKQDEEFFFQGPEWRKEKEGRADE
>A0119 Cluster_4 Partition_0
MGPPTKFSFSLFLVSVLVLCLGFALAKIDPELKQCKHQCKVQRQYDEQQKEQCVKECEKYYKEKKGREREHEEEEEEWGTGGVDEPSTHEPAEKHLSQCMRQCERQEGGQQKQLCRFRCQERYKKERGQHNYKREDDEDEDEDEAEEEDENPYVFEDEDFTTKVKTEQGKVVLLPKFTQKSKLLHALEKYRLAVLVANPQAFVVPSHMDADSIFFVSWGRGTITKILENKRESINVRQGDIVSISSGTPFYIANNDENEKLYLVQFLRPVNLPGHFEVFHGPGGENPESFYRAFSWEILEAALKTSKDTLEKLFEKQDQGTIMKASKEQIRAMSRRGEGPKIWPFTEESTGSFKLFKKDPSQSNKYGQLFEAERIDYPPLEKLDMVVSYANITKGGMSVPFYNSRATKIAIVVSGEGCVEIACPHLSSSKSSHPSYKKLRARIRKDTVFIVPAGHPFATVASGNENLEIVCFEVNAEGNIRYTLAGKKNIIKVMEKEAKELAFKMEGEEVDKVFGKQDEEFFFQGPEWRKEKEGRADE
...
```

