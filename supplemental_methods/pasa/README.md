Pipeline details for updating the house fly gene annotations
======

Notes: in what follows, the initial fastq files are assumed to be in a single directory names fem_inf1.fastq.gz, fem_unf2.fastq.gz, etc.

1. Generate a *de novo* transcriptome assembly from the combined infected and uninfected data using Trinity (version r20131110):

```
export SEQLIST=$(ls fem_*.fastq.gz)
Trinity.pl --seqType fq --JM 100G --single $SEQLIST --output trinout --CPU 6 --min_contig_length 150 --min_kmer_cov 2
```

2. Generate a new annotation for *M. domestica* using the Trinity assemblies and PASA. PASA (version 2) was run with the script launch_pasa.sh, which calls a number of other scripts.

