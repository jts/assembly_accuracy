# Tools for assessing the accuracy of a genome assembly

This repository contains a small tool to assess the accuracy of an assembly using minimap2.

It requires the `pyvcf` and `pysam` python packages, and `minimap2` to be on your PATH.

## Usage

```
fastmer.py --reference reference.fasta --assembly assembly.fasta --min-mapping-quality 10
```

## Example output:

```
assembly_name    percent_identity  qscore  num_matches  num_mismatches  num_insertions  num_deletions  4mer_acc  5mer_acc  6mer_acc  7mer_acc  8mer_acc  9mer_acc
nanopolish_r94   99.956154         33.58   4691618      163             730             1165           0.995     0.983     0.923     0.825     0.730     0.545
```
