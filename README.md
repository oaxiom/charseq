# charseq
Simple (relative term) ChAR-seq processer

## Why?
We used paired end reads, not long single end reads as suggested in the paper.
This has advantages:
1. We can get longer potential reads by using PE300

This has disadvantages:
1. The bridge may not be sequences making RNA/DNA ends unidentifiable.
2. This dramatically increases the complexity of the read processing.

## How?

python char_seq.py read1.fq.gz read2.fq.gz

## What?

