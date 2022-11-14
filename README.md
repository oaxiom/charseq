# charseq
Simple (relative term) ChAR-seq processor

## Why?
We used paired end reads, not long single end reads as suggested in the paper.
This has advantages:
1. We can get longer potential reads by using PE300
2. Mismatches can be fixed

This has disadvantages:
1. The bridge may not be sequences making RNA/DNA ends unidentifiable.
2. This dramatically increases the complexity of the read processing.

## How?

You should use cut_adapt (or equiavlent) to remove adaptors first!

`python char_seq.py read1.fq.gz read2.fq.gz`

## What?

Test on 1M reads:
`
test_split % python ~/Tools/charseq/charseq.py test_fq/Hs_char_all.rp1.trim.p1.fq.gz test_fq/Hs_char_all.rp1.trim.p2.fq.gz
1,000,000
Processed: 1,000,000 reads
Both pairs are too short to give a result: 0 (0.0%)
 At least one of the reads was a homopolymer: 13,507 (1.4%)
  Probable self prime 387 (0.0%)
   Can't overlap the reads: 0 (0.0%)
     PCR dupes: 96,769 (9.7%)
      [ Perfect bridge with 0 bp mismatch: 501,723 (50.2%)
      [ Rescued bridge with 1-2 bp mismatch: 26,553 (2.7%)
      [ Really no bridge: 363,830 (36.4%)
      [ >1 bridge: 65,456 (6.5%)
        With bridge: 525,506 (52.6%)
          [ Bridge on forward: 312,212 (59.4%)
          [ Bridge on reverse: 216,064 (41.1%)
           [ The bridge linker was not cut: 56,683 (5.7%)
           [ DNA sequence <20 bp too short: 322,458 (32.2%)
           [ RNA sequence <20 bp too short: 49,602 (5.0%)
Final number of reads kept: 96,763 (9.7%)
`
Yes, around 10% of valid reads is about right for ChAR-seq.

And it produces two new FASTA files:

Hs_char_all.rp1.dna.fq.gz
Hs_char_all.rp1.rna.fq.gz

## Who?

Andrew P. Hutchins
Southern University of Science and Technology,
Shenzhen, China
oaxiom@gmail.com


