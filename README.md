# charseq
Simple (relative term) ChAR-seq processor

## Why?
We used paired end reads, not long single end reads as suggested in the paper.
This has advantages:
1. We can get longer potential reads by using PE300
2. Mismatches can be fixed

This has disadvantages:
1. The bridge may not be sequenced making RNA/DNA ends unidentifiable.
2. This dramatically increases the complexity of the read processing.

## How?

You should use cutadapt (or equiavlent) to remove adaptors first!

```
python char_seq.py read1.fq.gz read2.fq.gz
python merge_rna_dna_bams.py dnareads.bam rnareads.bam
```

## What?

Test on 1M reads:

```
test_split % python charseq.py test_fq/Hs_char_all.rp1.trim.p1.fq.gz test_fq/Hs_char_all.rp1.trim.p2.fq.gz
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
```

Yes, around 10% of valid reads is about right for ChAR-seq.

And it produces two new FASTA files:

Hs_char_all.rp1.dna.fq.gz
Hs_char_all.rp1.rna.fq.gz

Then do you alignment, using your favourite aligners.
We use bowtie2 for the DNA reads and STAR for the RNA reads. Probably better to use one 
aligner for both sides though.

Then merge the resulting bams so that only reads where both ends aligned are retained
Run merge_rna_dna_bams.py to pair up the reads 

```
Process: ../align.dna/Hs_char_all.rp1.dna.bam
Rejected reads QC: 0 (0.0%)
Process: ../align.rna/Hs_char_all.rp1.rna.Aligned.out.bam
Rejected reads QC: 75,062 (9.9%)
Merging
Found 469,606 matching reads out of 686,251 DNA and 530,610 RNA reads
```

(This above example was done using the entire dataset, the split above was done using a 
subset of reads)

You are on your own from here. We suggest using MACS and/or your favourite read counter 
(we like scTE, but then we would). I think best practice and robust tools haven't 
really been developed yet (writing this: 2022).

## Who?

Andrew P. Hutchins
Southern University of Science and Technology,
Shenzhen, China
oaxiom@gmail.com


