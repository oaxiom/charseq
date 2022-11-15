
import sys, os
import pysam

if len(sys.argv) > 3:
    print('args: <dna_aling.bam> <rna_align.bam>')

paired_valid_read = {}

def proc_get_reads(bam_file):
    idx = 1 
    __rejected_reads = 0

    ret_reads = {}

    print('Process: {}'.format(bam_file))
    bam = pysam.AlignmentFile(bam_file, 'r')
    for seq in bam:
        idx += 1
        if (idx) % 1e6 == 0: 
            print('{:,}'.format(idx))
            #break

        if seq.is_unmapped or seq.is_duplicate or seq.is_qcfail or int(seq.mapping_quality) < 10:
            __rejected_reads += 1
            continue

        #print(seq)
        #print(seq.query_name)

        ret_reads[seq.query_name] = seq
    idx -= 1

    print('Rejected reads QC: {:,} ({:.1f}%)'.format(__rejected_reads, __rejected_reads/idx * 100))
    return ret_reads, bam

rna, rnatemplate = proc_get_reads(sys.argv[1])
dna, dnatemplate = proc_get_reads(sys.argv[2])

merged = set([])

#print(rna.keys())
#print(dna.keys())

print('Merging')
# merge valids;
for name in dna:
    if name not in rna:
        continue

    merged.add(name)

print('Found {:,} matching reads out of {:,} DNA and {:,} RNA reads'.format(len(merged), len(dna), len(rna)))
basename = os.path.split(sys.argv[1])[1].replace('.dna.bam', '')

dnaout = pysam.AlignmentFile("{}.dna.bam".format(basename), "wb", template=dnatemplate)
rnaout = pysam.AlignmentFile("{}.rna.bam".format(basename), "wb", template=rnatemplate)

for name in merged:
    rnaout.write(rna[name])
    #print(dna[name])
    dnaout.write(dna[name])

rnaout.close()
dnaout.close()


