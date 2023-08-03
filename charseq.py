
from itertools import product
from collections import deque
import sys, os, gzip
import regex

compdict = {'A': 'T',
            'C': 'G',
            'G': 'C',
            'T': 'A',
            'N': 'N',}

def needleman_wunsch(x, y):
    """Run the Needleman-Wunsch algorithm on two sequences.

    x, y -- sequences.

    Code based on pseudocode in Section 3 of:

    Naveed, Tahir; Siddiqui, Imitaz Saeed; Ahmed, Shaftab.
    "Parallel Needleman-Wunsch Algorithm for Grid." n.d.
    https://upload.wikimedia.org/wikipedia/en/c/c4/ParallelNeedlemanAlgorithm.pdf
    """
    N, M = len(x), len(y)
    s = lambda a, b: int(a == b)

    DIAG = -1, -1
    LEFT = -1, 0
    UP = 0, -1

    # Create tables F and Ptr
    F = {}
    Ptr = {}

    F[-1, -1] = 0
    for i in range(N):
        F[i, -1] = -i
    for j in range(M):
        F[-1, j] = -j

    option_Ptr = DIAG, LEFT, UP
    for i, j in product(range(N), range(M)):
        option_F = (
            F[i - 1, j - 1] + s(x[i], y[j]),
            F[i - 1, j] - 1,
            F[i, j - 1] - 1,
        )
        F[i, j], Ptr[i, j] = max(zip(option_F, option_Ptr))

    # Work backwards from (N - 1, M - 1) to (0, 0)
    # to find the best alignment.
    alignment = deque()
    i, j = N - 1, M - 1
    while i >= 0 and j >= 0:
        direction = Ptr[i, j]
        if direction == DIAG:
            element = i, j
        elif direction == LEFT:
            element = i, None
        elif direction == UP:
            element = None, j
        alignment.appendleft(element)
        di, dj = direction
        i, j = i + di, j + dj
    while i >= 0:
        alignment.appendleft((i, None))
        i -= 1
    while j >= 0:
        alignment.appendleft((None, j))
        j -= 1

    p1 = [None if i is None else x[i] for i, _ in alignment]
    p2 = [None if j is None else y[j] for _, j in alignment]
    return ''.join(i if i else j for i, j in zip(p1, p2))

    #return list(alignment)

def print_alignment(x, y, alignment):
    print("".join("-" if i is None else x[i] for i, _ in alignment))
    print("".join("-" if j is None else y[j] for _, j in alignment))

def get_consensus(x, y, alignment):
    p1 = [None if i is None else x[i] for i, _ in alignment]
    p2 = [None if j is None else y[j] for _, j in alignment]
    return ''.join(i if i else j for i, j in zip(p1, p2))

def score_matches(x, y, alignment):
    # Retuns the number of base pairs matching;
    p1 = [None if i is None else x[i] for i, _ in alignment]
    p2 = [None if j is None else y[j] for _, j in alignment]
    return sum(1 if i ==j else 0 for i, j in zip(p1, p2))

'''
# Test cases:
# Aligns,
p1 = 'GTTCCATATTAGCTAAACCGGCGTCCAAGGATCAGGCAGAAAAAAGGAAATGTTTGGAATTGCTGGTGTTTGAGTGGGGACTTGTGTAAACTGTTACTAAACTGGAGCAAGTATGTCCTGTCTAA'
p2 = 'CGTGAAGCGTTCCATATTAGCTAAACCGGCGTCCAAGGATCAGGCAGAAAACAGGAAATGTTTGGAATTGCTGGTGTTGGAGTGGGGACTTGTGTAAACTGTGACTAAACTGGAGCAAGTATGTCCTGTCTAA'

a = needleman_wunsch(p1, p2)
print_alignment(p1, p2, a)
print(get_consensus(p1, p2, a))
print(score_matches(p1, p2, a))

# No overlap;
p1 = 'TTTGCAAATAGTCATGCTTTCTCTCATAGTTGTAATGCTTATATCTGTTTAATGTCTTGTTGTTTTGTAAACCTGAAATAACATTTTTCGATTTTAATTTTGATACAGGTGCCTTCAGTTACTCACCATAAATAAAATGTTG'
p2 = 'GCCTTCAGTTACTCACCATAAATAAAATGTTGTCTAATTGAGTATGATAGGCTTCCTTACCCTAAATATTATTTATTTAAAGAGATAACTTATTCTTGATTTTCTAATTGGGTGGTATTGAAACCTATCATATAAATTTTCAATTCTGT'

a = needleman_wunsch(p1, p2)
print_alignment(p1, p2, a)
print(get_consensus(p1, p2, a))
print(score_matches(p1, p2, a))
'''


def rc(seq):
    return ''.join(reversed([compdict[i] for i in seq])) # new list

def fastqPE(filename1, filename2):
    oh1 = gzip.open(filename1, "rt")
    oh2 = gzip.open(filename2, "rt")

    name1 = "dummy"
    while name1 != "":
        name1 = oh1.readline().strip()
        seq1 = oh1.readline().strip()
        strand1 = oh1.readline().strip()
        qual1 = oh1.readline().strip()

        name2 = oh2.readline().strip()
        seq2 = oh2.readline().strip()
        strand2 = oh2.readline().strip()
        qual2 = oh2.readline().strip()

        yield ({"name": name1, "strand": strand1, "seq": seq1, "qual": qual1},
            {"name": name2, "strand": strand2, "seq": seq2, "qual": qual2})

    oh1.close()
    oh2.close()
    return

class stats:
    both_pairs_too_short = 0
    with_bridge = 0
    no_bridge = 0
    multi_bridge = 0
    perfect_bridge = 0
    cant_overlap_reads = 0
    f_strand = 0
    r_strand = 0
    dna_too_short = 0
    rna_too_short = 0
    rescued_bridge = 0
    pcr_dupe = 0
    homopolymer = 0
    probable_self_primes = 0
    no_bridge_linker_cut = 0


def find_bridge(seq, qual, stats=stats):
    find_f = regex.findall("(ACCGGCGTCCAAG)", seq)
    find_r = regex.findall("(CTTGGACGCCGGT)", seq)

    if len(find_f) > 1 or len(find_r) > 1:
        stats.multi_bridge += 1
        return None, full_seq, None, 'multi_bridge1'

    elif len(find_f) == 1 and len(find_r) == 1:
        stats.multi_bridge += 1
        return None, full_seq, None, 'multi_bridge2'

    elif len(find_f) > 0:
        stats.perfect_bridge += 1
        stats.f_strand += 1
        newseq = full_seq
        newqual = qual
        bridge_loc = newseq.find('ACCGGCGTCCAAG')

    elif len(find_r) > 0:
        stats.perfect_bridge += 1
        stats.r_strand += 1
        newseq = rc(full_seq)
        newqual = qual[::-1]
        bridge_loc = newseq.find('ACCGGCGTCCAAG') # Get the new F strand location;

    else:
        # slower mismatch path for the NW spliced read;
        hit_f = regex.findall("(ACCGGCGTCCAAG){s<=3}", seq)
        hit_r = regex.findall("(CTTGGACGCCGGT){s<=3}", seq)

        if len(hit_f) > 1 or len(hit_r) > 1:
            stats.multi_bridge += 1
            return None, full_seq, None, 'multi_bridge3'
        elif len(hit_f) >= 1 and len(hit_r) >= 1:
            stats.multi_bridge += 1
            return None, full_seq, None, 'multi_bridge4'
        elif hit_f or hit_r: # Found one
            stats.rescued_bridge += 1
            if hit_f:
                stats.f_strand += 1
                newseq = full_seq
                newqual = qual
                bridge_loc = seq.find(hit_f[0])
            elif hit_r:
                stats.r_strand += 1
                newseq = rc(full_seq)
                newqual = qual[::-1]
                bridge_loc = newseq.find(rc(hit_r[0]))
            #print(hit_f, hit_r, bridge_loc)
        else:
            return None, full_seq, None, 'no_bridge' # For debugging

    if bridge_loc == -1:
        # Should be impossible to reach here
        1/0

    return bridge_loc, newseq, newqual, None

def commonOverlapIndexOf(text1, text2):
    # Cache the text lengths to prevent multiple calls.
    text1_length = len(text1)
    text2_length = len(text2)
    # Eliminate the null case.
    if text1_length == 0 or text2_length == 0:
        return 0
    # Truncate the longer string.
    if text1_length > text2_length:
        text1 = text1[-text2_length:]
    elif text1_length < text2_length:
        text2 = text2[:text1_length]
    # Quick check for the worst case.
    if text1 == text2:
        return min(text1_length, text2_length)

    # Start by looking for a single character match
    # and increase length until no match is found.
    best = 0
    length = 1
    while True:
        pattern = text1[-length:]
        found = text2.find(pattern)
        if found == -1:
            return best
        length += found
        if text1[-length:] == text2[:length]:
            best = length
            length += 1

def commonOverlapIndexOf_mismatch(text1, text2):
    # Slower mismatch RE-based versions

    text1_length = len(text1)
    text2_length = len(text2)
    if text1_length == 0 or text2_length == 0:
        return 0

    if text1 == text2:
        return min(text1_length, text2_length)

    best = 0
    length = text1_length # must be more than the number of mismatches, otherwise just finds everything;
    all_hits = []
    while length > 3:
        pattern = text1[:-1]
        text2 = text2[:-1]
        length -= 1

        hit = regex.findall("(%s){s<=1}" % pattern, text2) # slow!!
        all_hits += hit

        #print(hit)
        if hit:
            best = text2.find(hit[0])
            #print('Hit!', hit, best)

    # fell off the end, nothing found;
    return best

def get_read_overlap(bc1, bc2, r1_seq, r2_seq, r1_qual, r2_qual):
    left1 = r1_seq[0:bc1]
    left2 = r2_seq[0:bc2]
    if len(left1) >= len(left2):
        left = left1
        left_qual = r1_qual[0:bc1]
    else: # left2 > left1
        left = left2
        left_qual = r2_qual[0:bc2]

    rite1 = r1_seq[bc1:]
    rite2 = r2_seq[bc2:]

    # TODO: Take the highest quality
    if len(rite1) >= len(rite2):
        rite = rite1
        rite_qual = r1_qual[bc1:]
    else: # left2 > left1
        rite = rite2
        rite_qual = r2_qual[bc2:]

    return left, rite, left_qual, rite_qual

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('charseq.py: PE1.fastq.gz PE2.fastq.gz')
        sys.exit()

    print('charseq.py')
    print('Inputs:' sys.argv[1], sys.argv[2])

    stub = os.path.split(sys.argv[1])[1].replace('.fastq', '').replace('.fq', '').replace('.gz', '').replace('.trim', '').replace('.p1','')

    min_read_size = 18
    rna_reads = []
    dna_reads = []
    pcr_dupes = set([])

    # Avoid global gubbins
    idx = 0
    for r1, r2 in fastqPE(sys.argv[1], sys.argv[2]):
        idx += 1
        if (idx) % 1e5 == 0:
            print(f'Processed {idx:,} read pairs'.format(idx))
            #break

        r1_seq = r1['seq']
        r2_seq = rc(r2['seq'])
        r1_qual = r1['qual']
        r2_qual = r2['qual'][::-1]

        if len(r1_seq) + len(r2_seq) < 55: # No way to get info out of this small
            stats.both_pairs_too_short += 1
            continue

        # Can reject the simple case when one read is a homopolymer
        if len(set(r1_seq)) == 1 or len(set(r2_seq)) == 1:
            stats.homopolymer += 1
            continue

        # See if there is a bridge in at least one read:
        b1f = regex.findall("(ACCGGCGTCCAAG)", r1_seq) # bridge must be on same strand
        b2f = regex.findall("(ACCGGCGTCCAAG)", r2_seq)
        b1r = regex.findall("(CTTGGACGCCGGT)", r1_seq) # Not optimised. Done for clarity
        b2r = regex.findall("(CTTGGACGCCGGT)", r2_seq)

        if len(b1f) > 1 or len(b2f) > 1 or len(b1f) > 1 or len(b2f) > 1: # Reject reads with Multibridges
            stats.multi_bridge += 1
            continue

        # This is written in a very verbose style...
        # TODO: Clean up, it's very fast so the method overhead is not a problem.
        if b1f and b2f: # Simple case, read 1 and 2 have a bridge on the top strand
            # Sometimes the rite is longer than the left... so take the longer of the two sides;
            bc1 = r1_seq.find(b1f[0])
            bc2 = r2_seq.find(b2f[0])

            left, rite, left_qual, rite_qual = get_read_overlap(bc1, bc2, r1_seq, r2_seq, r1_qual, r2_qual)

            #print(f'R1: {r1_seq}')
            #print(f'R2: {r2_seq}')
            #print(f'FS: {left}-{rite}')
            #print()

            full_seq = f'{left}{rite}'
            full_qual = f'{left_qual}{rite_qual}' # You should take the highest quality for each nucleotide...

        elif b1r and b2r: # Check the reverse:
            # Simple case, 1 bridge on each read, - strand;
            # Same code as above, just on the reversed  seq:
            # Sometimes the rite is longer than the left... so take the longer of the two sides;
            bc1 = r1_seq.find(b1r[0])
            bc2 = r2_seq.find(b2r[0])

            left, rite, left_qual, rite_qual = get_read_overlap(bc1, bc2, r1_seq, r2_seq, r1_qual, r2_qual)

            #print(f'R1: {r1_seq}')
            #print(f'R2: {r2_seq}')
            #print(f'FS: {left}-{rite}')
            #print()

            # reverse them here;
            full_seq = rc(f'{left}{rite}')
            full_qual = f'{left_qual}{rite_qual}'[::-1] # You should take the highest quality for each nucleotide...
        else:
            # Bridge not found in both reads, see if we can rescue the bridge
            # by looking for mismatch

            if not b1f: b1f = regex.findall("(ACCGGCGTCCAAG){s<=2}", r1_seq) # Using more than 3 is not a good idea.
            if not b2f: b2f = regex.findall("(ACCGGCGTCCAAG){s<=2}", r2_seq)
            if not b1r: b1r = regex.findall("(CTTGGACGCCGGT){s<=2}", r1_seq) # Not optimised. Done for clarity
            if not b2r: b2r = regex.findall("(CTTGGACGCCGGT){s<=2}", r2_seq)

            if len(b1f) > 1 or len(b2f) > 1 or len(b1f) > 1 or len(b2f) > 1: # Reject reads with Multibridges
                stats.multi_bridge += 1
                continue

            if b1f and b2f:
                bc1 = r1_seq.find(b1f[0])
                bc2 = r2_seq.find(b2f[0])

                left, rite, left_qual, rite_qual = get_read_overlap(bc1, bc2, r1_seq, r2_seq, r1_qual, r2_qual)
                # Fix for perfect barcode
                rite = rite[13:]
                rite_qual = rite_qual[13:]

                #print(f'R1: {r1_seq}')
                #print(f'R2: {r2_seq}')
                #print(f'FS: {left}-ACCGGCGTCCAAG-{rite}')
                #print()

                full_seq = f'{left}ACCGGCGTCCAAG{rite}'
                full_qual = f'{left_qual}{rite_qual}' # You should take the highest quality for each nucleotide...
                stats.rescued_bridge += 1

            elif b1r and b2r:
                bc1 = r1_seq.find(b1r[0])
                bc2 = r2_seq.find(b2r[0])

                left, rite, left_qual, rite_qual = get_read_overlap(bc1, bc2, r1_seq, r2_seq, r1_qual, r2_qual)
                # Fix for perfect barcode
                rite = rite[13:]
                rite_qual = rite_qual[13:]

                #print(f'R1: {r1_seq}')
                #print(f'R2: {r2_seq}')
                #print(f'FS: {left}-CTTGGACGCCGGT-{rite}')
                #print()

                full_seq = rc(f'{left}CTTGGACGCCGGT{rite}')
                full_qual = f'{left_qual}{rite_qual}'[::-1] # You should take the highest quality for each nucleotide...
                stats.rescued_bridge += 1

            else:
                # If there is one bridge, then just use that read, and it seems reasonable to add the opposite read onto the 3' end;
                if len(b1f) == 1:

                    full_seq = f'{r1_seq}{r2_seq}'
                    full_qual = f'{r1_qual}{r2_qual}'
                elif len(b2f) == 1:

                    full_seq = f'{r1_seq}{r2_seq}'
                    full_qual = f'{r1_qual}{r2_qual}'
                elif len(b1r) == 1:

                    full_seq = f'{rc(r1_seq)}{rc(r2_seq)}'
                    full_qual = f'{r1_qual[::-1]}{r2_qual[::-1]}'
                elif len(b2r) == 1:

                    full_seq = f'{rc(r1_seq)}{rc(r2_seq)}'
                    full_qual = f'{r1_qual[::-1]}{r2_qual[::-1]}'
                else: # i.e. all zero
                    stats.no_bridge += 1
                    continue
                    # Seems not worth even trying if no bridge with a mismatch >=3

                '''

                # Did the read have one bridge on any side?
                print(b1f, b2f, b1r, b2r)

                print(r1_seq)
                print(r2_seq)

                # No obvious bridges. Can we stitch together anyway?
                # try for the very slow path...
                full_qual = f"{r1_qual}{r2_qual}"
                full_seq = needleman_wunsch(r1_seq, r2_seq) # see in NW does a better job;
                print(full_seq)
                full_qual = full_qual[0:len(full_seq)]
                #stats.no_bridge += 1
                #continue
                '''

        if full_seq in pcr_dupes:
            stats.pcr_dupe += 1
            continue
        pcr_dupes.add(full_seq)

        bridge_loc, full_seq, full_qual, fail_type = find_bridge(full_seq, full_qual)

        if not bridge_loc:
            stats.no_bridge += 1
            continue

        stats.with_bridge += 1

        rna = full_seq[0:bridge_loc-1] # still has UMI, -7 cuts of AANNNAA
        dna = full_seq[bridge_loc+17:] # It is definitly 17...

        rna = rna.rstrip("A") # Seems this is useful according the the char-seq code;
        # Do twice for speed;
        if len(dna) <= min_read_size:
            stats.dna_too_short += 1
            continue
        if len(rna) <= min_read_size:
            stats.rna_too_short += 1
            continue

        if dna.startswith('TTTAATTAA'): # No cut, the PacI site is still intact. no genomic DAN on this frag
                          #TTTAATTAAGTCGGAGATCA
            stats.no_bridge_linker_cut += 1
            continue

        #print(rna)

        bridge_seq = full_seq[bridge_loc:bridge_loc+17]
        #print(f'Passed read: {rna} {bridge_seq} {dna}')

        #print(dna[0:30])
        # DNA need to find
        #if dna[0:4] == 'TAAT':
        #    dna = dna[4:]

            #dna = dna[4:]

        if len(dna) < min_read_size:
            stats.dna_too_short += 1
            continue

        rna_reads.append({'name': r1['name'], 'seq': rna, 'qual': full_qual[0:len(rna)]})
        dna_reads.append({'name': r2['name'], 'seq': dna, 'qual': full_qual[len(full_qual)-len(dna):]})
        # if PacI not done, then

    # Make sure these are in the same order as processed above, so they cascade to the failure point!
    print(f'Processed: {idx:,} reads'.format(idx))
    print(f"Both pairs are too short to give a result: {stats.both_pairs_too_short:,} ({stats.both_pairs_too_short / idx *100:.1f}%)")
    print(f' At least one of the reads was a homopolymer: {stats.homopolymer:,} ({stats.homopolymer / idx *100:.1f}%)')
    print(f"  Probable self prime {stats.probable_self_primes:,} ({stats.probable_self_primes / idx *100:.1f}%)")
    print(f"   Can't overlap the reads: {stats.cant_overlap_reads:,} ({stats.cant_overlap_reads / idx *100:.1f}%)")
    print('     PCR dupes: {:,} ({:.1f}%)'.format(stats.pcr_dupe, stats.pcr_dupe / idx * 100))
    print('      Perfect bridge with 0 bp mismatch: {:,} ({:.1f}%)'.format(stats.perfect_bridge, stats.perfect_bridge / idx * 100))
    print('      Rescued bridge with 1-2 bp mismatch: {:,} ({:.1f}%)'.format(stats.rescued_bridge, stats.rescued_bridge / idx * 100))
    print('      Really no bridge: {:,} ({:.1f}%)'.format(stats.no_bridge, stats.no_bridge / idx * 100))
    print('      >1 bridge: {:,} ({:.1f}%)'.format(stats.multi_bridge, stats.multi_bridge / idx * 100))
    print('        With bridge: {:,} ({:.1f}%)'.format(stats.with_bridge, stats.with_bridge / idx * 100))
    #print('          Bridge on forward: {:,} ({:.1f}%)'.format(stats.f_strand, stats.f_strand / stats.with_bridge * 100))
    #print('          Bridge on reverse: {:,} ({:.1f}%)'.format(stats.r_strand, stats.r_strand / stats.with_bridge * 100))
    print(f'           The bridge linker was not cut: {stats.no_bridge_linker_cut:,} ({stats.no_bridge_linker_cut / idx * 100:.1f}%)')
    print(f'           DNA sequence <{min_read_size} bp too short: {stats.dna_too_short:,} ({stats.dna_too_short / idx * 100:.1f}%)')
    print(f'           RNA sequence <{min_read_size} bp too short: {stats.rna_too_short:,} ({stats.rna_too_short / idx * 100:.1f}%)')
    print('Final number of reads kept: {:,} ({:.1f}%)'.format(len(dna_reads), len(dna_reads) / idx * 100))

    file_rna_reads = gzip.open('{}.rna.fq.gz'.format(stub), 'wt')
    for r in rna_reads:
        file_rna_reads.write('{}\n{}\n+\n{}\n'.format(r['name'], r['seq'], r['qual']))
    file_rna_reads.close()

    file_dna_reads = gzip.open('{}.dna.fq.gz'.format(stub), 'wt')
    for d in dna_reads:
        file_dna_reads.write('{}\n{}\n+\n{}\n'.format(d['name'], d['seq'], d['qual']))
    file_dna_reads.close()


