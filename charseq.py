
import sys, os, gzip
import regex

compdict = {'A': 'T',
            'C': 'G',
            'G': 'C',
            'T': 'A',
            'N': 'N',}

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

if len(sys.argv) < 3:
    print('simple_split.py: PE1.fastq.gz PE2.fastq.gz')
    sys.exit()

stub = os.path.split(sys.argv[1])[1].replace('.fastq', '').replace('.fq', '').replace('.gz', '').replace('.trim', '').replace('.p1','')

rna_reads = []
dna_reads = []
pcr_dupes = set([])

min_read_size = 20

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

bridge_f = 'ACCGGCGTCCAAG' # for searching, actual bridge is a bit longer;
bridge_r = 'CTTGGACGCCGGT'

def find_bridge(seq, qual, stats=stats, mismatch=2):
    find_f = regex.findall("(ACCGGCGTCCAAG)", seq)
    find_r = regex.findall("(CTTGGACGCCGGT)", seq)

    if len(find_f) > 1 or len(find_r) > 1:
        stats.multi_bridge += 1
        return None, None, None

    elif len(find_f) == 1 and len(find_r) == 1:
        stats.multi_bridge += 1
        return None, None, None

    elif len(find_f) > 0:
        stats.perfect_bridge += 1
        stats.f_strand += 1
        newseq = full_seq
        newqual = qual
        bridge_loc = newseq.find(bridge_f)

    elif len(find_r) > 0:
        stats.perfect_bridge += 1
        stats.r_strand += 1
        newseq = rc(full_seq)
        newqual = qual[::-1]
        bridge_loc = newseq.find(bridge_f) # Get the new F strand location;

    else:
        # slower mismatch path:
        hit_f = regex.findall("(ACCGGCGTCCAAG){s<=3}", seq)
        hit_r = regex.findall("(CTTGGACGCCGGT){s<=3}", seq)

        if len(hit_f) > 1 or len(hit_r) > 1:
            stats.multi_bridge += 1
            return None, None, None
        elif len(hit_f) >= 1 and len(hit_r) >= 1:
            stats.multi_bridge += 1
            return None, None, None
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
            return None, None, None

    if bridge_loc == -1:
        # Should be impossible to reach here
        1/0

    return bridge_loc, newseq, newqual

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

#print(commonOverlapIndexOf_mismatch('ACTCGTTGCT',  'ACTCGTTGCTTTACAA'))
#print(commonOverlapIndexOf_mismatch('CACTCCTTGCT',  'ACTCGTTGCTTTACAA'))
# CACTCCTTGCT
#  ACTCGTTGCTTTACAA [expected]
#   ACTCGTTGCT
#1/0

idx = 0
for r1, r2 in fastqPE(sys.argv[1], sys.argv[2]):
    idx += 1
    if (idx) % 1e5 == 0:
        print('{:,}'.format(idx))
        #break

    if len(r1['seq']) + len(r2['seq']) < 55: # No way to get info out of this small
        stats.both_pairs_too_short += 1
        continue

    # Can reject the simple case when one read is a homopolymer
    if len(set(r1['seq'])) == 1 or len(set(r2['seq'])) == 1:
        stats.homopolymer += 1
        continue

    # This is wrong, as reads may go into each other.
    # I can't just try to overlap as they may have mismatches;
    # Also the r2 is the RC anyway...
    #full_seq = r1['seq'].upper() + r2['seq'].upper()
    #qual_seq = r1['qual'].upper() + r2['qual'].upper()

    r1_seq = r1['seq']
    r2_seq = rc(r2['seq'])
    r1_qual = r1['qual']
    r2_qual = r2['qual'][::-1]
    kmpa = commonOverlapIndexOf(r1_seq, r2_seq)
    kmpb = commonOverlapIndexOf(r2_seq, r1_seq) # Need to do twice for convergent;

    # get the best overlap
    if kmpa < 3 and kmpb < 3: # no great overlap found, so assume they can be concatenated;
        # seems many of the errors are on the 3' end, trim it and see if you get an overlap;
        r1_seq = r1_seq[5:-5]
        r1_qual = r1_qual[5:-5]
        r2_seq = r2_seq[5:-5]
        r2_qual = r2_qual[5:-5]
        kmpa = commonOverlapIndexOf(r1_seq, r2_seq)
        kmpb = commonOverlapIndexOf(r2_seq, r1_seq) # Need to do twice for convergent;
        #print('R1', r1_seq)
        #print('R2', r2_seq)
        #print(kmpa, kmpb)
        if kmpa < 3 and kmpb < 3:
            # I think there is a reasonable argument for pasting them together
            # and just letting the bridge finder see if it looks valid.
            full_seq = f'{r1_seq}{r2_seq}'
            full_qual = f"{r1_qual}{r2_qual}"
            kmpa = -1
            kmpb = -2
            #kmpa = commonOverlapIndexOf(r1_seq, r2_seq)
            #kmpb = commonOverlapIndexOf(r2_seq, r1_seq) # Need to do twice for convergent;
            # Can't find an overlap even after trimming
            # Try a mismatch version here, about ~28%. This part is very slow
            #kmpa = commonOverlapIndexOf_mismatch(r1_seq, r2_seq)
            #kmpb = commonOverlapIndexOf_mismatch(r2_seq, r1_seq) # Need to do twice for convergent;
            #print('R1', r1_seq)
            #print('R2', r2_seq)
            #print(kmpa, kmpb)
            #print(' O', full_seq)
            #stats.cant_overlap_reads += 1
            #continue

    #print('R1', r1_seq)
    #print('R2', r2_seq)
    #print('KMPA', kmpa, 'KMPB', kmpb)

    if kmpa == kmpb:
        # seems these are mainly palindromes, self amplifications?
        stats.probable_self_primes += 1
        continue
    if kmpa < 3 and kmpb < 3:
        # Can't find an overlap even after trimming and mismatch (above?)
        #print('R1', r1_seq)
        #print('R2', r2_seq)
        #print(' O', full_seq)
        #stats.cant_overlap_reads += 1
        #continue
        full_seq = f'{r1_seq}{r2_seq}'
        full_qual = f"{r1_qual}{r2_qual}"
    elif kmpa > kmpb: # correct;
        full_seq = f'{r1_seq}{r2_seq[kmpa:]}'
        full_qual = f"{r1_qual}{r2_qual[kmpa:]}"
    elif kmpb > kmpa: # correct;
        full_seq = f'{r2_seq}{r1_seq[kmpb:]}'
        full_qual = f'{r2_qual}{r1_qual[kmpb:]}'
    else:
        1/0 # Impossible! (Maybe)

    #print(' O', full_seq)

    if full_seq in pcr_dupes:
        stats.pcr_dupe += 1
        continue
    pcr_dupes.add(full_seq)

    bridge_loc, full_seq, full_qual = find_bridge(full_seq, full_qual)
    if not bridge_loc:
        stats.no_bridge += 1
        continue

    stats.with_bridge += 1

    rna = full_seq[0:bridge_loc-7] # still has UMI, -7 cuts of AANNNAA
    dna = full_seq[bridge_loc+17:]

    rna = rna.rstrip("A")
    # Do twice for speed;
    if len(dna) < min_read_size:
        stats.dna_too_short += 1
        continue
    if len(rna) < min_read_size:
        stats.rna_too_short += 1
        continue

    if dna.startswith('TTTAATTAAGTCG'): # No cut, no genomic DAN on this frag
                       #TTTAATTAAGTCGGAGATCA
        stats.no_bridge_linker_cut += 1
        continue

    bridge_seq = full_seq[bridge_loc:bridge_loc+17]
    #print(f'read: {rna} {bridge_seq} {dna}')
    #print(f'{rna} {bridge_seq}')
    #print(f'{bridge_seq} {dna}')

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
print('      [ Perfect bridge with 0 bp mismatch: {:,} ({:.1f}%)'.format(stats.perfect_bridge, stats.perfect_bridge / idx * 100))
print('      [ Rescued bridge with 1-2 bp mismatch: {:,} ({:.1f}%)'.format(stats.rescued_bridge, stats.rescued_bridge / idx * 100))
print('      [ Really no bridge: {:,} ({:.1f}%)'.format(stats.no_bridge, stats.no_bridge / idx * 100))
print('      [ >1 bridge: {:,} ({:.1f}%)'.format(stats.multi_bridge, stats.multi_bridge / idx * 100))
print('        With bridge: {:,} ({:.1f}%)'.format(stats.with_bridge, stats.with_bridge / idx * 100))
print('          [ Bridge on forward: {:,} ({:.1f}%)'.format(stats.f_strand, stats.f_strand / stats.with_bridge * 100))
print('          [ Bridge on reverse: {:,} ({:.1f}%)'.format(stats.r_strand, stats.r_strand / stats.with_bridge * 100))
print(f'           [ The bridge linker was not cut: {stats.no_bridge_linker_cut:,} ({stats.no_bridge_linker_cut / idx * 100:.1f}%)')
print(f'           [ DNA sequence <{min_read_size} bp too short: {stats.dna_too_short:,} ({stats.dna_too_short / idx * 100:.1f}%)')
print(f'           [ RNA sequence <{min_read_size} bp too short: {stats.rna_too_short:,} ({stats.rna_too_short / idx * 100:.1f}%)')
print('Final number of reads kept: {:,} ({:.1f}%)'.format(len(dna_reads), len(dna_reads) / idx * 100))

file_rna_reads = gzip.open('{}.rna.fq.gz'.format(stub), 'wt')
for r in rna_reads:
    file_rna_reads.write('{}\n{}\n+\n{}\n'.format(r['name'], r['seq'], r['qual']))
file_rna_reads.close()

file_dna_reads = gzip.open('{}.dna.fq.gz'.format(stub), 'wt')
for d in dna_reads:
    file_dna_reads.write('{}\n{}\n+\n{}\n'.format(d['name'], d['seq'], d['qual']))
file_dna_reads.close()


