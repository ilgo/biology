from collections import defaultdict
from datetime import datetime
from kmerset import PatCount, KmerSet
import util

TABLE = {'A':'0', 'C':'1', 'G':'2', 'T':'3'}
COMPLEMENT = {'A':'T', 'T':'A', 'C':'G', 'G':'C'} 

def pp(result, join_char=' '):
    '''
    get the result as a one-liner w/o commas

    returns the list as a string

    :param result: a list of objects
    :type result: [ obj, ...]
    :param join_char: the character to use for joining the strings
    :type join_char: `str` of length 1

    :rtype: `str`
    '''
    return join_char.join([str(obj) for obj in result])


def kmer(text, size, min_count=0):
    '''
    rosalind_1a.
    find the k-mers of given size
    
    returns a list of the most frequent k-mers of lenght 'size'

    :param text: the DNA string to search
    :type text: `str`
    :param size: the lenght of the k-mers to discover
    :type size: `int`
    :param min_count: the kmer must appear at least this many times
    :type min_count: `int`

    :rtype: [`str`, ...]
    '''
    table = defaultdict(int)
    for n in range(len(text)-size-1):
        table[text[n:n+size]] += 1
    max_len = max(table.values()) if min_count == 0 else min_count
    return [k for k, v in table.items() if v >= max_len]


def complement(text):
    '''
    rosalind_1b.
    transform this text into its DNA complement {A <-> T, C <-> G}

    return a complementary string

    :param text: a DNA string
    :type text: `str`

    :rtype: `str`
    '''
    compl = [COMPLEMENT.get(c) for c in text[::-1]]
    return ''.join(compl)


def pattern_start(genome, pattern, mismatch=0):
    '''
    rosalind_1c.
    find all starting positions of the pattern in the genome, with at most d mismatches

    returns a list of starting positions

    :param genome: the genome that contains the pattern
    :type genome: `str`
    :param pattern: the pattern for matching
    :type pattern: `str`
    :param mismatch: the allowable mismatches
    :type mismatch: `int`

    :rtype: [`int`, ...]
    '''
    pat_len = len(pattern)
    size = len(genome) - pat_len + 1
    return [i for i in range(size) if is_approximate(genome[i:i+pat_len], pattern, mismatch)]

def clumps(genome, k, clump, count):
    '''
    rosalind_1d.
    find clumps of kmers inside the genome. A clump is a pattern of given lenght that appears at least count times.

    returns a list of kmers that fullfill the clump condition

    :param genome: the genome to analyze
    :type genome: `str`
    :param k: the lenght of the kmer that are clumping.
    :type k: `int`
    :param clump: the length of the clump
    :type clump: `int`
    :param count: the minmimum appearance of the kmer in the clump
    :type count: `int`

    :rtype: [ `str`, ...]
    '''
    kmers = set() 
    clmps = (genome[i:i+clump] for i in range(len(genome) - clump + 1))
    for km in (kmer(clmp, k, min_count=count) for clmp in clmps):
        kmers.update(km)
    return list(kmers)

def skew(genome):
    '''
    rosalind_1e.
    skew is the difference between C & G bases in a genome.

    returns a list of prefix_i(genome) C/G difference counts

    :param genome: calculate the skew from this
    :type genome: `str`

    :rtype: [ `int`, ...]
    '''
    sk = 0
    skews = []
    skews.append(sk)
    for base in genome:
        if base == 'C':
            sk -= 1
        elif base == 'G':
            sk += 1
        skews.append(sk)   
    return skews

def min_skew(genome):
    '''
    rosalind_1e.
    fird the minimum skews in the genome

    returns a list of indeces from th egeome where skew is a minimum

    :param genome: the genome string
    :type genome: `str`
    
    :rtype: [`int`, ...]
    '''
    _skew = skew(genome)
    _skew_min = min(_skew)
    return [i for i, sk in enumerate(_skew) if sk == _skew_min]


def is_approximate(genome, pattern, mismatch):
    '''
    rosalind_1.f
    find all approximate matches of pattern in genome with a maximum of 'mismatch' count

    returns True if the genome matches the pattern with at most d mismatches

    :param genome: the genome that contains the pattern
    :type genome: `str`
    :param pattern: the pattern for matching
    :type pattern: `str`
    :param mismatch: the allowable mismatches
    :type mismatch: `int`

    :rtype: `bool`
    '''
    misses = 0
    for g, p in zip(genome, pattern):
        if g != p:
            misses += 1
            if misses > mismatch:
                return False
    return True

#-----------------------------------------
#stepic

def pat2num(pat):
    '''
    transform pattern to number

    :param pat: a ACGT dna pattern
    :type pat: `str`
    
    :rtype: `int`
    '''
    return int(''.join([TABLE[c] for c in pat]) , 4)

def num2pat(n, k):
    '''
    return the dna represented by the number for kmer k

    :param n: the pattern number
    :type n: `int`
    :param k: the number belongs to kmer k
    :type k: `int`

    :rtype: `str`
    '''
    nucleotides = []
    for i in range(k-1, -1, -1):
        n, r = divmod(n, 4**i)
        nucleotides.append('ACGT'[n])
        n = r
    return ''.join(nucleotides) 

def frequency_array(dna, k):
    '''
    generate a frequency array for the kmers in the dna string

    returns a list of counts for all kmers, non exisitng kmers must be filled with 0's

    :param str dna: the string to examine
    :param int k: the kmer length
    :rtype: [`int`, ...]
    '''
    freqs = {num2pat(n, k):0 for n in range(4**k)}
    for s in util.kmer_gen(dna, k):
        freqs[s] += 1
    return [v for k,v in sorted(freqs.items(), key=lambda t: t[0])] 

def pattern_count(text, pattern, d=0):
    return sum([util.hamming(kmer, pattern) <= d for kmer in util.kmer_gen(text, len(pattern))])


def approximate_pattern(text, pattern, d):
    '''
    get the starting positions of all patterns with at most d mismatches in the text

    :param str text: the genome to be analyzed
    :param str pattern: the pattern to be detected
    :param int d: mismatches allowed
    :return: a list of starting positions
    :rtype: [`int`, ...]
    '''
    pos = []
    for i, kmer in enumerate(util.kmer_gen(text, len(pattern))):
        if util.hamming(pattern, kmer) <= d:
            pos.append(i)
    return pos


def early_hamming(s1, s2, d):
    count = 0
    for c1, c2 in zip(s1, s2):
        if c1 != c2:
            count += 1
            if count > d:
                break
    return count


def mismatch_kmers(text, k, d, compl=False):
    '''
    find all most frequent kmers in text that have at most d mismatches
    
    seems to be a variation of lesson3 motif search

    :param str text: the genome to be analyzed
    :param int k: patterns of lenght k 
    :param int d: mismatches allowed
    :return: most frequent kmer with max d mismatches
    :rtype: [`str`, ...]
    '''
    kmers = defaultdict(int)
    for kmer in util.kmer_gen(text, k):
        kmers[kmer] += 1

    motifs = defaultdict(int)
    for kmer, freq in kmers.items():
        for motif in util.motif_generator(kmer, d):
            motifs[motif] += freq
        
        if compl:
            for motif in util.motif_generator(complement(kmer), d):
                motifs[motif] += freq

    max_score = max(motifs.values())
    return sorted([kmer for kmer, freq in motifs.items() if freq == max_score])

