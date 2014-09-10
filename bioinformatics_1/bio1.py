from collections import defaultdict
from kmerset import PatCount, KmerSet

COMPLEMENT = {'A':'T', 'T':'A', 'C':'G', 'G':'C'} 

def pp(result):
    '''
    get the result as a one-liner w/o commas

    returns the list as a string

    :param result: a list of objects
    :type result: [ obj, ...]

    :rtype: `str`
    '''
    return ' '.join([str(obj) for obj in result])


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
    rev_range = range(len(text)-1, -1, -1)
    compl = [COMPLEMENT.get(text[n]) for n in rev_range]
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

def approximate_kmer(genome, kmer, min_count=0, mismatch=0):
    '''
    rosalind_1a.
    find the k-mers of given size
    
    returns a list of the most frequent k-mers of lenght 'size'

    :param genome: the DNA string to search
    :type genome: `str`
    :param kmer: the lenght of the k-mers to discover
    :type kmer: `int`
    :param min_count: the kmer must appear at least this many times
    :type min_count: `int`
    :param mismatch: the allowable mismatches
    :type mismatch: `int`

    :rtype: [`str`, ...]
    '''
    #kmersets = []
    #
    #for n in range(len(genome)-kmer-1):
    #    done = False
    #    pc = PatCount(genome[n:n+kmer])
    #    for kmerset in kmersets:
    #        if pc in kmerset:
    #            kmerset.insert(pc)
    #            done = True
    #            break
    #    if not done:
    #        for kmerset in kmersets:
    #            if is_approximate(pc.pattern, kmerset.head().pattern, mismatch):
    #                kmerset.insert(pc)
    #                done = True
    #                break
    #    if not done:
    #        kmerset = KmerSet()
    #        kmerset.insert(pc)
    #        kmersets.append(kmerset) 

    #print('\n'.join(str(kms) for kms in kmersets))
    #max_len = max([kmerset.head().count for kmerset in kmersets]) if min_count == 0 else min_count
    #return [kmerset.head().pattern for kmerset in kmersets if kmerset.head().count == max_len]
         

    table = {}
    max_len = 0
    for n in range(len(genome)-kmer-1):
        pattern = genome[n:n+kmer]
        v = 0 
        for m in range(len(genome)-kmer-1):
            p2 = genome[m:m+kmer]
            v += is_approximate(pattern, p2, mismatch)
        if v >= max_len:
            max_len = v
            table[pattern] = v

            
    #print(table)
    #max_len = max(table.values()) if min_count == 0 else min_count
    return [k for k, v in table.items() if v >= max_len]
