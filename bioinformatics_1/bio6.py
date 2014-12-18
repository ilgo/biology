import util
import bio1

from bisect import bisect_left

def nest_list_print(nested_list):
    '''
    print a nested list formatted

    :param nested_list: a list of lists
    :type nested_list: [[`int`, ...], [`int`, ...], ...]

    :rtype: `str`
    '''
    inner = [' '.join(['%+d' %k for k in l]) for l in nested_list]
    return '\n'.join(['(%s)' %i for i in inner])

def greedy_reversal(signed_perm):
    '''
    a greedy reversal sorting

    returns each step of the sorting

    :param signed_perm: a permutation of n signed integers
    :type signed_perm: [`int`, ...]
   
    :rtype: [[`int`, ...], [`int`, ...], ...]
    ''' 
    steps = []
    for n in range(len(signed_perm)):
        idx = signed_perm.index(n+1) if n+1 in signed_perm else signed_perm.index(-(n+1))
        if idx != n:
            rev = [-k for k in signed_perm[n:idx+1]][::-1]
            signed_perm = signed_perm[:n] + rev + signed_perm[idx+1:]
            steps.append(signed_perm[:])
        if signed_perm[n] < 0:
            signed_perm[n] *= -1
            steps.append(signed_perm[:])
    return steps
    

def breakpoint_count(signed_perm):
    '''
    counts the breakpoints in the signed permutation

    returns the count 

    :param signed_perm: a permutation of n signed integers
    :type signed_perm: [`int`, ...]

    :rtype: `int`
    ''' 
    c = 0
    if signed_perm[0] != 1:
        c += 1
    if signed_perm[-1] != len(signed_perm):
        c += 1
    for n,m  in zip(signed_perm[:-1], signed_perm[1:]):
        if n+1 != m:
            c += 1
    return c



def shared_kmers(k, s1, s2):
    '''
    shared kmers can appear on more than one location as well as complements.

    :param `int` k: len of shared kmers
    :param `str` s1: the first string
    :param `str` s2: the second string
    :rtype: [(`int`, `int`), ...]
    :returns: a list of indeces from both strings where kmers match
    '''
    data = {}
    for idx, kmer in enumerate(util.kmer_gen(s2, k)):
        idxs = data.get(kmer, [])
        idxs.append(idx)
        data[kmer] = idxs
    data = [(k,v) for k, v in data.items()]
    data.sort(key=lambda r: r[0])
    keys = [r[0] for r in data]   
    shared = []
    for idx, kmer in enumerate(util.kmer_gen(s1, k)):
        for s in [kmer, bio1.complement(kmer)]:
            pos = bisect_left(keys, s)
            match = data[pos][1] if pos != len(keys) and data[pos][0] == s else []
            if match:
                shared.extend([(idx, m) for m in match])
    return shared

def signed_permutation_string(seq):
    inner = ' '.join(['%+d' %n for n in seq])   
    return '(' + inner + ')'

