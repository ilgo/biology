from itertools import permutations, combinations, product

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


def kmer_gen(dna, k):
    '''
    return all kmers of length k form dna string

    returns a string generator:

    :param dna: a dna string
    :type dna: `str`
    :param k: lenght of kmers to be returned
    :type k: `int`

    :rtype: (`str`, ...)
    '''
    for pos in range(len(dna)-k+1):
        yield dna[pos:pos+k]


def hamming(s1, s2):
    '''
    calculate the Hamming distance between string s1 & s2

    returns the distance

    :param s1: the first string
    :type s1: `str`
    :param s2: the second string
    :type s2: `str`

    :rtype: `int`
    '''
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))



def motif_generator(kmer, d):

    yield kmer
    alt_bases = {'A':'CGT', 'C':'AGT', 'G':'ACT', 'T':'ACG'}
    for dd in range(1, d+1): 
        for cols in combinations(range(len(kmer)), dd):
            for repl in product(*[alt_bases[kmer[i]] for i in cols]):
                km = list(kmer)
                for col, char in zip(cols, repl):
                    km[col] = char
                yield ''.join(km)
        

