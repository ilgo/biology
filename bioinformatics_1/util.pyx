from itertools import permutations, combinations

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



def motif_generator2(char *kmer, int d):

    cdef int col, car

    yield kmer
    km = tuple(list(kmer))
    for cols in combinations(range(len(km)), d):
        for repl in permutations([65, 67, 71,84], d):
            km3 = list(km)
            for col, car in zip(cols, repl):
                km3[col] = car
            yield ''.join([chr(n) for n in km3])
        

def early_hamming(char *s1, char *s2, int d):

    cdef int count
    cdef char c1, c2

    count = 0
    for c1, c2 in zip(s1, s2):
        if c1 != c2:
            count += 1
            if count > d:
                break
    return count
