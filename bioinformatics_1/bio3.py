import util
import math, sys
from random import randint, random
import operator as op
from itertools import permutations, combinations, product, chain, repeat, accumulate

def motif_enumerate(dnas, k, d):
    '''
    Brute forcing motif enumeration

    Returns the (k,d) motifs in the dna strings

    :param dnas: a list of dna strings
    :type dnas: [`str`, ...]
    :param k: lenght of motifs
    :type k: `int`
    :param d: number of mutations
    :type d: `int`

    :rtype: [`str`, ...]
    '''
    motifs = []
    done = set()
    for dna in dnas:
        for pos in range(len(dna)-k+1):
            kmer = dna[pos:pos+k]
            gen_motifs = [motif for motif in util.motif_generator(kmer, d)]
            for gen_motif in gen_motifs:
                if gen_motif not in done:
                    done.add(gen_motif)
                    gen_motifs2 = [motif for motif in util.motif_generator(gen_motif, d)]
                    if motif_in_dna(gen_motifs2, dnas):
                        motifs.append(gen_motif)
                            
    return motifs


def motif_in_dna(motifs, dnas):
    '''
    determines if at least one of the motifs appears in each dna

    returns True if at least one motif appears in each dna

    :param motifs: a list of motifs of equal length
    :type motifs: [`str`, ...]
    :param dnas: a list of dna strings
    :type dnas: [`str`, ...]

    :rtype: `bool`
    '''
    for dna in dnas:
        contained = False
        for motif in motifs:
            if motif in dna:
                contained = True
                break

        if not contained:
            return False
    return True


def entropy(a, c , g , t):
    '''
    calculate the information entropy

    :param a,c,g,t: the 4 nucleotides
    :type a,c,g,t: `float` [0.0, 1.0]

    :rtype: `float`
    '''
    en_a = a*math.log2(a) if a != 0.0 else 0.0
    en_c = c*math.log2(c) if c != 0.0 else 0.0
    en_g = g*math.log2(g) if g != 0.0 else 0.0
    en_t = t*math.log2(t) if t != 0.0 else 0.0
    return -(en_a + en_c + en_g + en_t)



def motif_ham_score(pattern, dnas):
    '''

    '''
    scores = [] 
    k = len(pattern)
    for dna in dnas:
        scores.append(min([util.hamming(pattern, s) for s in util.kmer_gen(dna, k)]))
    return sum(scores)   


def median_string(dna, k):
    '''
    find the k-mer in dna which is the median string
    
    returns a string of len k

    :param dna: a list of dna strings
    :type dna: [`str`, ...]
    :param k: the len of the pattern to discover
    :type k: `int`

    :rtype: `str`
    '''
    best_pat = ''
    best_score = sys.maxsize
    for pattern in product(*['ACGT'] * k):
        pattern = ''.join(pattern)
        score = motif_ham_score(pattern, dna)
        if score <= best_score:
            best_pat = pattern
            best_score = score
    return best_pat


def most_probable(dna, k, profile):
    '''
    Find the kmer most probable string in dna, from a nucleotide profile

    :param dna: a dna string
    :type dna: `str`
    :param k: the kmer length
    :type k: `int`
    :param profile: a profile matrix '4 x k'
    :type profile: [[`float`, ...], [`float`, ...], [`float`, ...], [`float`, ...]]

    :rtype: `str`
    '''
    #print(k, len(profile[0]))
    assert k == len(profile[0])
    n_pos = {k:v for k,v in zip('ACGT', range(4))}
    best_kmer = ''
    best_score = -1.0
    for kmer in util.kmer_gen(dna, k):
        score = 1 
        for i, c in enumerate(kmer):
            val = profile[n_pos[c]][i]
            score *= val if val != 0 else 1.0
        if score > best_score:
            best_score = score
            best_kmer = kmer
    return best_kmer


def motif_matrix(dnas, pseudo=False):
    '''
    calculate a profile from dnas

    returns a matrix 4 x k, the 4 rows are the A, C, G and T rows.

    :param dnas: a list of kmers
    :type dnas: [`str`, ...]

    :rtype: [[`float`, ...], [`float`, ...], [`float`, ...], [`float`, ...]]
    '''
    div = len(dnas)+4 if pseudo else len(dnas)
    init = 1 if pseudo else 0

    k = len(dnas[0])
    profile = [[0]*k for _ in range(4)]

    for n in range(k):
        counts = {k:v for k,v in zip('ACGT', repeat(init))}
        for row in dnas:
            counts[row[n]] += 1
        counts = {k:v/div for k,v in counts.items()}
        for i, c in enumerate('ACGT'):
            profile[i][n] = counts[c]
    return profile
        

def consensus(motifs):
    '''
    return the consensus for the motifs

    returns a kmer representing the consensus string 

    :param motifs: a list of kmers
    :type motifs: [`str`, ...]

    :rtype: `str`
    '''
    consens = []
    for i in range(len(motifs[0])):
        count = {k:v for k,v in zip('ACGT', [0]*4)}
        for j in range(len(motifs)):
            c = motifs[j][i]
            count[c] += 1
        c1 = 'A' if count['A'] > count['C'] else 'C'
        c2 = 'G' if count['G'] > count['T'] else 'T'
        c = c1 if count[c1] > count[c2] else c2
        consens.append(c)
    return ''.join(consens)


def greedy_motif(k, t, dnas, pseudo=False):
    '''

    '''
    best_motifs = [dna[:k] for dna in dnas]
    best_score = sys.maxsize
    motifs = best_motifs[:]

    for kmer in util.kmer_gen(dnas[0], k):
        motifs = [kmer]
        for i in range(1, t):
            profile = motif_matrix(motifs, pseudo=pseudo)
            probable = most_probable(dnas[i], k, profile)
            motifs.append(probable)
        score = motif_ham_score(consensus(motifs), motifs)
        if score <= best_score:
            best_score = score
            best_motifs = motifs[:]
    return best_motifs


def _randomized_motif_search(dna, k, t, pseudo=False):
    '''

    :param dna: a list of dna strings
    :type dna: [`str`, ...]
    :param k: the len of the pattern to discover
    :type k: `int`
    :param t: how many motifs are there
    :type t: `int`

    :rtype: [`str`, ...]
    '''
    motifs = []
    for i in range(t):
        pos = randint(0, len(dna[0])-k)
        motifs.append(dna[i][pos:pos+k])
    best_motifs = motifs[:]
    best_score = motif_ham_score(consensus(best_motifs), best_motifs)

    while True:
        profile = motif_matrix(best_motifs, pseudo=pseudo)
        motifs = [most_probable(d, k, profile) for d in dna] 
        score = motif_ham_score(consensus(motifs), motifs)
        if score < best_score:
            best_motifs = motifs[:]
            best_score = score
        else: 
            break
    return (best_motifs, best_score)

def randomized_motif_search(dna, k, t, N):
    
    best_motifs = None
    best_score = sys.maxsize
    for n in range(N):
        motifs, score = _randomized_motif_search(dna, k, t, pseudo=True)
        if score < best_score:
            best_motifs = motifs
            best_score = score
    return best_motifs, best_score



def _gibbs_sampler(dna, k, t, N):
    motifs = []
    for i in range(t):
        pos = randint(0, len(dna[0])-k)
        motifs.append(dna[i][pos:pos+k])
    best_motifs = motifs[:]
    best_score = motif_ham_score(consensus(best_motifs), best_motifs)

    for _ in range(N):
        i = randint(0, t-1)
        gibbs_motifs = motifs[:i] + motifs[i+1:]
        profile = motif_matrix(gibbs_motifs, pseudo=True)
        motifs[i] = profile_random_motif(dna[i], k, profile)        
        score = motif_ham_score(consensus(motifs), motifs)
        if score < best_score:
            best_motifs = motifs[:]
            best_score = score
    return (best_motifs, best_score)
    

def profile_random_motif(dna, k, profile):
    motifs = {}
    index = {'A':0, 'C':1, 'G':2, 'T':3}
    for kmer in util.kmer_gen(dna, k):
        values = [profile[index[c]][i] for i,c in enumerate(kmer)]
        motifs[kmer] = [val for val in accumulate(values, func=op.mul)][-1]
    denom = sum([val for val in motifs.values()])
    dice = []
    val = 0.0
    for k,v in motifs.items():
        val += v/denom
        dice.append((val, k))
    choice = random()
    for v,k in dice:
        if choice < v:
            return k
    assert 0==1
 


def gibbs_sampler(dna, k, t, N):
    best_motifs = None
    best_score = sys.maxsize
    for _ in range(20):
        motifs, score = _gibbs_sampler(dna, k, t, N)
        if score < best_score:
            best_motifs = motifs
            best_score = score
    return best_motifs, best_score
    
