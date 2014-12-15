from itertools import product, combinations
from random import randint
import bio2


def delta_A(seq):
    '''
    create the pairwise differences for a given integer sequence
    or create data for the turnpike problem.

    The following equality must hold.
        seq == turnpike(delta_A(seq))    

    :param [`int`, ...] seq: a sequence of integers
    :returns: the pairwise differences of all elements in the sequence
    :rtype: [`int, ...]
    '''
    res = [n-m for n, m in product(seq, repeat=2)]
    return sorted(res)


def turnpike_1(seq):
    '''
    first try, does not work
    '''
    points = [0]
    diffs = [n for n in sorted(seq) if n > 0]
    for k in range(1,seq.count(0)-1):
        end = diffs[-1]
        i = 0
        j = len(diffs) -1
    
        while True:
            if i >= j:
                print('!!!!! overflow')
                break

            diff_sum = diffs[i] + diffs[j]
            if diff_sum == end:
                break
            elif diff_sum > end:
                j -= 1
            else:
                i += 1
        if i < j:
            print(k, i, diffs[i], j, diffs[j])
            points.append(diffs[i] + points[k-1])
            diffs.pop(i)
            diffs.pop(j)
    points.append(max(seq))
    return points
           


def rand_seq(n, m):
    '''
    a sequence with n+2 values taken from 1..m
    return [0, n1, n2, ..., m]
    ''' 
    seq = []
    for i in range(n):
        r = randint(1, m)
        while r in seq:
            r = randint(1, m)
        seq.append(r)    
    seq.append(0)
    seq.append(m)
    return sorted(seq)


def turnpike_2(seq):
    sols = [[0, max(seq)]]
    diffs = [n for n in sorted(seq) if n > 0]
    seq_spectrum = bio2.spectrum_dict([n for n in sorted(seq) if n >= 0])

    for n,m in pairwise(diffs):
        new_sols = []
        for k in set([n,m]):
            for sol in sols:
                new_sol = sol[:]
                new_sol.append(k)
                new_sols.append(sorted(new_sol)) 
        sols = []
        for sol in new_sols:
            delta = [n for n in delta_A(sol) if n >= 0]
            sol_spectrum = bio2.spectrum_dict(delta)
            if bio2.is_inner_spectrum(sol_spectrum, seq_spectrum):
                sols.append(sol)
        
    return sols

def pairwise(diffs):
    length = max(diffs)
    for n in set(diffs):
        if n <= (length+1)/2:
            m = length - n
            if m != n and m in diffs:
                yield (n,m)
            elif diffs.count(m) >= 2:
                yield(n, m)
            
# next try. create all possible pairwise solution, so that there will be an 
# exact amount of slots to be filled.
# possible slots come from the pairwaise summands

def turnpike(seq):
    zeros = seq.count(0)
    print(zeros)
    top = max(seq)
    diffs = [n for n in sorted(seq) if n > 0]
    pair_diff = [pd for pd in pairwise(diffs)]
    for prod in product(*pair_diff):
        print(len(prod))
        for slots in combinations(prod, zeros - 2):
            pike = [0, top]
            pike.extend(slots)
            if delta_A(pike) == seq:
                return sorted(pike)
 
