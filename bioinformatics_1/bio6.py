

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
