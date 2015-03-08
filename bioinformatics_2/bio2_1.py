import util

def trie(patterns):
    '''
    create an adjacency list of the trie for the given patterns
    each node = {i: [[chars], [pos]]}
    returns tuples with (start, end, character)

    :param patterns: a list of strings
    :type patterns: [`str`, ...]
    :rtype: [(`int`, `int`, `char`), ...) 
    '''
    nodes = {1:[[],[]]}
    count = 1
    for s in patterns:
        idx = 1
        for c in s:
            if c in nodes[idx][0]:
                pos = nodes[idx][0].index(c)
                idx = nodes[idx][1][pos]
            else:
                count += 1
                nodes[idx][0].append(c)
                nodes[idx][1].append(count)
                nodes[count] = [[],[]]
                idx = count
    return nodes

def trie_repr(trie):
    '''
    Rosalind wants the trie in a special format
    '''
    edges = []
    for k, vs in trie.items():
        for c, t in zip(vs[0], vs[1]):
            edges.append('%d %d %s' %(k, t, c))
    return edges 


def trie_matching(dna, patterns):
    '''
    find all starting points in dna, where a pattern is found
    returns a list of the startig positions.

    :param `str` dna: here the patterns are to be discovered
    :param patterns: a set of strings
    :rtype patterns: [`str`, ...]
    '''
    matches = []
    tree = trie(patterns)
    for i in range(len(dna)):
        node = tree[1]
        pos = i
        while dna[pos] in node[0]:
            idx = node[0].index(dna[pos])
            node = tree[node[1][idx]]
            pos += 1
            if pos >= len(dna):
                break
        if not node[0]:
            matches.append(i)
    return matches


def longest_repeat(s):
    '''
    longest repated sequence in a string
    '''
    best = ''
    subfix_sort = sorted([s[i:] for i in range(len(s))])
    for lex_less, lex_more in zip(subfix_sort[:-2], subfix_sort[1:]):
        if len(lex_less) > len(lex_more):
            lex_less, lex_more = lex_more, lex_less
        pos = 0
        while pos < len(lex_less) and lex_less[pos] == lex_more[pos]:
            pos += 1
        if pos > len(best):
            best = lex_less[:pos]
            #print('>>> ', best) 
    return best


def bwt(s):
    '''
    burrow-whellers transform of a string
    '''
    ss = s + s[:-1]
    bwt_mat = sorted([p for p in util.kmer_gen(ss, len(s))])
    #print (bwt_mat)
    return ''.join([p[-1] for p in bwt_mat])  


def tag_tuple(s):
    '''
    tag each entry with is character occurrence count
    '''
    tt = []
    cs = {'A':0, 'C':0, 'G':0, 'T':0, '$':0}
    for c in s:
        cs[c] += 1
        tt.append((c, cs[c]))
    return tt


def inv_bwt(inv_s):
    '''
    Reconstruct a string from its BWT transform
    '''
    s = []
    # create columns with char index
    last = tag_tuple(inv_s)
    first = tag_tuple(sorted(inv_s))

    #print(last)
    #print(first)

    # loop back to beginning
    l_pos = last.index(('$',1))
    while l_pos != 0:
        entry = first[l_pos]
        s.append(entry[0])
        l_pos = last.index(entry)
    s.append('$')
    return ''.join(s)

def bw_match(bwt, patterns):
    '''

    '''
    res = []
    last = tag_tuple(bwt)
    first = tag_tuple(sorted(bwt)) 

    for p in patterns:
        pos = first.index((p[0], 1))
        while(first[pos][0] == p[0]):
            c = 0
            pos = last.index(first[pos])
            
