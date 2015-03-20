from collections import namedtuple
from pprint import pprint

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
    nodes = {0:[[],[]]}
    count = 0
    for s in patterns:
        idx = 0
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
            edges.append('%d->%d:%s' %(k, t, c))
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
        node = tree[0]
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
          

Node = namedtuple('Node', ['value', 'nodes', 'leaves'])
Leaf = namedtuple('Leaf', 'pos sid')
#class Leaf(namedtuple('Leaf', 'pos sid, suf')):
#    __slots__ = ()
#    def __str__(self):
#        return self.suf


class STree:
    def __init__(self, patterns):
        self.root = Node(value=[], nodes=[], leaves=[])
        self.pat_map = self.make_pattern_mapping(patterns)
        self.construct()

    def make_pattern_mapping(self, patterns):
        '''
        make an id:pattern+$ mapping
        return the created dictionary
        ''' 
        pat_map = {}
        if not isinstance(patterns, list):
            patterns = [patterns]
        for pat in patterns:
            if not pat.endswith('$'): 
                pat = pat + '$'
            pat_map[id(pat)] = pat
        return pat_map

    def construct(self):
        for sid, pat in self.pat_map.items():
            for pos, suffix in enumerate(util.suffixes(pat)):
                #print (self.edges())
                #print (pos, suffix)
                node = self.root
                inserted = False
                while not inserted:
                    prefix_leaf, prefix = self.same_prefix_leaf(suffix, node.leaves)
                    if prefix_leaf:
                        leaf1 = Leaf(pos=pos+prefix, sid=sid)
                        leaf2 = Leaf(pos=prefix_leaf.pos+prefix, sid=prefix_leaf.sid)
                        node.leaves.remove(prefix_leaf)
                        subnode = Node(value=[suffix[:prefix]], nodes=[], leaves=[leaf1, leaf2])
                        node.nodes.append(subnode)
                        inserted = True
                        #print ('New Node: ', subnode.value, pat[leaf1.pos:], pat[leaf2.pos:])
                        continue                       
                    prefix_node, prefix = self.same_prefix_node(suffix, node.nodes)
                    if prefix_node:
                        if prefix != len(prefix_node.value[0]):
                            #print (prefix_node, suffix)
                            leaf = Leaf(pos=pos+prefix, sid=sid)
                            prefix_val = prefix_node.value[0]
                            subnode = Node(value=[prefix_val[prefix:]], nodes=prefix_node.nodes[:], leaves=prefix_node.leaves[:])
                            prefix_node.value[0] = prefix_val[:prefix]
                            prefix_node.nodes.clear()
                            prefix_node.leaves.clear()
                            prefix_node.nodes.append(subnode)
                            prefix_node.leaves.append(leaf)
                            inserted = True
                            #print ('New Node: ', subnode.value[0], pat[leaf.pos:])
                        else:
                            node = prefix_node
                            pos += prefix
                            suffix = suffix[prefix:]
                    else:
                        # prefix is 0
                        leaf = Leaf(pos=pos, sid=sid)
                        node.leaves.append(leaf)
                        inserted = True     

    def same_prefix_node(self, suffix, nodes):
        # check all suffixes of node.value and get longest common one.
        node, prefix = None, 0
        for n in nodes:
            limit = min(len(suffix), len(n.value[0]))
            c = 0
            while c < limit and suffix[c] == n.value[0][c]:
                c += 1
            if c > prefix:
                node = n
                prefix = c
        return node, prefix
            

    def same_prefix_leaf(self, suffix, leaves):
        leaf, prefix = None, 0       
        for l in leaves:
            pat = self.pat_map[l.sid][l.pos:]
            limit = min(len(suffix), len(pat))
            c = 0
            while c < limit and pat[c] == suffix[c]:
                c+= 1
            if c > prefix:
                leaf = l
                prefix = c
        return leaf, prefix  


    def edges(self, node=None, content=None):
        if not content:
            content = []
        if not node:
            node = self.root
            
        for leaf in node.leaves:
            pattern = self.pat_map[leaf.sid]
            content.append(pattern[leaf.pos:])    
        for subnode in node.nodes:
            content.append(subnode.value[0])
            self.edges(subnode, content)
        return content
        


def suffix_tree_pieces(patterns):
    '''
    creating a suffix tree from all the patterns

    returns the suffix tree

    :param patterns: the strings with which to create the suffix tree
    :type patterns: [`str`, ...]

    :rtype: STree
    '''
    stree = STree(patterns)
    return stree.edges()
