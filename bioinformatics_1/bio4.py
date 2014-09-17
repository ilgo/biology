from bio2 import cyclic_cutter
from collections import defaultdict
from graph import GraphFactory, euler_cycle, is_balanced

concat = lambda ps: ','.join([p for p in ps]) if isinstance(ps, list) else ps 


def str_composition(text, n):
    '''
    get  all substrings of length n

    retuen a list of substrings len(n)

    :param text: the string to analyse
    :type text: `str`
    :param n: lenght of substrings
    :type n: `int`

    :rtype: [`str`, ...]
    '''
    pos = 0
    while pos < len(text) - n + 1:
        yield text[pos:pos+n]
        pos += 1


def format_overlaps(adjacency_dict):
    '''
    format an overlap dictionary for Rosalind

    returns a formateted list of lexicographically sorted overlaps. eg:'ATCGG' -> 'TCGGT'

    :param adjacency_dict: a dict of overlaps
    :type adjacency_dict: {'source:'target', ...}
    
    :rtype: `str`
    '''
    out = ['%s -> %s' %(src, concat(adjacency_dict[src])) for src in sorted(adjacency_dict.keys())]
    return '\n'.join(out)


def str_overlap(kmers):
    '''
    create an adjacency graph for the given overlapping kmers

    return a dict of overlapping kmers. eg:'ATCGG' -> 'TCGGT'
    
    :param kmers: the overlapping kmers
    :type kmers: [`str`, ...]

    :rtype: {'source:'target', ...}
    '''
    overlaps = {}
    for kmer in kmers:
        suffix = kmer[1:]
        if suffix not in overlaps:
            overlaps[suffix] = [kmer, []]   
        else:
            overlaps[suffix][0] = kmer

        prefix = kmer[:-1]
        if prefix not in overlaps:
            overlaps[prefix] = ['', [kmer]]
        else:
            overlaps[prefix][1].append(kmer)
                
    return {k:v for k, v in overlaps.values() if k and v}


def deBruijn(kmers):
    '''
    constructs a deBruijn graph from a list of kmers

    returns an adjacency deBruijn graph 

    :param kmers: an iterable returning kmers
    :type k: `iter`

    :rtype: {'source: ['target'...], ...}
    '''
    overlaps = {}
    for kmer in kmers:
        suffix = kmer[1:]
        prefix = kmer[:-1]
        if prefix not in overlaps:
            overlaps[prefix] = [suffix]   
        else:
            overlaps[prefix].append(suffix)
    return {k:sorted(list(set(v))) for k, v in overlaps.items()}


def eulerCycle(edge_iter, start_node=0):
    '''
    find an euler cycle from a directed graph

    returns a list of nodes that form an euler cycle

    :param edge_iter: an iterable for edges
    :type edge_iter: `iter` -> (`int, [`int`, ...])
    :param start_node: for debugging purpose the node to start the tour
    :type start_node: `int`

    :rtype: [`int`, ...]
    '''
    graph = GraphFactory(edge_iter, directed=True)
    if not is_balanced(graph):
        return []
    return euler_cycle(graph, start_node)
    

def eulerPath(edge_iter):
    '''
    find an euler path in a directed graph

    returns a list of nodes that form an euler path

    :param edge_iter: an iterable for edges
    :type edge_iter: `iter` -> (`int, [`int`, ...])

    :rtype: [`int`, ...]
    '''
    graph = GraphFactory(edge_iter, directed=True)
    in_outs = {}
    for k,vs in graph.items():
        entry = in_outs[k] if k in in_outs else [0,0]
        entry[1] += len(vs)
        in_outs[k] = entry
        for v in vs:
            entry2 = in_outs[v] if v in in_outs else [0,0]
            entry2[0] += 1
            in_outs[v] = entry2
   
    start, end = None, None
    for k, inout in in_outs.items():
        diff = inout[0] - inout[1]
        if diff < 0:
            start = k
        elif diff > 0:
            end = k
        if start and end:
            break
    graph[end].append(start)
    path = euler_cycle(graph)
    start_idx = path.index(end) + 1
    return path[start_idx:-1] + path[:start_idx]
    

def str_reconstruct(kmers):
    '''

    '''
    edge_iter = ((k,v) for k, v in deBruijn(kmers).items())
    kmers = eulerPath(edge_iter)
    return kmers[0][:-1] + ''.join([kmer[-1] for kmer in kmers])
    
