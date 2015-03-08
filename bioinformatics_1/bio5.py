from itertools import product
from bio1 import pp
from collections import deque
import sys, math, util
from array import array

prefixes = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']


def pairwise(t):
    it = iter(t)
    return zip(it,it)


def change(value, coins):
    '''
    Calculate the minimum number of coins needed to sum up to a given value.

    returns the minimum number of coins

    :param value: the value the coins need to be summed to.
    :type value: `int`
    :param coins: a list of possible coin values
    :type coins: [`int`, ...]

    :rtype: `int`
    '''
    count = 1
    changes = {coin:count for coin in coins} 
    data = coins[:]
    while value not in changes:
        sums = []
        count += 1
        for d, c in product(*[data, coins]):
            coin_sum = d + c
            if coin_sum not in changes:
                changes[coin_sum] = count
                sums.append(coin_sum)
        data = sums[:]
        
    return changes[value]   

def read_manhatten_data(file_path):
    '''
    read a manhatten data file (rosalind_5b)

    returns tuple of grid, down, right arrays

    :param file_path: the location of the data file
    :type file_path: `str`

    :rtype: ([[0, ...],[0, ...], ...], [[`int`, ...],[`int`, ...], ...], [[`int`, ...],[`int`, ...], ...])
    '''
    with open(file_path, 'r') as f:
        X, Y = (int(n) for n in f.readline().split(' '))
        grid = [[0] * (Y+1) for _ in range(X+1)]
        down = [] 

        for _ in range(X):
            down.append([int(n) for n in f.readline().split(' ')])

        # drop the '-\n' line
        f.readline()

        right = []
        for _ in range(X+1):
            right.append([int(n) for n in f.readline().split(' ')])
        return (grid, tuple(down), tuple(right))


def manhatten_path(grid, down, right):
    '''
    calculate the longest path from (0,0) -> (x,y)

    returns longest path

    :param grid: the grid initalized to all zeros
    :type grid: [[`int`, ...],[`int`, ...], ...]
    :param down: the path-length from up to below nodes
    :type down: [[`int`, ...],[`int`, ...], ...]
    :param right: the path lenght from left to right nodes
    :type right: [[`int`, ...],[`int`, ...], ...]

    :rtype: `int`
    '''
    x_len = len(grid)
    y_len = len(grid[0])

    for n in range(1, y_len):
        grid[0][n] = grid[0][n-1] + right[0][n-1]
    for n in range(1, x_len):
        grid[n][0] = grid[n-1][0] + down[n-1][0]
   
    for x in range(1, x_len):
        for y in range(1, y_len):
            vertical_sum   = grid[x-1][y] + down[x-1][y]
            horizontal_sum = grid[x][y-1] + right[x][y-1]
            grid[x][y] = max(vertical_sum, horizontal_sum)    

    return grid[x_len-1][y_len-1]



def LCS(dna1, dna2):
    '''
    calculates a Longest Common Subsequence (LCS) for 2 dna strings

    :param `str` dna1: a dna string
    :param `str` dna2: another dna string
    :return: a LCS string
    '''
    dag = _build(dna1, dna2)
    lcs = _track(dag, dna1)
    return lcs


def _build(s1, s2):
    '''
    build a DAG data structure, that also can hold the backtrack pointers
    the nodes in the returned list are already topologically sorted.

    :param `str` s1: a string
    :param `str` s2: another string
    :returns: [{'val':0, 'src':0}, ...] 
    :return: the DAG graph from the 2 strings 
    '''
    nodes = (len(s1)+1) * (len(s2)+1)
    dag = {n:{'val':0, 'src':-1} for n in range(nodes)}
    #print(len(dag), len(s1), len(s2))    

    # fill top row, has only horizontal edge 
    for i in range(1, len(s1)+1):
        node = dag[i]
        node['src'] = i-1

    # fill all other rows
    row_len = len(s1) + 1
    for i in range(len(s1)+1, nodes):
        div, rem = divmod(i, row_len)
        node = dag[i]
        if rem == 0:
            v_char = s2[div-1]
            node['src'] = i - row_len
        else:
            diag = 1 if s1[rem-1] == v_char else 0 
            edge1 = (dag[i-1]['val'], i-1)
            edge2 = (dag[i-row_len]['val'], i-row_len)
            edge3 = (dag[i-row_len-1]['val'] + diag, i-row_len-1)
            edge = max(edge1, edge2, edge3) 
            node['val'] = edge[0]
            node['src'] = edge[1]

    return dag


def _track(dag, s1):
    '''
    find a longest common subsequence

    :param dag: the DAG Nodes
    :type dag: [:class:`DagNode`, ...]
    :rtype: `str`
    :returns: a longest common subsequence LCS
    '''
    lcs = []
    row_len = len(s1) + 1
    idx = max(dag.keys())
    while idx > 0:
        node = dag[idx]
        #print(idx, node, lcs)
        #print('------------')
        if idx - node['src'] == row_len + 1:
            j = idx % row_len
            lcs.append(s1[j-1])
        idx = node['src']

    return ''.join(lcs[::-1])


# DAG longest path
#--------------------

def make_edge(edge_str):
    path, weight = edge_str.strip().split(':')
    src, tar = path.split('->')
    return {'f':int(src), 't':int(tar), 'w':int(weight)}    


def DAG(sink, edges, levels=1):
    '''
    create a DAG representation from the edges
    start and end nodes are also known

    :param `int` levels: equals 3 for affine_gap triple plane structure
    
    assumes that sink is the highest numbered node
    '''
    dag = {i:array('i') for i in range((sink+1)*levels)}
    for edge in edges:
        dag[edge['f']].append(edge['t'])
        dag[edge['f']].append(edge['w'])
        #print("f:%d\tt:%d\tw:%d" %(edge['f'], edge['t'], edge['w']))
    return dag

def longest_path(src, sink, dag, mins=-sys.maxsize):
    done = {}
    lps = {src:(0, None)} # val, from
    q = deque([src])
    while len(q) > 0:
        node = q.popleft()
        if node in dag:
            edges = dag[node]
            for edge_0, edge_1 in pairwise(edges):
                #if edge_0 not in done:
                #q.append(edge_0)
                #done[edge_0] = None
                val = lps.get(node, (mins, None))
                val_to  = lps.get(edge_0, (mins, None))
                if val[0] + edge_1 > val_to[0]:
                    lps[edge_0] = (val[0] + edge_1, node)  
                    q.append(edge_0)
    #print ('\n'.join(['%r : %r ' %(k,v) for k,v in lps.items()]))
    del dag
    back = [sink]
    edge = lps[sink]
    while edge[1]:
        back.append(edge[1])
        edge = lps[edge[1]]
    back.append(0)
    return lps[sink][0], back[::-1]


def longest_path_fmt(res):
    out = [str(res[0])]
    out.append('->'.join([str(p) for p in res[1]]))
    return '\n'.join(out)
   


#-----------------------------------------------------------
# global align with score matrix

def read_score_matrix(matrix_name):
    with open('data/' + matrix_name + '.txt', 'r') as f:
        matrix = []
        for line in f.readlines():
            entry = [int(i) for i in line.split(',')]
            matrix.append(entry)
        return matrix


def _make_edges(s1, s2, matrix, sigma, source_connect=False, fit_connect=False, overlap_connect=False):
    '''
    :param `str` s1: the first string
    :param `str` s2: the second string
    :param matrix: the scoring matrix, scores and penalties for matches and mismatches
    :type matrix: [[`int,...][`int`,..], ...]
    :param `int` sigma: penalty for indels
    :param `bool` source_connect: if True, connect every node with source and sink.
    :param `bool` fit_connect: connect first row nodes to source and last row nodes to sink. 
    :param `bool` overlap_connect: connect first row nodes to source and last column nodes to sink. 
    '''
    print(s1, s2)
    n = (len(s1)+1) * len(s2)
    for i in range(n):
        tar = i + len(s1)+1
        x,y = divmod(i, len(s1)+1)
        if source_connect and i > 0 and i < n:
            yield {'f':0, 't':i, 'w':0}  
            yield {'f':i, 't':n+len(s1), 'w':0} 

        elif (fit_connect or overlap_connect) and i > 0 and i < len(s1)+1:
            yield {'f':0, 't':i, 'w':0}  
 
        if overlap_connect and i % (len(s1)+1) == len(s1):
            yield {'f':i, 't':n+len(s1), 'w':0} 

        if y == len(s1):
            # last column
            yield {'f':i, 't':tar, 'w':sigma}  
        else:
            yield {'f':i, 't':i + 1, 'w':sigma}
            #if not fit_connect:  
            yield {'f':i, 't':tar,   'w':sigma}
            c1 = prefixes.index(s1[y])
            c2 = prefixes.index(s2[x])
            yield {'f':i, 't':tar + 1, 'w':matrix[c1][c2]} 
    # last row
    for i in range(n, n+len(s1)):
        yield {'f':i, 't':i + 1, 'w':sigma}
        if source_connect or fit_connect:
            yield {'f':i, 't':n+len(s1), 'w':0}  


def align_strings(s1, s2, path):
    # s1 is vertical
    # s2 is horizontal
    as1, as2 = [], []
    rev_path = path[::-1]
    for i, p1 in enumerate(rev_path[:-1]):
        p2 = rev_path[i+1]
        x1, y1 = divmod(p1, len(s1)+1)
        x2, y2 = divmod(p2, len(s1)+1)
        #print(x1,y1, x2, y2)
    
        if abs(x1 - x2) == 1 and abs(y1 - y2) == 1:
            as1.append(s1[y2])
            as2.append(s2[x2])

        elif x1 == x2:
            as2.append('-')
            as1.append(s1[y2])

        elif y1 == y2:
            as2.append(s2[x2])
            as1.append('-')

    return [''.join([c for c in s[::-1]]) for s in (as1, as2)]
  


def global_align(s1, s2, sigma=-5):
    blosum = read_score_matrix('blosum62')
    edges = _make_edges(s1, s2, blosum, sigma)
    #print ('\n'.join([str(e) for e in edges]))
    sink = (len(s1) + 1) * (len(s2) + 1) - 1 
    dag = DAG(sink, edges)
    #print('\n'.join(['%d: %r' %(k,v) for k,v in dag.items()]))
    length, path = longest_path(0, sink, dag)
    return (length, path)


def local_align(s1, s2, sigma=-5, matrix='pam250'):
    pam250 = read_score_matrix(matrix)
    edges = _make_edges(s1, s2, pam250, sigma, source_connect=True)
    #print ('\n'.join([str(e) for e in edges]))
    sink = (len(s1) + 1) * (len(s2) + 1) - 1 
    dag = DAG(sink, edges)
    #print('\n'.join(['%d: %r' %(k,v) for k,v in dag.items()]))
    length, path = longest_path(0, sink, dag)
    if path[1] < len(s1)+1:
        path = path[1:]
    return (length, path)


def levenstein(s1,s2):
    matrix = [[-1 for i in range(len(prefixes))] for n in range(len(prefixes))]
    for n in range(len(prefixes)):
        matrix[n][n] = 0
    edges = _make_edges(s1, s2, matrix, -1)
    sink = (len(s1) + 1) * (len(s2) + 1) - 1 
    dag = DAG(sink, edges)
    length, path = longest_path(0, sink, dag)
    return -length
        

def fit_align(s1,s2):
    matrix = [[-1 for i in range(len(prefixes))] for n in range(len(prefixes))]
    for n in range(len(prefixes)):
        matrix[n][n] = 1
    edges = _make_edges(s1, s2, matrix, -1, fit_connect=True)
    sink = (len(s1) + 1) * (len(s2) + 1) - 1 
    dag = DAG(sink, edges)
    length, path = longest_path(0, sink, dag)
    if path[1] < len(s1)+1:
        path = path[1:]
    if path[-1] - path[-2] < len(s1)+1:
        path = path[:-1] 
    return (length, path)

def overlap_align(s1, s2):
    matrix = [[-2 for i in range(len(prefixes))] for n in range(len(prefixes))]
    for n in range(len(prefixes)):
        matrix[n][n] = 1
    edges = _make_edges(s1, s2, matrix, -2, overlap_connect=True)
    sink = (len(s1) + 1) * (len(s2) + 1) - 1 
    dag = DAG(sink, edges)
    length, path = longest_path(0, sink, dag)
    if path[1] < len(s1)+1:
        path = path[1:]
    #if path[-1] - path[-2] < len(s1)+1:
    path = path[:-1] 
    return (length, path)


def affine_gap(s1, s2, gap_penalty=-11 , extend_penalty=-1 , matrix_name='blosum62', matrix=None):

    if not matrix:
        matrix = read_score_matrix(matrix_name)
    edges = _make_affine_edges(s1, s2, matrix, gap_penalty, extend_penalty)
    #print ('\n'.join([str(e) for e in edges]))
    sink = (len(s1) + 1) * (len(s2) + 1) - 1 
    dag = DAG(sink, edges, levels=3)
    #print('\n'.join(['%d: %r' %(k,v) for k,v in dag.items()]))
    length, path = longest_path(0, sink, dag)
    #if path[1] < len(s1)+1:
    #    path = path[1:]
    return (length, path)


def _make_affine_edges(s1, s2, matrix, gap_penalty, extend_penalty):
    '''
    Using mod #node_count to create the upper and lower planes 
    central plane -> diagonals
    upper plane   -> columns
    lower plane   -> rows

    :param `str` s1: the first string
    :param `str` s2: the second string
    :param matrix: the scoring matrix, scores and penalties for matches and mismatches
    :type matrix: [[`int,...][`int`,..], ...]
    :param `int` gap_penalty: penalty for starting an opening
    :param `int` extend_penalty: penalty for extending a gap
    '''
    n = (len(s1)+1) * len(s2)
    upper = n + len(s1) + 1
    lower = 2 * upper
    for i in range(n):
        x,y = divmod(i, len(s1)+1)
        tar = i + len(s1)+1
        if y == len(s1):
            # last column
            yield {'f':i, 't':tar+upper, 'w':gap_penalty}

            yield {'f':i+upper, 't':tar+upper, 'w':extend_penalty}

            yield {'f':i+upper, 't':i, 'w':0}
            yield {'f':i+lower, 't':i, 'w':0} 

        else:
            c1 = prefixes.index(s1[y])
            c2 = prefixes.index(s2[x])
            yield {'f':i, 't':tar + 1, 'w':matrix[c1][c2]}
 
            yield {'f':i, 't':tar+upper, 'w':gap_penalty}
            yield {'f':i, 't':i+lower+1, 'w':gap_penalty}

            yield {'f':i+upper, 't':tar+upper, 'w':extend_penalty}
            yield {'f':i+lower, 't':i+lower+1, 'w':extend_penalty} 
            
            if i > 0:
                yield {'f':i+upper, 't':i, 'w':0}
                yield {'f':i+lower, 't':i, 'w':0} 
             
    # last row
    for i in range(n, n+len(s1)):
        yield {'f':i, 't':i+lower+1, 'w':gap_penalty}

        yield {'f':i+lower, 't':i+lower+1, 'w':extend_penalty} 
        yield {'f':i+lower, 't':i, 'w':0} 
                
        yield {'f':i+upper, 't':i, 'w':0}
      
    i = n + len(s1)  
    yield {'f':i+upper, 't':i, 'w':0}
    yield {'f':i+lower, 't':i, 'w':0} 


def align_affine_strings(s1, s2, path):
    # s1 is vertical
    # s2 is horizontal
    
    upper = (len(s1) + 1) * (len(s2) + 1)
    lower = 2 * upper

    as1, as2 = [], []
    rev_path = path[::-1]
    for i, p1 in enumerate(rev_path[:-1]):
        p2 = rev_path[i+1]

        if p1 >= upper and p2 < upper:
            # returning to central plane is a NOOP
            #print('noop')
            continue

        x1, y1 = divmod(p1 % upper, len(s1)+1)
        x2, y2 = divmod(p2 % upper, len(s1)+1)
    
        if p1 < upper and p2 < upper:
            # on central plane
            #print('central')
            as1.append(s1[y2])
            as2.append(s2[x2])

        elif (p1 < upper or p1 >= lower) and p2 >= lower:
            # on lower plane
            #print('lower')
            as2.append('-')
            as1.append(s1[y2 % len(s1) - 1])

        elif p1 < lower and p2 >= upper and p2 < lower:
            # on upper plane
            #print('upper')
            as2.append(s2[x2 % len(s2) - 1])
            as1.append('-')

        #print(as1, as2)

    return [''.join([c for c in s[::-1]]) for s in (as1, as2)]
  

def middle_edge(s1, s2, sigma=-5, matrix='blosum62'):
    
    half = math.floor(len(s2)) / 2)
    mat = read_score_matrix(matrix)
    sink = (len(s1) + 1) * (half + 1) - 1 

    edges = _make_edges(s1, s2[:half], mat, sigma)
    dag = DAG(sink, edges)
    middle_nodes = [(k,v) for k,v in dag.items() if (k+1) % (half+1) == 0]
    print(util.pp(middle_nodes, join_char='\n'))
 
    #print('\n'.join(['%d: %r' %(k,v) for k,v in dag.items()]))
    lg_dict = longest_paths(0, sink, dag)
    #m_nodes = [k for k, v in middle_nodes]
    nodes = [(k,lg_dict[k]) for k,v in middle_nodes]
    print(nodes)
    print('---')

    r_s1, r_s2 = s1[::-1], s2[::-1]
    edges = _make_edges(r_s1, r_s2[:half], mat, sigma)
    dag = DAG(sink, edges)
    r_middle_nodes = [(k,v) for k,v in dag.items() if (k+1) % (half+1) == 0]
    r_lg_dict = longest_paths(0, sink, dag)
    rm_nodes = [k for k, v in r_middle_nodes]
    sm = rm_nodes[0] + rm_nodes[-1]
    r_nodes = [(sm-k,lg_dict[k]) for k,v in middle_nodes[::-1]]
    #r_nodes = [sm-n for n in r_nodes]

    print(util.pp(middle_nodes, join_char='\n'))
    print(r_nodes)


    return #(length, path)

def longest_paths(src, sink, dag, mins=-sys.maxsize):
    done = {}
    lps = {src:(0, None)} # val, from
    q = deque([src])
    while len(q) > 0:
        node = q.popleft()
        if node in dag:
            edges = dag[node]
            for edge_0, edge_1 in pairwise(edges):
                #if edge_0 not in done:
                #q.append(edge_0)
                #done[edge_0] = None
                val = lps.get(node, (mins, None))
                val_to  = lps.get(edge_0, (mins, None))
                if val[0] + edge_1 > val_to[0]:
                    lps[edge_0] = (val[0] + edge_1, node)  
                    q.append(edge_0)
    #print ('\n'.join(['%r : %r ' %(k,v) for k,v in lps.items()]))
    del dag
    back = [sink]
    edge = lps[sink]
    while edge[1]:
        back.append(edge[1])
        edge = lps[edge[1]]
    back.append(0)
    return lps
