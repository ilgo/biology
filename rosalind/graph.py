from collections import defaultdict
from random import randint


def GraphFactory(edge_iter, directed=False):
    '''
    creates an adjacency datastructure from an edge-iterator

    returns a dict, representing the adjacency structure

    :param edge_iter: an iterable that returns edges as from,[to,..]
    :type edge_iter: `iter` -> (`int`, [`int`, ...])
    :param directed: is the input a directed graph
    :type directed: `bool`

    :rtype: `dict`
    '''
    data = defaultdict(list)
    for source, targets in edge_iter:
        data[source].extend(targets)
        if not directed:
            for target in targets:
                data[target].append(source)                                       
    return data
  

def g_read(lines):
    for line in lines:
        src, pointer, targets = line.strip().split(' ')
        targets = [int(target) for target in targets.split(',')]
        yield (int(src), targets) 


def connected_components(graph):
    visited = []
    components = 0
    for k, v in graph.items():
        if k not in visited:
            visited.append(k)
            todo = [vv for vv in v if vv not in visited]
            while todo:
                for t in todo:
                    visited.append(t)
                    todo.remove(t)
                    todo.extend([vv for vv in graph[t] if vv not in visited])
            components += 1
    return components


def is_balanced(graph):
    '''
    Checks if the graph is balanced with respect to in/out bound edges for each node.
    Assumes that node enumeration starts at 0.

    returns True if the graph is balanced 

    :param graph: a graph data structure
    :type graph: {`int` : [`int`, ...], ...}

    :rtype: `bool`
    '''
    table = [0] * len(graph)
    for k, vs in graph.items():
        table[k] += len(vs)
        for v in vs:
            table[v] += 1
    for n in table:
        if n % 2 == 1:
            return False    
    return True


def euler_cycle(graph, start_node=None):
    '''
    find an euler cycle in a directed graph

    returns a list of nodes that form an euler cycle, first and last node id are the same

    :param graph: a directed adjacency graph
    :type graph: {`int` : [`int`, ...], ...}
    :param start_node: for debugging purpose the node to start the tour
    :type start_node: `int`

    :rtype: [`int`, ...]
    '''
    path = []
    euler = {k:v for k,v in graph.items()}
    max_path_len = sum([len(v) for v in graph.values()]) + 1

    node = start_node if start_node else list(graph.keys())[0]
    path.append(node)
    while len(path) < max_path_len:
        try:
            next_node = euler[node].pop()
            node = next_node
            path.append(node)
        except:
            combined = list(zip(range(len(path)), path))[:-1]
            for i,last_node in combined[::-1]:
                if euler[last_node]:
                    path = path[i:] + path[1:i+1]
                    node = last_node
                    break 
    return path 
