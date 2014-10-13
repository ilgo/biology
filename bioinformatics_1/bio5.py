from itertools import product
from bio1 import pp

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


def LCS(data):
    '''
    calculates longest common subsequence (LCS) for 2 strings

    returns the LCS

    :param data:
    :type data:

    :rtype: `str`
    '''
    h_len = len(data[0])
    v_len = len(data[1])

    horizontal = [c for c in data[0]]
    vertical = [c for c in data[1]]
    grid = []
    for _ in range(v_len+1):
        grid.append([0] * (h_len+1))
    
    for v in range(v_len):
        for h in range(h_len):
            match = 1 if horizontal[h] == vertical[v] else 0
            #print(v,h,'match:', match)
            #print(v,h+1,'left')
            #print(v+1,h,'above')
            grid[v+1][h+1] = max([grid[v][h] + match, grid[v][h+1], grid[v+1][h]])
            #print(v+1,h+1,'point:',grid[v+1][h+1])
            #print('...')
        #print('>>>>>>>>>>>>>')

    #print(pp(grid, join_char='\n'))
    #return grid[v_len-1][h_len-1]

    lcs = []
    n = grid[v_len-1][h_len-1]
    v, h = v_len-1, h_len-1
    while n > 0:
        dia = grid[v-1][h-1]
        up  = grid[v-1][h]
        left = grid[v][h-1]
        if dia >= up and dia >= left:
            if dia < n:
                lcs.append(vertical[v-1])
            n = dia
            v = v-1
            h = h-1
            #print('dia')
        elif up >= left:
            v = v-1
            #print('up')
        else:
            h = h-1
            #print('left')

    return ''.join(lcs[::-1]) 
