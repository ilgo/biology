import heapq
from functools import total_ordering

@total_ordering
class PatCount(object):
    def __init__(self, pattern):
        self.pattern = pattern 
        self.count = 0
    def __eq__(self, other):
        if other != None and isinstance(other, PatCount):
            return self.pattern == other.pattern
        return False
    def __lt__(self, other):
        if other != None and isinstance(other, PatCount):
            return self.count < other.count
        return False
    def __str__(self):
        return 'p:{pattern};  c:{count}'.format(pattern=self.pattern, count=self.count) 

class KmerSet(object):
    def __init__(self):
        self.patCnts = []

    def __contains__(self, pc):
        for kmers in self.patCnts:
            if kmers.pattern == pc.pattern:
                return True
        return False

    def insert(self, pc):
        idx = -1
        try:
            idx = self.patCnts.index(pc)
        except(ValueError):
            pass

        if idx >= 0:
            pc = self.patCnts.pop(idx)
        pc.count += 1
        heapq.heappush(self.patCnts, pc)

    def head(self):
        return self.patCnts[0]

    def __str__(self):
        return '[' + ' '.join([str(pc) for pc in self.patCnts]) + ']'   


