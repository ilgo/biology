#!/usr/bin/env python
from algo.find import seq, find_all_re
from datetime import datetime, timedelta


if __name__ == '__main__':

    begin = datetime.now()
    g = seq(10000000)

    for _ in range(10):
        s = seq(12)
        #print(g)
        print(s)
        start = datetime.now()
        matches = find_all_re(s, g)
        secs = (datetime.now() - start).total_seconds()
        print(secs, matches)
   
    print('Run for %0.5f sec' %(datetime.now() - begin).total_seconds())

