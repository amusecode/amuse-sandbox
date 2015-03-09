import profile
import pstats
import sys

if __name__ == '__main__':
    stats = pstats.Stats(sys.argv[1])
    stats.strip_dirs()
    stats.sort_stats('cumulative')
    
    print 'INCOMING CALLERS:'
    stats.print_callers()
    
    print 'OUTGOING CALLEES:'
    stats.print_callees()
import profile

import pstats

import sys


