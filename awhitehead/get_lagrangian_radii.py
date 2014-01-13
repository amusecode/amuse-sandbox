import sys
import re

time = []
lag_rad = [
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    []]


for fname in sys.argv[1:]:
    f = open(fname, "r")
    for line in f:
        m = re.search("%%% time= (\d+\.\d+)", line)
        if m:
            time.append(float(m.group(1)))
        m = re.search("Rlagr\[9\]= (\d+\.\d+) (\d+\.\d+) (\d+\.\d+) "+\
                      "(\d+\.\d+) (\d+\.\d+) (\d+\.\d+) (\d+\.\d+) "+\
                      "(\d+\.\d+) (\d+\.\d+)", line)
        if m:
            for i in range(0,9):
                lag_rad[i].append(float(m.group(i+1)))

for j in range(0, len(time)):
    s = "%.2f," % time[j]
    for i in range(0,9):
        s += "%10e," % lag_rad[i][j]
    print s

