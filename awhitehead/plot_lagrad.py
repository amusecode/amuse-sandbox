import sys
import matplotlib as mpl
mpl.use('PDF')
from matplotlib.pyplot import *

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

outfname = sys.argv[1].split(".")[0] + ".pdf"
with open(sys.argv[1], "r") as f:
    for line in f:
        cols = line.split(",")
        if len(cols) > 1:
            time.append(float(cols[0]))
            for i in range(0,9):
                lag_rad[i].append(float(cols[i+1]))

for i in range(0,9):
    semilogy(time, lag_rad[i])
savefig(outfname, format="pdf")

