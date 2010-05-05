import numpy as np
import interface as paco

def count_pattern_length(ACdata) :
  count = 1
  while ACdata[count] == 0 and count<len(ACdata) :
    count += 1
  if count==1 or count == len(ACdata) :
    count = 0
  return count

def check_resulting_pattern(pattern, unique_pattern) :
  lp = len(unique_pattern)
  il=0
  ir=lp
  count_pattern = 0
  for i in range(int(len(pattern)/lp)) :
    if ir < len(pattern) and pattern[il:ir] == unique_pattern :
      count_pattern += 1
    il = ir
    ir += lp
  return count_pattern

f = open("src/signal.dat")
lines = f.readlines()

t = [] # not used
x = []
y = []
vx = []
vy = []
for line in lines:
  l = line.strip().split()
  if len(l) >= 5:
    t.append(float(l[0]))
    x.append(float(l[1]))
    y.append(float(l[2]))
    vx.append(float(l[3]))
    vy.append(float(l[4]))

p = paco.PACOInterface()
print "Processing pattern..."
pattern = p.PatternConstruct(x, y, vx, vy)

treshold = 0.95
print "Autocorrelate..."
ACdata = p.auto_correlate(pattern, treshold)
pattern_length = count_pattern_length(ACdata)
unique_pattern = pattern[0:pattern_length]

print "Check consistency..."
count_pattern = check_resulting_pattern(pattern, unique_pattern)
print "number of repeated patterns: ", count_pattern

print "PACO output (treshold=", treshold, "): "
print "   pattern: ",  unique_pattern 
print "   length:  ",  pattern_length
print "   repeats: ",  count_pattern
