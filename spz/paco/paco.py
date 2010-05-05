import numpy as np

def auto_correlate(pattern, treshold) : 
  size = len(pattern)
  data = np.zeros(size)
  for ip in range(size) :
    for id in range(size) :
      if pattern[id] == pattern[(ip+id)%size] :
        data[ip] += 1
    data[ip] = data[ip]/size
  for id in range(len(data)):
    if data[id] >= treshold:
      data[id] = 1
    else :
      data[id] = 0
  return data

def count_pattern_length(ACdata) :
  count = 1
  while ACdata[count] == 0 and count<len(ACdata) :
    count += 1
  if count==1 or count == len(ACdata) :
    count = 0
  return count

# Cumbing fonction delta which identifies sections of the signal that
# contain the repeating pattern p of length lp.
def calculate_delta_function(ACdata, pattern, lp) :
  for id in range(len(pattern-2*lp)):
    if pattern[id]==pattern[lp+1]:
      ACdata[i] = 1
  return ACdata

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

t = []
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

import interface as paco

p = paco.PACOInterface()
print "Processing pattern..."
pattern = p.PatternConstruct(x, y, vx, vy)
#print pattern

treshold = 0.95
print "Autocorrelate..."
ACdata = auto_correlate(pattern, treshold)
pattern_length = count_pattern_length(ACdata)
unique_pattern = pattern[0:pattern_length]
#print unique_pattern

print "Check consistency..."
count_pattern = check_resulting_pattern(pattern, unique_pattern)
print "number of repeated patterns: ", count_pattern

print "PACO output (treshold=", treshold, "): "
print "   pattern: ",  unique_pattern 
print "   length:  ",  pattern_length
print "   repeats: ",  count_pattern

