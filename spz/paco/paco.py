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

def count_zeros(data) :
  count = 0
  for id in data :
    if id == 0 :
      count += 1
      if count==1 or count == len(data) :
        lp = 0
      else :
        lp = count
  return count

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
pattern = p.PatternConstruct(x, y, vx, vy)
print pattern

treshold = 0.95
data = auto_correlate(pattern, treshold)
lp = count_zeros(data)




