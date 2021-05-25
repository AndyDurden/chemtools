#!/bin/python
import sys
import math


# Usage: python activechoice.py cis.out threshold homo maxroot_to_include1 maxroot_to_include2 ...  

thresh = float(sys.argv[2])
f = open(sys.argv[1], 'r')
currline = f.readline().split()
roots = []
currroot = 0
while (len(currline) < 2 or currline[0]+currline[1] != "FinalExcited"):
  if ((not len(currline)<5) and currline[0]+currline[2] == "RootLargest"):
    newroot = []
    currline = f.readline().split()
    while currline != []:
      if math.fabs(float(currline[8])) >= thresh:
        newroot.append(int(currline[0]))
        newroot.append(int(currline[2]))
      currline = f.readline().split()
    roots.append(newroot)
  currline = f.readline().split()  
  
print(roots)
homo = int(sys.argv[3])
for root in sys.argv[4:]:
  mino = 999999
  maxo = 0
  for i in roots[:int(root)]:
    mino = min(mino, min(i))
    maxo = max(maxo, max(i))
  total_space = (maxo-mino)+1
  occupied = (homo-mino)+1
  ndet = (math.factorial(total_space)/(math.factorial(occupied)*math.factorial(total_space-occupied)))**2
  print("For root: "+str(root))
  print("min: "+str(mino)+", max: "+str(maxo))
  print(str(total_space)+" orbital, "+str(occupied*2)+" electron active space")
  print("N Determinants: "+str(ndet)+"\n")








  
