#!/bin/python

import sys

import numpy as np
f = np.math.factorial
def ch(n,k):
  #print(n,k)
  return f(n)/(f(k)*f(n-k))


# returns the average energy of the n lowest states
def get_avge(filename, n):
	engs = []
	f = open(filename, 'r')
	l = f.readline()
	while l:
		if (len(l.split()) > 4):
			if l.split()[0]+l.split()[1]+l.split()[3]=="Singletstateenergy:":
				engs.append(float(l.split()[4]))
				l=f.readline()
			else: l=f.readline()
		else: l=f.readline()
	return( sum(engs[0:n])/float(n) )


jobtypes = ["cisno_casci", "sacasscf", "ivo_casci"]
xyzs = ["pyrazine"]
occ = ['2','4','5','6','7']
unocc = ['2','3','4','5']
jobdir = "/home/adurden/883/final/1/"
avgn = 4

def compile_data(jobtypes, xyzs, occ, unocc, jobdir, avgn):
  data = {}
  for jobtype in jobtypes:
    data[jobtype] = {}
    for xyz in xyzs:
      avge = []
      for i in range(0,len(occ)): avge.append([])
      for i in range(0,len(occ)):
        for j in range(0,len(unocc)):
          if ( ch( int(occ[i])+int(unocc[j]) , int(occ[i])  )**2 ) > 500000: continue
          # my clearly well thought out naming scheme
          jobpath = jobdir+jobtype+occ[i]+"_"+unocc[j]+"/"+xyz+"/"
          outname = xyz+"_"+jobtype+occ[i]+"_"+unocc[j]+".out"
          avge[i].append(get_avge(jobpath+outname, avgn))
      print(avge)
      data[jobtype][xyz] = avge
  return data

def write_csv(data, occ, unocc):
  for jobtype in data:
    for xyz in data[jobtype]:
      csvs = ',' # first cell blank
      for i in unocc: csvs=csvs+i+',' # headers
      csvs = csvs+"\n"
      for i in range(0,len(occ)):
        csvs = csvs+occ[i]+','
        for j in range(0,len(unocc)):
          if ( ch( int(occ[i])+int(unocc[j]) , int(occ[i])  )**2 ) > 500000: continue
          print(jobtype,xyz,i,j,data)
          csvs = csvs+str(data[jobtype][xyz][i][j])+','
        csvs = csvs+'\n'
    f = open(xyz+'_'+jobtype+'.csv','w')
    f.write(csvs)
    f.close()


data = compile_data(jobtypes, xyzs, occ, unocc, jobdir, avgn)
write_csv(data, occ, unocc)










