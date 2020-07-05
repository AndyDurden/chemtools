#!/bin/python

import sys
import itertools
import numpy as np
f = np.math.factorial
def ch(n,k):
  return f(n)/(f(k)*f(n-k))

occ = [2,4,5,6,7]
unocc = [2,3,4,5]
homo = 21



jp = open("../job_params.ini", 'w')
xp = open("../xyz_params.ini", 'w')

xp.write("[pyrazine.xyz]\ncisnostates=4\ncasstates=4\nbasis=6-31gss\ncassinglets=4\n")

jobtypes = ["cisno_casci.in", "sacasscf.in", "ivo_casci.in"]

for jt in jobtypes:
  for e in itertools.product(occ,unocc):
    o, u = str(e[0]), str(e[1])
    c, a = str(homo-e[0]), str(e[0]+e[1])
    # Make sure our active space isnt too big
    if ( ch(int(a),int(o))**2 ) > 500000: continue
    xp.write("closed_"+o+"_"+u+"="+c+"\n")
    xp.write("active_"+o+"_"+u+"="+a+"\n")
    jp.write("["+jt.split(".")[0]+o+"_"+u+"]\n")
    if jt == "sacasscf.in": jp.write("casstates=casstates\n")
    if jt == "cisno_casci.in": jp.write("cisnostates=cisnostates\n")
    jp.write("template="+jt+"\n")
    jp.write("basis=basis\n")
    jp.write("cassinglets=cassinglets\n")
    jp.write("closed=closed_"+o+"_"+u+"\n")
    jp.write("active=active_"+o+"_"+u+"\n")


xp.close()
jp.close()
