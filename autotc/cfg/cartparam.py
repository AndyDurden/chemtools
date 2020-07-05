#!/bin/python

import sys
import itertools
import numpy as np

Fmin = 0.005
Fmax = 0.300
Fn = 30

Imin=1
Imax=8

jp = open("../job_params.ini", 'w')
xp = open("../xyz_params.ini", 'w')

xp.write("[ag4.xyz]\nclosed=36\nactive=9\n")

for e in itertools.product(list(np.linspace(Fmin,Fmax,Fn)),range(Imin,Imax+1)):
  f = "{0:.3f}".format(e[0])
  i = str(e[1])
  xp.write("pnof_"+f+"_"+i+"="+f+"\n")
  xp.write("pnoi_"+f+"_"+i+"="+i+"\n")
  jp.write("[pno_"+f+"_"+i+"]\n")
  jp.write("template=pno_casci.in\n")
  jp.write("pno_polar_z=pnof_"+f+"_"+i+"\n")
  jp.write("pno_iters=pnoi_"+f+"_"+i+"\n")
  jp.write("closed=closed\n")
  jp.write("active=active\n")


xp.close()
jp.close()
