#!/bin/python

e01_eV = 14.24336102
detune_eV = [-5.,-2.5, -1, -.5, -.25, -.1, -0.01, 0.01, 0.1, 0.25, 0.5, 1, 2.5, 5]
eng_au = []
h_eVs = 4.13566766225*10**(-15)
eV_to_au = 27.21138602

for delta in detune_eV:
  eng_au.append((e01_eV+delta)/eV_to_au  )


nsteps = [10000, 25000, 50000, 100000, 200000]


jp = open("../job_params.ini", 'w')
xp = open("../xyz_params.ini", 'w')

xp.write("[ethylene.xyz]\n")
for nstep in nsteps: xp.write(str(nstep)+"="+str(nstep)+"\n")
for i in range(0,len(eng_au)):
  xp.write( str(detune_eV[i]).replace("-","n")+"="+"{:.9E}".format(eng_au[i])+"\n")
  for nstep in nsteps:
    jp.write("["+str(detune_eV[i]).replace("-","n")+"_"+str(nstep)+"]\ntemplate=hf_floquet_tdcasci.in\n")
    jp.write("photoneng="+str(detune_eV[i]).replace("-","n")+"\nnstep="+str(nstep)+"\n")
    
xp.close()
jp.close()
