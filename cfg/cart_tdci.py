#!/bin/python

e01_eV = 14.24336102
detune_eV = [-5.,-2.5, -1, -.5, -.25, -.1, -0.01, 0.01, 0.1, 0.25, 0.5, 1, 2.5, 5]
freq_hz = []
h_eVs = 4.13566766225*10**(-15)

for delta in detune_eV:
  freq_hz.append((e01_eV+delta)/h_eVs  )


nsteps = [10000, 25000, 50000, 100000, 200000]


jp = open("../job_params.ini", 'w')
xp = open("../xyz_params.ini", 'w')

xp.write("[ethylene.xyz]\n")
for nstep in nsteps: xp.write(str(nstep)+"="+str(nstep)+"\n")
for i in range(0,len(freq_hz)):
  xp.write( str(detune_eV[i]).replace("-","n")+"="+"{:.9E}".format(freq_hz[i])+"\n")
  for nstep in nsteps:
    jp.write("["+str(detune_eV[i]).replace("-","n")+"_"+str(nstep)+"]\ntemplate=hf_tdcasci.in\n")
    jp.write("laser_freq="+str(detune_eV[i]).replace("-","n")+"\nnstep="+str(nstep)+"\n")
    
xp.close()
jp.close()
