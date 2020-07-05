#!/bin/python
pnosteps = ['1','2','3']
dft_types = ['svwn', 'svwn1', 'svwn3', 'pz81', 'p86', 'pw91', 'pw91c', 'blyp', 'lypc', 
             'pbe', 'pbec', 'revpbe', 'wpbe', 'bop', 'mubop', 'pbe0', 'revpbe0',
             'b3lyp', 'b3lyp3', 'b3lyp5', 'b3p86', 'b3pw91', 'bhandhlyp', 'camb3lyp',
             'rcam-b3lyp', 'wpbeh', 'b97', 'wb97', 'wb97x']

jp = open("job_params.ini", 'w')
xp = open("xyz_params.ini", 'w')

xp.write("[ag4.xyz]\n1=1\n2=2\n3=3\n")
for dft in dft_types:
  xp.write(dft+"="+dft+"\n")
  for s in pnosteps:
    jp.write("["+dft+"_"+s+"]\ntemplate=pno_casci.in\n")
    jp.write("method="+dft+"\npno_iters="+s+"\n")
    
xp.close()
jp.close()
