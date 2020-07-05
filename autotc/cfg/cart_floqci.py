#!/bin/python

#e01_au = 0.11403342059492161 # ag2
e01_au = 0.06361450956239878 # ag4
d_au = 0.001 # equal to 0.068 eV
n=50



jp = open("../job_params.ini", 'w')
xp = open("../xyz_params.ini", 'w')

# s = "ag2"
s = "ag4"

xp.write("["+s+".xyz]\n")
xp.write(str(e01_au)+"="+str(e01_au)+"\n")
jp.write("["+"{0:.4}".format(e01_au)+"]\ntemplate="+s+"_floqci.in\nphotoneng="+str(e01_au)+"\n")

for i in range(1,n+1):
  photoneng_p = e01_au + i*d_au
  xp.write(str(photoneng_p)+"="+str(photoneng_p)+"\n")
  jp.write("["+"{0:.4}".format(photoneng_p)+"]\ntemplate="+s+"_floqci.in\nphotoneng="+str(photoneng_p)+"\n")
  photoneng_m = e01_au - i*d_au
  xp.write(str(photoneng_m)+"="+str(photoneng_m)+"\n")
  jp.write("["+"{0:.4}".format(photoneng_m)+"]\ntemplate="+s+"_floqci.in\nphotoneng="+str(photoneng_m)+"\n")

xp.close()
jp.close()
