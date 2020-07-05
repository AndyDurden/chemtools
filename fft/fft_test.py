#!/bin/python

import numpy as np
import subprocess
# time-series data should be in attoseconds

# frequency-series data should be in eV


c=2.99792458*10**10 # in cm/s
h=4.1357*10**-15 # in eV*s

# in fft we take attosecond in -> inv attosecond out
# then multiply by h*10**18 to get eV

# frequencies should be in eV i guess
def gen_complex_exponential_data(t, amplist, freqlist, phaselist):
  n = len(amplist)
  freqlist_ias = []
  for f in freqlist:
    freqlist_ias.append((f/h)*10**-18)
  signal = np.zeros(len(t), dtype=np.complex128)
  for i in range(0,n):
    signal += amplist[i] * np.exp(-1j*2*np.pi*freqlist_ias[i] *t + phaselist[i])
  return signal


def compare_fft_result(t,amplist,freqlist,phaselist, fft_file, tol=4):
  fil = open(fft_file,'r')
  print(fft_file+"\n")
  T = []
  F = []
  Fr = []
  Fi = []
  for line in fil:
    T.append(float(line.split(",")[0]))
    F.append(float(line.split(",")[3]))
    Fr.append(float(line.split(",")[1]))
    Fi.append(float(line.split(",")[2]))
  fil.close()
  delT = abs(T[2]-T[1])
  print("Frequency Resolution: "+str(delT)+" eV\n")
  n = len(amplist)
  for i in range(0,n):
    j = 0
    while (T[j] < freqlist[i]-tol*delT): j+=1
    start = j
    Fabssum = 0
    Frsum = 0
    Fisum = 0
    print("\nPeak with A: "+str(amplist[i])+", f: "+str(freqlist[i])+", phi: "+str(phaselist[i])+"\n")
    print("    f (eV), Fr,   Fi,   abs(F)\n")
    while (j< 1+start+2*tol):
      print("    {:0.2f}, {:0.2f} + {:0.2f}i,  {:0.2f}".format(T[j],Fr[j],Fi[j],F[j]))
      #print("    "+str(T[j])+", "+str(F[j])+"\n")
      Fabssum+=F[j]
      Frsum+=Fr[j]
      Fisum+=Fi[j]
      j+=1
    Fsumabs = np.sqrt( Frsum**2 + Fisum**2)
    print("\n  AbsSum: {:0.2f}, A error: {:0.2f} ({:0.2f}%)".format( Fabssum, (Fabssum-amplist[i]), 100.*abs(Fabssum-amplist[i])/amplist[i]) )
    print("  Sum: {:0.2f} + {:0.2f}i, Sum norm: {:0.2f}, A error: {:0.2f} ({:0.2f}%)\n".format( Frsum, Fisum, Fsumabs, (Fsumabs-amplist[i]), 100.*abs(Fsumabs-amplist[i])/amplist[i]) )
    #print("  Sum: "+str(Fsum)+", A error: "+ str( Fsum - amplist[i] )+" ("+str( abs(Fsum-amplist[i])/amplist[i])+"%)\n")
  return 0


def static_test():
  t = np.linspace(0,100000,100000)
  amplist = np.array([1.,2.,2.,2.])
  freqlist = np.array([1.,2.,3.,4.])
  #phaselist = np.array([0.,0.,np.pi/2.,0.])
  #phaselist = np.array([0.,np.pi,0,0.])
  phaselist = np.array([0.,0.,0.,0.])
  sig = gen_complex_exponential_data(t,amplist,freqlist,phaselist)
  f = open('test_sig','w')
  f.write("  time (as)   real   imaginary  \n")
  for i in range(0,len(t)):
    f.write("".join([str(t[i]),",",str(np.real(sig[i])),",",str(np.imag(sig[i])),"\n"]))
  f.close()
  p = subprocess.Popen("fft -z 200000 --xmin 0 --xmax 6 --ymin 0 --ymax 4 --outname statictest_raw test_sig",shell=True)
  p.wait()
  print("fft return code: "+str(p.poll()))
  compare_fft_result(t, amplist, freqlist, phaselist, "statictest_raw_fft.csv")
  print("\n--------\n")
  p = subprocess.Popen("fft -z 200000 -m --xmin 0 --xmax 6 --ymin 0 --ymax 4 --outname statictest_m test_sig",shell=True)
  p.wait()
  print("fft return code: "+str(p.poll()))
  compare_fft_result(t, amplist, freqlist, phaselist, "statictest_m_fft.csv")
  print("\n--------\n")
  p = subprocess.Popen("fft -z 200000 -w hamming --xmin 0 --xmax 6 --ymin 0 --ymax 4 --outname statictest_ham test_sig",shell=True)
  p.wait()
  print("fft return code: "+str(p.poll()))
  compare_fft_result(t, amplist, freqlist, phaselist, "statictest_ham_fft.csv")
  print("\n--------\n")
  p = subprocess.Popen("fft -z 200000 -m -w hamming --xmin 0 --xmax 6 --ymin 0 --ymax 4 --outname statictest_wham test_sig",shell=True)
  p.wait()
  print("fft return code: "+str(p.poll()))
  compare_fft_result(t, amplist, freqlist, phaselist, "statictest_wham_fft.csv")
  print("\n--------\n")

  return 0


static_test()
