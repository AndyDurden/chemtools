#!/bin/python

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import numpy as np
import math
import sys
from optparse import OptionParser
usage = "usage: %prog [options] file1 file2 ..."
parser=OptionParser(usage=usage)
# General
parser.add_option("--outname", action="store", default=False, dest="outname", help="Set a name for associated output files like [outname]_fft.csv and [outname].png. Defaults to the input filename.")

# Pre-processing
parser.add_option("--nskip", action="store", default=1, dest="nskip", help="only include every Nth datapoint in input")
parser.add_option("--truncate", action="store", default=999999999, dest="truncate", help="dont include data after Nth datapoint")

parser.add_option("-c", action="store", default="none", dest="center_mode", help="Centering mode. Options: s - shift the first datapoint to zero. u - add a (0,1) point. none - Don't center the data.")
parser.add_option("-m", "--mirror", action="store_true", default=False, dest="mirror", help="mirror data before FFT")
parser.add_option("-w", "--windowfunc", action="store", default=False, dest="windowfunc", help="Select windowing function ") # hanning, hamming, bh, exp
parser.add_option("--rawproc", action="store_true", default=False, dest="rawproc", help="Write a csv with pre-processed input data before FFT")
parser.add_option("-z", action="store", default=0, dest="zeropad", help="Pad this many zeros to the end of the data (applied after window)")
#parser.add_option("-p", "--pulse", action="store_true", default=False, dest="pulse", help="For TDCI TCorrFn with pulses, data starts after the 'pulse ended' line.")

# Alternate FT algos
parser.add_option("--fftw", action="store_true", default=False, dest="fftw", help="Use pyfftw instead of numpy.fft.fft")
parser.add_option("--pade", action="store_true", default=False, dest="use_pade", help="Perform FT using method of Pade approximants.")
parser.add_option("--psig", action="store", default=False, dest="psig", help="'signal damp factor [windowing?], yield peaks with [sic: at maximum(?)] FWHM of 2/sigma' (???) Default: no signal dampening/windowing. Original code had default of 100.")
parser.add_option("--pmin", action="store", default=0.0,   dest="pmin", help="Pade: specify frequency calculation lower bound (in eV)")
parser.add_option("--pmax", action="store", default=15.0,  dest="pmax", help="Pade: specify frequency calculation upper bound (in eV)")
parser.add_option("--pstep", action="store", default=False, dest="pstep", help="Pade: specify number of points to generate in spectrum. Default: same as input data length")
# Plotting 
#parser.add_option("-l", "--limits", action="store", default=False, dest="limstr", help="Specify graph window in the format: [xmin],[xmax],[ymin],[ymax]   do not include spaces. (broken right now)")
parser.add_option("-t", "--title", action="store", default="", dest="title", help="Specify graph title")
parser.add_option("--no-plot", action="store_true", default=False, dest="no_plot", help="Skip generating a plot/graph.")
parser.add_option("--no-fft", action="store_true", default=False, dest="no_fft", help="Only plot data. [file1, file2, ...] are accepted as the output files (output_fft.csv) from a previous run")
parser.add_option("--xmin", action="store", default=False, dest="xmin", help="Specify minimum x of graph window")
parser.add_option("--xmax", action="store", default=False, dest="xmax", help="Specify maximum x of graph window")
parser.add_option("--ymin", action="store", default=False, dest="ymin", help="Specify minimum y of graph window")
parser.add_option("--ymax", action="store", default=False, dest="ymax", help="Specify maximum y of graph window")
parser.add_option("-e", "--excitations", action="store", default=False, dest="tcoutfile", help="Specify a TeraChem outfile to plot vertical excitations (with significant transition dipole) from.")
parser.add_option("--ag3ndata", action="store_true", default=False, dest="ag3ndata", help="Draw vertical lines at experimental Ag3- energies.")

# auto windowing stuff that never worked
#parser.add_option("-w", "--window-thresh", action="store", type="float", dest="window_thresh", help="truncate tail of data in plot where all Y below this value. Default 0.025", default=0.025 )
#parser.add_option("-W", "--window-manual", action="store", type="int", dest="window_manual", help="Crop window to this many points on X axis", default=0)
(options, args) = parser.parse_args()


#print(options.mirror)
#print(options.window_thresh)

c=2.99792458*10**10 # in cm/s
h=4.1357*10**-15 # in eV*s

td_thresh = 1.0
pmin, pmax = 0,0

truncate = int(options.truncate)

if options.nskip: options.nskip = int(options.nskip)
if options.pstep: options.pstep = int(options.pstep)
if options.pmin:
  pmin = int(options.pmin)
  pmin = (pmin/h)*10**-18
if options.pmax:
  pmax = int(options.pmax)
  pmax = (pmax/h)*10**-18
  #print("au range: "+str((pmin, pmax)))


#X_s=np.array(freqs)*10**18 # from inv attosecond to inv seconds
#X_eV = h*X_s # energy in eV




# Alternate FT algo from Lapota paper using Pade approximants
# Code stolen from Joshua Goings
# https://github.com/jjgoings/pade
# Related links:
# https://joshuagoings.com/assets/real-time-electronic-structure.pdf
# https://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00511
# https://en.wikipedia.org/wiki/Pade_approximant
def pade(time,signal,sigma,max_len=None,w_min=0.0,w_max=10.0,w_step=0.01):
  """ Routine to take the Fourier transform of a time signal using the method
    of Pade approximants.

    Inputs:
  time:      (list or Numpy NDArray) signal sampling times
  signal:    (list or Numpy NDArray) 

  Optional Inputs:
  sigma:     (float) signal damp factor, yields peaks with 
  FWHM of 2/sigma 
  max_len:   (int) maximum number of points to use in Fourier transform
  w_min:     (float) lower returned frequency bound
  w_max:     (float) upper returned frequency bound
  w_step:    (float) returned frequency bin width

  Returns:
  fsignal:   (complex NDArray) transformed signal
  frequency: (NDArray) transformed signal frequencies

  From: Bruner, Adam, Daniel LaMaster, and Kenneth Lopata. "Accelerated 
  broadband spectra using transition signal decomposition and Pade 
  approximants." Journal of chemical theory and computation 12.8 
  (2016): 3741-3750. 
  """

  # center signal about zero
  #signal = np.asarray(signal) - signal[0]

  stepsize = time[1] - time[0]

  # Damp the signal with an exponential decay. 
  if sigma:
    damp = np.exp(-(stepsize*np.arange(len(signal)))/float(sigma))
    signal *= damp

  M = len(signal)
  N = int(np.floor(M / 2))

  # Check signal length, and truncate if too long
  if max_len:
    if M > max_len:
      N = int(np.floor(max_len / 2))

  # G and d are (N-1) x (N-1)
  # d[k] = -signal[N+k] for k in range(1,N)
  d = -signal[N+1:2*N]

  try:
    from scipy.linalg import toeplitz, solve_toeplitz
    # Instead, form G = (c,r) as toeplitz
    #c = signal[N:2*N-1]
    #r = np.hstack((signal[1],signal[N-1:1:-1]))
    b = solve_toeplitz((signal[N:2*N-1],\
          np.hstack((signal[1],signal[N-1:1:-1]))),d,check_finite=False)
  except (ImportError,np.linalg.linalg.LinAlgError) as e:  
    # OLD CODE: sometimes more stable
    # G[k,m] = signal[N - m + k] for m,k in range(1,N)
    G = signal[N + np.arange(1,N)[:,None] - np.arange(1,N)]
    b = np.linalg.solve(G,d)

  # Now make b Nx1 where b0 = 1
  b = np.hstack((1,b))

  # b[m]*signal[k-m] for k in range(0,N), for m in range(k)
  a = np.dot(np.tril(toeplitz(signal[0:N])),b)
  p = np.poly1d(a)
  q = np.poly1d(b)

  # choose frequencies to evaluate over 
  frequency = np.arange(w_min,w_max,w_step)

  W = np.exp(-1j*frequency*stepsize)

  fsignal = p(W)/q(W)

  return fsignal, frequency 



# For some reason, the scipy FFT algo sucks for arrays of prime length
def is_prime(n):
  #"""https://stackoverflow.com/questions/18833759/"""
  if n % 2 == 0 and n > 2: 
    return False
  return all([bool(n % i) for i in range(3, int(math.sqrt(n)) + 1, 2)])



tds = []
def get_tds(filename):
  f = open(filename,'r')
  l = f.readline()
  while l:
    if l.split() != ["Singlet","state", "electronic", "transitions:"]: l = f.readline()
    else:
      #print("found tds")
      f.readline();f.readline() # headers not used
      f.readline();l=f.readline() # loads first row
      while l.split()[0] != "2":
        tds.append(l.split()[3:7])
        l = f.readline()
  f.close()

  E = []
def getE(filename):
  f = open(filename,'r')
  l = f.readline()
  while l:
    if l.split()!=['Root', 'Mult.', 'Total', 'Energy', '(a.u.)', 'Ex.', 'Energy', '(a.u.)', 'Ex.', 'Energy', '(eV)', 'Ex.', 'Energy', '(nm)']: l=f.readline()
    else:
      #print("found Es")
      l = l.split() #headers
      f.readline();l=f.readline().split()
  #rows.append([l[2], "0.0"]) # dont care about ground state
  l=f.readline()
  while l!="\n":
    E.append(l.split()[4])
    l = f.readline()
  f.close()

# RED = X, BLUE = Y, GREEN = Z
# format is (R,G,B), normalize to max of 1
def color_from_tds(td):
  tdf = [abs(float(td[0])),abs(float(td[1])),abs(float(td[2]))]
  #print(str(tdf)+"\n")
  #print(str(max(tdf))+"\n")
  if max(tdf) < 0.02: return((0,0,0)) # just in case
  return( (tdf[0]/max(tdf),tdf[2]/max(tdf),tdf[1]/max(tdf))  )

if options.tcoutfile:
  getE(options.tcoutfile)
  get_tds(options.tcoutfile)


def write_fft(infile):
  X = []
  Yr = []
  Yi = []
  with open(infile, 'r') as f:
    #pulse_ended = 0
    i = 0
    for line in f:
      #if options.pulse and pulse_ended == 0:
      #  if line == "======Pulse Ends=====\n": pulse_ended = 1
      #  else: continue
      line = line.split(',') 
      if len(line) != 3: pass # ignore blank or header lines
      elif i%options.nskip != 0: i+=1 # only include every Nth datapoint
      elif i > truncate: i+=1
      else:
        #X.append(np.float64(line[0])-1.0)
        X.append(np.float64(line[0]))
        Yr.append(np.float64(line[1]))
        Yi.append(np.float64(line[2]))
        i+=1


  #print("First few X: ", X[0],X[1],X[2])


  #Yr = np.array(Yr)
  #Yi = np.array(Yi)

  Y = np.vectorize(complex)(Yr, Yi)

  if options.center_mode == "s":
    X = np.array(X)-(X[0])

  if options.mirror:
    dif = X[1]-X[0]
    if options.center_mode == "u":
      X = np.array(list(np.linspace(-X[-1],-X[0],len(X))) + [0.0] + X)
      Y = np.array(list(np.conj(Y[::-1])) + list([np.complex(1.+0j)]) + list(Y)) 
    else:
      X = np.array(list(np.linspace(-X[-1],-X[1],len(X)-1)) + list(X))
      Y = np.array(list(np.conj(Y[1:][::-1])) + list(Y)) 
  else:
    if options.center_mode == "u":
      X.insert(0,0.0)
      Y = np.array( [np.complex(1.+0j)] + list(Y) )

  #print("First few X after shift/mirror: ", X[0],X[1],X[2])
  #print("len(X): "+str(len(X)))
  #print("len(Y) before primecheck: "+str(len(Y)))

  # Since prime length numpy.fft sucks so much,
  # just remove points until we get a non-prime length
  while is_prime(len(Y)):
    if options.mirror:
      Y = Y[1:-1]
      X = X[1:-1]
    else:
      Y = Y[:-1]
      X = X[:-1]

  #print("Len(Y) after primecheck: "+str(len(Y)))

  Yi = Y.imag
  Yr = Y.real

  W = []

  #print(options.windowfunc)

  if options.windowfunc == "hanning":
   N = len(Yr)
   for i in range(0,len(Yr)):
    W.append(0.5-0.5*np.cos(2*np.pi*i/N))
    Yr[i] = Yr[i]*(0.5 - 0.5*np.cos(2*np.pi*i/N))
    Yi[i] = Yi[i]*(0.5 - 0.5*np.cos(2*np.pi*i/N))

  if options.windowfunc == "hamming":
   N = len(Yr)
   for i in range(0,len(Yr)):
    W.append((0.54 - 0.46*np.cos(2*np.pi*i/N)))
    Yr[i] = Yr[i]*(0.54 - 0.46*np.cos(2*np.pi*i/N))
    Yi[i] = Yi[i]*(0.54 - 0.46*np.cos(2*np.pi*i/N))
    
  if options.windowfunc == "bh":
   N = len(Yr)
   for i in range(0,len(Yr)):
    W.append((0.355768-0.487396*np.cos(2*np.pi*i/N)+0.144232*np.cos(4*np.pi*i/N)-0.012604*np.cos(6*np.pi*i/N)))
    Yr[i] = Yr[i]*(0.355768-0.487396*np.cos(2*np.pi*i/N)+0.144232*np.cos(4*np.pi*i/N)-0.012604*np.cos(6*np.pi*i/N))
    Yi[i] = Yi[i]*(0.355768-0.487396*np.cos(2*np.pi*i/N)+0.144232*np.cos(4*np.pi*i/N)-0.012604*np.cos(6*np.pi*i/N))

  if options.windowfunc == "exp":
    N = len(Yr)
    tau = (N/2.)*(8.69/60)
    for i in range(0,len(Yr)):
      W.append(np.exp(-abs(i-((N-1)/2))/tau))
      Yr[i] = Yr[i]*(np.exp(-abs(i-((N-1)/2))/tau))
      Yi[i] = Yi[i]*(np.exp(-abs(i-((N-1)/2))/tau))


  windowf = open("window.csv",'w')
  for w in W: windowf.write(str(w)+"\n")
  windowf.close()

  delX = abs(X[4] - X[3])
  X = list(X)
  Yr = list(Yr)
  Yi = list(Yi)
  for _ in range(0, int(options.zeropad)):
    X.append( X[-1]+ delX )
    Yr.append(0.0)
    Yi.append(0.0)

  X = np.array(X)
  Yr = np.array(Yr)
  Yi = np.array(Yi)

  Y = np.vectorize(complex)(Yr, Yi)


  if options.rawproc:
    if options.outname:
      fraw = open(options.outname+"_rawproc.csv",'w')
    else:
      fraw = open(infile.split(".")[0]+"_rawproc.csv",'w')
    for j in range(0,len(X)): fraw.write(str(X[j])+","+str(Y[j].real)+","+str(Y[j].imag)+","+str(abs(Y[j]))+"\n")
    fraw.close()


  FFs = None
  #print(len(X))
  if options.use_pade:
    if not options.pstep:
      options.pstep = len(X)
    #import pdb; pdb.set_trace()
    FFs, freqs = pade(X, Y, sigma=options.psig, w_min=pmin*2*np.pi, w_max=pmax*2*np.pi, w_step=((float(pmax) - float(pmin))/(options.pstep))  )
    #print( (abs(FFs),freqs) )
    freqs = freqs/(2.*np.pi)

  elif options.fftw:
    import pyfftw
    FFs = pyfftw.interfaces.numpy_fft.fft(Y.real)
    freqs = np.fft.fftfreq(len(X), X[1]-X[0])
 
  else:
    FFs = np.fft.fft(Y.real)
    freqs = np.fft.fftfreq(len(X), X[1]-X[0])
  #print("finished fft")

  FFs = FFs/float(len(X)) # scale by number of DFFT bins, sort of normalization

  X_s=np.array(freqs)*10**18 # from inv attosecond to inv seconds
  X_eV = h*X_s # energy in eV
  #X_icm = X_s/c # energy in inverse centimeters

  if options.outname:
    f = open(options.outname+"_fft.csv",'w')
  else:
    f = open(infile.split(".")[0]+"_fft.csv",'w')
  for j in range(0,len(X_eV)): f.write(str(X_eV[j])+","+str(FFs[j].real)+","+str(FFs[j].imag)+","+str(abs(FFs[j]))+"\n")
  f.close()

  # Find range for plot
  #i = len(Y)-1
  #while Y[i] < options.window_thresh: i-=1
  #print( (len(X_eV), len(FFs)) )
  #print( (X_eV, abs(FFs)) )
  return (X_eV, FFs)
  
def import_from_csv(filename):
  X_eV = []
  FFs = []
  r = []
  i = []
  # makes flag interchangable with normal call
  if filename.split("/")[-1] == "TCorrFn": f = open(filename+"_fft.csv",'r')
  else: f=open(filename,'r')
  for line in f:
    linesplit = line.split(',')
    if len(linesplit) < 2:
      pass
      #print(linesplit)
    X_eV.append(float(linesplit[0]))
    r.append(float(linesplit[1]))
    i.append(float(linesplit[2]))
  FFs = np.vectorize(complex)(r,i)
  return (X_eV, FFs)


def drawgraph(X_ev, FFs,infile):
  # GRAPH RANGES
  xlim = [0.0,10]
  ylim = [-0.01,2]
  #aspectr = (xlim[1]-xlim[0])*9./((ylim[1]-ylim[0])*16)
  #aspectr = (xlim[1]-xlim[0])*0.45/((ylim[1]-ylim[0])*1)

  # not sure why but this doesnt work
  #if options.limstr:
  #  xlim = options.limstr.split(",")[0:2]
  #  ylim = options.limstr.split(",")[2:4]


  if options.xmin: xlim[0] = float(options.xmin)
  if options.xmax: xlim[1] = float(options.xmax)
  if options.ymin: ylim[0] = float(options.ymin)
  if options.ymax: ylim[1] = float(options.ymax)

  aspectr = (xlim[1]-xlim[0])*0.45/((ylim[1]-ylim[0])*1)

  plt.tight_layout()
  #plt.xticks(np.arange(xlim[0], xlim[1]+1, 1.0))
  #print(np.arange(xlim[0], xlim[1]+1, 1.0))

  fig = plt.figure()
  ax = fig.add_subplot(111,aspect=aspectr)
  #ax = fig.add_subplot(111)
  ax.set_ylabel('R')
  ax.set_xlabel('Energy (eV)')
  ax.set_ylim(ylim)
  ax.set_xlim(xlim)
  ax.xaxis.set_major_locator(plticker.MultipleLocator(1.0))
  # h/w = aspect
  #ax.set_aspect( aspectr)
  #ax.set_aspect()
  #ax.set_aspect(5)


  if options.tcoutfile:
    for i in range(0,len(E)):
      if float(tds[i][3]) >= td_thresh:
        #print((E[i],tds[i]))
        col=color_from_tds(tds[i])
        ax.axvline(x=float(E[i]), color=col, linewidth=0.5, label="v.ex. at {0:.2f} eV".format(float(E[i]), c=col))

  if options.ag3ndata:
    for e in (np.array([ 3.62, 4.83, 5.57, 5.80, 6.10, 6.3])-2.43 ):
      ax.axvline(x=e, color="black", linewidth=0.5, linestyle='--', label="v.ex. at {0:.2f} eV".format(e, c="black", linestyle='--'))

  try:
    if options.title.split()[-1] == "MRCIS":
      for e in [1.72196183080346, 1.73877592080593, 2.91311642895416, 3.28976082034876]:
        ax.axvline(x=e, color="black", linewidth=0.5, linestyle='--', label="v.ex. at {0:.2f} eV".format(e, c="black", linestyle='--'))
  except: pass


  if options.title: plt.title(options.title)
  ax.legend()
  #ax.set_aspect( (xlim[1]-xlim[0])*0.45/((ylim[1]-ylim[0])*1))
  ax.plot(X_eV, abs(FFs),linewidth=0.4,color="#000000")
  #ax.set_aspect(aspect=0.5)
  if options.outname:
    plt.savefig(options.outname+".png",dpi=900,aspect=0.5)
  else:
    plt.savefig(infile.split(".")[0]+".png",dpi=900,aspect=0.5)
  return 0


  
for datafile in args:
  X_eV = []
  FFs = [] 
  if options.no_fft:
    X_eV, FFs = import_from_csv(datafile)
  else:
    X_eV, FFs = write_fft(datafile)
  #print("max abs: "+str( max( abs(FFs[1:]))))
  if not options.no_plot:
    drawgraph(X_eV, FFs, datafile)


