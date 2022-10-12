#!/usr/bin/python2.7

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


import matplotlib.ticker as plticker
import numpy as np
import math
import sys
from optparse import OptionParser
usage = "usage: %prog [options] file1 file2 ..."
parser=OptionParser(usage=usage)
# General
parser.add_option("--outname", action="store", default=False, dest="outname", help="Set a name for associated output files like [outname]_fft.csv and [outname].png. Defaults to the input filename.")
parser.add_option("--floquet", action="store", default="No", dest="floquet", help="Input E in a.u. Will fft the sum of file1, file2*e^(-i*n*t), file3*e^(-i*2*n*t), etc, which are the input files shifted by multiples of E in frequency. (Use this for Floquet-TDCI individual photon correlation functions)")
#parser.add_option("--combine-ffts", action="store_true", default=False, dest="combine_ffts", help="file1, file2, file3, etc.'s output fft's will be added together for the output graph. If combined with --floquet [E], inputs will NOT be shifted before fft, instead they will be shifted afterwards. Expect 0-symmetry artifacts to show up at the shifted centers.")
#parser.add_option("--floquetscale", action="store", default=1.0, dest="floquetscale", help="The whole --floquet thing decreases output amplitudes by a lot. Not sure why. Use this to scale up the output. Default is 1.0 (no scaling)")


# Pre-processing
parser.add_option("--nskip", action="store", default=1, dest="nskip", help="only include every Nth datapoint in input")
parser.add_option("--truncate", action="store", default=999999999, dest="truncate", help="dont include data after Nth datapoint")

parser.add_option("-c", action="store", default="none", dest="center_mode", help="Centering mode. Options: s - shift the first datapoint to zero. u - add a (0,1) point. none - Don't center the data.")
parser.add_option("-m", "--mirror", action="store_true", default=False, dest="mirror", help="mirror data before FFT")
parser.add_option("-w", "--windowfunc", action="store", default=False, dest="windowfunc", help="Select windowing function ") # hanning, hamming, bh, exp
parser.add_option("--rawproc", action="store_true", default=False, dest="rawproc", help="Write a csv with pre-processed input data before FFT")
parser.add_option("-z", action="store", default=0, dest="zeropad", help="Pad this many zeros to the end of the data (applied after window)")

# Alternate FT algos
parser.add_option("--fftw", action="store_true", default=False, dest="fftw", help="Use pyfftw instead of numpy.fft.fft")
parser.add_option("--pade", action="store_true", default=False, dest="use_pade", help="Perform FT using method of Pade approximants.")
parser.add_option("--psig", action="store", default=False, dest="psig", help="'signal damp factor [windowing?], yield peaks with [sic: at maximum(?)] FWHM of 2/sigma' (???) Default: no signal dampening/windowing. Original code had default of 100.")
parser.add_option("--pmin", action="store", default=0.0,   dest="pmin", help="Pade: specify frequency calculation lower bound (in eV)")
parser.add_option("--pmax", action="store", default=15.0,  dest="pmax", help="Pade: specify frequency calculation upper bound (in eV)")
parser.add_option("--pstep", action="store", default=False, dest="pstep", help="Pade: specify number of points to generate in spectrum. Default: same as input data length")

# Post-FFT Processing
parser.add_option("--window-scale", action="store_true", default=False, dest="window_scaling", help="Scale output FT by the frequency response of the window function.")
parser.add_option("-u", "--arbitrary-units", action="store_true", default=False, dest="yunits", help="Scale Y to arbitrary units (max signal=1) ignores points within [-0.1 eV, 0.1 eV] to miss the zero-peak.")


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


parser.add_option("--xshift", action="store", default=False, dest="xshift", help="shift x values by constant")

# Analyses afterwards
parser.add_option("--popcheck", action="store_true", dest="popcheck", help="Integrate peaks at energies from terachem output, put them in a .csv with transition dipoles and populations from Pop file.")


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
#options.floquetscale = float(options.floquetscale)
if options.floquet != "No":
  options.floquet = float(options.floquet)
  print("f: "+str(options.floquet))

# pade options
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
    print("Toeplitz solver failed. Using general case solver for Pade-- careful not to run out of memory!")
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


def get_tds(filename):
  tds = []
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
  return tds

def getE(filename):
  E = []
  f = open(filename,'r')
  l = f.readline()
  while l:
    if l.split()!=['Root', 'Mult.', 'Total', 'Energy', '(a.u.)', 'Ex.', 'Energy', '(a.u.)', 'Ex.', 'Energy', '(eV)', 'Ex.', 'Energy', '(nm)']: l=f.readline()
    else:
      print("found Es")
      #l = l.split() #headers
      f.readline();l=f.readline()
      break
  #rows.append([l[2], "0.0"]) # dont care about ground state
  l=f.readline()
  while l!="\n":
    E.append(l.split()[4])
    l = f.readline()
  f.close()
  return E

# Used for color coding the direction of the transition dipole of states which -e represents as vertical lines
# RED = X, BLUE = Y, GREEN = Z
# format is (R,G,B), normalize to max of 1
def color_from_tds(td):
  tdf = [abs(float(td[0])),abs(float(td[1])),abs(float(td[2]))]
  #print(str(tdf)+"\n")
  #print(str(max(tdf))+"\n")
  if max(tdf) < 0.02: return((0,0,0)) # just in case
  return( (tdf[0]/max(tdf),tdf[2]/max(tdf),tdf[1]/max(tdf))  )


# We'll return the populations as averaged from step 3-13
# since populations change in steps 0,1 for delta pulse.
def getDeltaPop(filename):
  f = open(filename, 'r')
  f.readline();f.readline();f.readline() # get to third step
  n = 10
  float10 = lambda x: float(x)/n
  avg = map(float10, (f.readline()).split(',')[1:])
  for _ in range(0,n-1):
    avgt = map(float10, (f.readline()).split(',')[1:])
    avg = [sum(x) for x in zip(avg, avgt)]
  return avg
    

# inputs:
#  engs: list of energies in eV
#  fftX: the frequency values in eV of fft data
#  fftY: fft amplitudes
#
# Kinda crappy, but we average points that are less than 1/10th of the maximum to try to avoid peaks
# how does nmr plotting software auto-zero work?
def get_baseline(engs, fftX, fftY):
  baseline = 0j
  maxpeak = abs(max(fftY[5:]))
  lowerN = indexOfX(fftX, -50.0)
  upperN = indexOfX(fftX, 50.0)
  N = upperN - lowerN
  #n = len(fftX)-5
  print("maxpeak: "+str(maxpeak))
  print("init n: "+str(N))
  #for i in range(5, len(fftX)):
  for i in range(lowerN, upperN):
    #baseline += abs(fftY[i])
    if (abs(fftX[i]) > 0.3):  # don't integrate center peak into the zero, its huge.
      baseline += fftY[i]
    else:
      N = N-1
    #if (abs(fftY[i])< maxpeak/5):
    #  baseline += abs(fftY[i])
    #else:
    #  n = n-1
  print("end n: "+str(N))
  #baseline = baseline/(2.0*n)
  baseline = baseline/float(N)
  return baseline
  
# the same as list.index(x), but for ordered float lists where entries
# may not be exactly equal to x.
def indexOfX(fftX, x):
  i = 0
  while fftX[i]< x:
    i+=1
  # if the index on the greater-than side is closer to engs[i], use that one instead
  if ( abs(fftX[i]-x) > abs(fftX[i+1]-x)):
    i+=1
  return i

def IntegratePeaks(engs, fftX, fftY):
  baseline = get_baseline(engs, fftX, fftY)
  print("\nintegration baseline: "+str(baseline)+"\n")
  maxwidth = 0.4 # in eV
  integral = np.zeros(len(engs), dtype=np.cdouble)
  widths = np.zeros(len(engs), dtype=np.cdouble)
  dx = abs(fftX[10]-fftX[11])
  for i in range(0,len(engs)):
    peakj = indexOfX(fftX, engs[i])
    # set lower-bound for integration
    lowerbound = fftX[peakj] - maxwidth
    # don't overlap with other peaks, so go halfway between peaks
    if (i>0):
      if (abs(engs[i-1]-engs[i])/2.0 < maxwidth):
        lowerbound = engs[i] - abs(engs[i-1]-engs[i])/2.0
    if (lowerbound < dx):
      lowerbound = 2*dx # set to second point so we don't get the zero-peak
    # set upper-bound
    upperbound = fftX[peakj] + maxwidth
    if (i<len(engs)-1):
      if (abs(engs[i+1]-engs[i])/2.0 < maxwidth):
        upperbound = engs[i] + abs(engs[i+1]-engs[i])/2.0
    if (upperbound > fftX[len(fftX)-1]):
      upperbound = fftX[len(fftX)-3]
    # get indices of bounds
    lowerboundj = indexOfX(fftX, lowerbound)
    upperboundj = indexOfX(fftX, upperbound)
    #import pdb; pdb.set_trace()
    widths[i] = upperboundj - lowerboundj
    for j in range(lowerboundj, upperboundj+1):
      integral[i] += fftY[j]
      #if (abs(fftY[j]) >= baseline):
      #  integral[i] += abs(fftY[j])
    print("peak "+str(i)+"("+str(engs[i])+"eV): integral="+str(integral[i])+", lowerbound="+str(lowerbound)+
          " ("+str(lowerboundj)+"), upperbound="+str(upperbound)+" ("+str(upperboundj)+
          "), width="+str(widths[i])+", scaled="+str(abs(integral[i] - widths[i]*baseline))+"\n")
  return integral, widths, baseline
      
# Write a .csv that compares peak integrals to TDCI projected populations
def PopCheck(X_eV, FFs, E, tds):
  pops = getDeltaPop("Pop0")
  engs = []
  states = []
  for i in range(0, len(E)):
    if float(tds[i][3]) > 0.5:
      engs.append(float(E[i]))
      states.append(i)
  for i in range(0, len(tds)):
    for j in range(0,4):
      tds[i][j] = float(tds[i][j])

  print("X_eV: "+str(X_eV[0:10]))
  integrals, widths, baseline = IntegratePeaks(engs, X_eV, FFs)
  
  f = open("fft_analysis.csv",'w')
  f.write("State,Ex.Eng(eV),Tx,Ty,Tz,|T|,Pop,SumRe,SumIm,SumAbs\n")
  for i in range(0,len(states)):
    scaledsum = integrals[i] - widths[i]*baseline
    f.write(str(states[i]+1)+","+str(engs[i])+","+str(tds[states[i]][0])+","+
    str(tds[states[i]][1])+","+str(tds[states[i]][2])+","+str(tds[states[i]][3])+
    ","+str(pops[states[i]+1])+","+str(integrals[i].real)+","+str(integrals[i].imag)+","+str(abs(integrals[i]))+"\n" )
  f.close()


# SHIFTS EACH CORRFN BY e^(-i n*photoneng)
# Ylist should be ordered [0photons, 1photon, ...]
# In contrast to combine_FFT_outputs(), this function requires the same X range for all Y data.
def combine_corrfns(X, Ylist, photoneng):
  # sanity check
  for i in range(0,len(Ylist)):
    if len(X) != len(Ylist[i]):
      print(str(i)+"th element of Ylist has length "+str(len(Ylist[i]))+" which doesnt match the length of X, "+str(len(X)))
      return None
  #options.floquetscale = float(options.floquetscale)
  shiftY = [Ylist[0]]
  for i in range(1,len(Ylist)):
    # from hartree to J with as, kg*m^2/as^2
    # 1 H * (4.359744*10**-18 J / 1 H) *  1/(6.626*10**-34 J*s) * (1 s/ 10**18 as) = 1/as
    #import pdb; pdb.set_trace()
    photon_eng_ias = photoneng * 4.3597447222071*(10**-18) * 1./(6.62607015*10**-34) * (1./(10**18))
    # ugh X is in attoseconds so i have to do some kind of conversion
    #shifter = np.exp(-1j* i*photon_eng_ias*X)
    #shifter = options.floquetscale*np.exp(-1j*2*np.pi*i*photon_eng_ias*X)
    shifter = np.exp(-1j*2*np.pi*i*photon_eng_ias*X)
    # component-wise multiplication
    shiftY.append( Ylist[i]*shifter )
    #import pdb; pdb.set_trace()
    #shiftY.append( np.multiply(Ylist[i],shifter) )
    #shiftY.append( shifter )
  # combine shifted Y's
  outputY = sum(shiftY)
  print("type outputY: "+str(type(outputY)))
  return outputY

# SHIFTS FFT OUTPUTS AND THEN SUMS THEM
# THIS IS DUMB AND YOU SHOULD DO combine_corrnfs INSTEAD UNLESS YOU HAVE A GOOD REASON
# This function accepts separate X ranges, because if we shifted floquet data
# after fft, X range will be different, and we'll need to sort the Y datapoints into
# the nearest frequency bins of Xlist[0]
# Offset should be list of floats to shift the X of that output by.
def combine_FFT_output(Xlist, Ylist, offset):
  i = 1
  d = Xlist[0][1]-Xlist[0][0]
  outputY = Ylist[0][:]
  for i in range(0,len(Ylist)):
    for xi in range(0,len(Xlist[i])):
      # Closest index to offset X value
      # this requires that Xlists are linspaces that start at 0
      outputY[ int( (Xlist[i][xi]-offset[i])/d ) ] += Ylist[i][xi]
  return outputY


def import_corrfn(infile):
  X = []
  Yr = []
  Yi = []
  with open(infile, 'r') as f:
    i = 0
    for line in f:
      line = line.split(',') 
      if len(line) < 3: pass # ignore blank or header lines
      elif i%options.nskip != 0: i+=1 # only include every Nth datapoint
      elif i > truncate: i+=1
      else:
        X.append(np.float64(line[0]))
        Yr.append(np.float64(line[1]))
        Yi.append(np.float64(line[2]))
        i+=1
  Y = np.vectorize(complex)(Yr, Yi)
  X = np.array(X)
  return(X, Y)


def process_corrfn(X,Y):
  if options.center_mode == "s":
    X = np.array(X)-(X[0])

  if options.mirror:
    dif = X[1]-X[0]
    
    if options.center_mode == "u": # add a (0,1) point in the middle, mirror other points around it
      X = np.array(list(np.linspace(-X[-1],-X[0],len(X))) + [0.0] + X)
      Y = np.array(list(np.conj(Y[::-1])) + list([np.complex(1.+0j)]) + list(Y)) 
    else: # looks like this is broke, got X: ..., -2, 1, 2, ...
          # This should probably only be used with shift.
      X = np.array(list(np.linspace(-X[-1],-X[1],len(X)-1)) + list(X))
      Y = np.array(list(np.conj(Y[1:][::-1])) + list(Y)) 
    
  else: # no mirror
    if options.center_mode == "u":
      X.insert(0,0.0)
      Y = np.array( [np.complex(1.+0j)] + list(Y) )

  Yi = Y.imag
  Yr = Y.real

  W = []

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

  if options.zeropad:
    delX = abs(X[4] - X[3])
    X = list(X)
    Yr = list(Yr)
    Yi = list(Yi)
    for _ in range(0, int(options.zeropad)):
      X.append( X[-1]+ delX )
      Yr.append(0.0)
      Yi.append(0.0)
      W.append(0.0)
      if options.mirror:
        X.insert(0, X[0]- delX )
        Yr.insert(0, 0.0)
        Yi.insert(0, 0.0)
        W.insert(0, 0.0)

    X = np.array(X)
    Yr = np.array(Yr)
    Yi = np.array(Yi)

  Y = np.vectorize(complex)(Yr, Yi)

  #print("len(Y) before primecheck: "+str(len(Y)))
  # Since prime length numpy.fft sucks so much,
  # just remove points until we get a non-prime length
  while is_prime(len(Y)):
    print("Input was prime ("+str(len(Y))+", removing endpoint(s) to improve performance")
    if options.mirror:
      Y = Y[1:-1]
      X = X[1:-1]
      W = W[1:-1]
    else:
      Y = Y[:-1]
      X = X[:-1]
      W = W[:-1]

  return (X, Y, W)


# write_complex_csv( options.outname+"_rawproc.csv" , X, Y)
# X, regular array
# Y, array of complex numbers
def write_complex_csv(outfile, X, Y):
  f = open(outfile, 'w')
  for j in range(0,len(X)): f.write(str(X[j])+","+str(Y[j].real)+","+str(Y[j].imag)+","+str(abs(Y[j]))+"\n")
  f.close()
  return None


def calculate_fft(X, Y, infile):

  FFs = None
  W_FF = np.array([1])
  if options.use_pade:
    if not options.pstep:
      options.pstep = len(X)
    FFs, freqs = pade(X, Y, sigma=options.psig, w_min=pmin*2*np.pi, w_max=pmax*2*np.pi, w_step=((float(pmax) - float(pmin))/(options.pstep))  )
    freqs = freqs/(2.*np.pi)

  elif options.fftw:
    import pyfftw
    FFs = pyfftw.interfaces.numpy_fft.fft(Y)
    freqs = np.fft.fftfreq(len(X), X[1]-X[0])
    FFs = np.fft.fftshift(FFs)
    freqs = np.fft.fftshift(freqs)
    if ((options.window_scaling) and len(W) > 1):
      if len(W) != len(X): print("W and X dimension mismatch")
      W_FF = pyfftw.interfaces.numpy_fft.fft(np.array(W))
      FFs = FFs/W_FF
 
  else:
    print("y[1]: "+str(Y[1])+"\n")
    FFs = np.fft.fft(Y.real)
    freqs = np.fft.fftfreq(len(X), X[1]-X[0])
    FFs = np.fft.fftshift(FFs)
    freqs = np.fft.fftshift(freqs)
    if ((options.window_scaling) and len(W) > 1):
      if len(W) != len(X): print("W and X dimension mismatch")
      W_FF = np.fft.fft(np.array(W))
      FFs = FFs/W_FF

  print("Scaling by number of DFFT bins")
  FFs = FFs/(float(len(X))) # scale by number of DFFT bins, sort of normalization

  X_s=np.array(freqs)*10**18 # from inv attosecond to inv seconds
  X_eV = h*X_s # energy in eV
  #X_icm = X_s/c # energy in inverse centimeters

  if options.outname:
    f = open(options.outname+"_fft.csv",'w')
    if options.windowfunc and options.window_scaling:
      g = open(options.outname+"_window_fft.csv",'w')
  else:
    f = open(infile.split(".")[0]+"_fft.csv",'w')
    if options.windowfunc and options.window_scaling:
      g = open(infile.split(".")[0]+"_window_fft.csv",'w')
  for j in range(0,len(X_eV)):
    f.write(str(X_eV[j])+","+str(FFs[j].real)+","+str(FFs[j].imag)+","+str(abs(FFs[j]))+"\n")
  if options.windowfunc and options.window_scaling:
    for j in range(0,len(X_eV)):
      g.write(str(X_eV[j])+","+str(W_FF[j].real)+","+str(W_FF[j].imag)+","+str(abs(W_FF[j]))+"\n")
  f.close()
  if options.windowfunc and options.window_scaling:
    g.close()

  # Find range for plot
  #i = len(Y)-1
  #while Y[i] < options.window_thresh: i-=1
  #print( (len(X_eV), len(FFs)) )
  #print( (X_eV, abs(FFs)) )
  mmm = list(abs(FFs)[1:]).index(max(abs(FFs)[1:]))
  print("max: ")
  print( (X_eV[1:][mmm], abs(FFs)[1:][mmm]) )
  maxf = abs(FFs)[1:][mmm]
  return (X_eV, FFs, W_FF, maxf)




# for fft output data
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

def windowgraph(X_ev, W_FFs):
  xlim = [0.0,10]
  ylim = [-0.01,2]
  if options.xmin: xlim[0] = float(options.xmin)
  if options.xmax: xlim[1] = float(options.xmax)
  if options.ymin: ylim[0] = float(options.ymin)
  if options.ymax: ylim[1] = float(options.ymax)
  aspectr = (xlim[1]-xlim[0])*0.45/((ylim[1]-ylim[0])*1)
  plt.tight_layout()

  fig = plt.figure()
  ax = fig.add_subplot(111,aspect=aspectr)
  ax.set_ylabel('R (Arbitrary Units)')
  ax.set_xlabel('Energy (eV)')
  ax.set_ylim(ylim)
  ax.set_xlim(xlim)
  ax.xaxis.set_major_locator(plticker.MultipleLocator(1.0))
  # h/w = aspect
  #ax.set_aspect( aspectr)
  #ax.set_aspect()
  #ax.set_aspect(5)


  if options.title: plt.title(options.title+" window")
  #ax.legend()
  #ax.set_aspect( (xlim[1]-xlim[0])*0.45/((ylim[1]-ylim[0])*1))
  ax.plot(X_eV, abs(W_FFs),linewidth=0.4,color="#000000")
  #ax.set_aspect(aspect=0.5)
  if options.outname:
    plt.savefig(options.outname+"_window.png",dpi=900,aspect=0.5)
  else:
    plt.savefig(infile.split(".")[0]+"_window.png",dpi=900,aspect=0.5)
  return 0

def scale_to_arbitrary_units(X_eV, FFs):
  maxF = 0.0
  maxi = 0
  for i in range(0, len(X_eV)):
    if ((X_eV[i] > -0.1) and (X_eV[i] < 0.1)): # do not include center peak
      continue
    else:
      if abs(FFs[i]) > maxF:
        maxF = abs(FFs[i])
        maxi = i
  # max found
  for i in range(0, len(X_eV)):
    FFs[i] = FFs[i]/maxF
  # return the max, FFs array should be changed.
  print("scaling to arbitrary units:")
  print("maxi , X_eV[maxi]: "+str(maxi)+", "+str(X_eV[maxi]))
  print(str(maxF)+"->"+str(abs(FFs[maxi])))
  return FFs, maxi, maxF 


def drawgraph(X_eV, FFs,infile, maxf=None):
  # GRAPH RANGES
  xlim = [0.0,10]
  ylim = [-0.01,2]
  # x index of the max value (that isnt at x=0)
  mmm = list(abs(FFs)[1:]).index(max(abs(FFs)[1:]))
  maxf = abs(FFs)[1:][mmm]
  #aspectr = (xlim[1]-xlim[0])*9./((ylim[1]-ylim[0])*16)
  #aspectr = (xlim[1]-xlim[0])*0.45/((ylim[1]-ylim[0])*1)

  # not sure why but this doesnt work
  #if options.limstr:
  #  xlim = options.limstr.split(",")[0:2]
  #  ylim = options.limstr.split(",")[2:4]


  if options.xmin: xlim[0] = float(options.xmin)
  if options.xmax: xlim[1] = float(options.xmax)
  if options.ymin: ylim[0] = float(options.ymin)
  if options.ymax:
    ylim[1] = float(options.ymax)
  elif maxf != None:
    ylim[1] = float(maxf)

  aspectr = (xlim[1]-xlim[0])*0.45/((ylim[1]-ylim[0])*1)

  plt.tight_layout()
  #plt.xticks(np.arange(xlim[0], xlim[1]+1, (xlim[1]+1-xlim[0])/10.))
  #print(np.arange(xlim[0], xlim[1]+1, 1.0))

  fig = plt.figure()
  ax = fig.add_subplot(111,aspect=aspectr)
  #ax = fig.add_subplot(111)
  ax.set_ylabel('R')
  ax.set_xlabel('Energy (eV)')
  ax.set_ylim(ylim)
  ax.set_xlim(xlim)
  ax.xaxis.set_major_locator(plticker.MultipleLocator(0.05))
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

  if options.title: plt.title(options.title)
  ax.legend()
  #ax.set_aspect( (xlim[1]-xlim[0])*0.45/((ylim[1]-ylim[0])*1))
  ax.plot(X_eV, abs(FFs),linewidth=0.4,color="#000000")
  #ax.plot(X_eV, FFs.real,linewidth=0.4,color="#000000")
  #ax.plot(X_eV, FFs.real,linewidth=1.0,color="#000000")
  #ax.set_aspect(aspect=0.5)
  if options.outname:
    #plt.savefig(options.outname+".png",dpi=900,aspect=0.5)
    plt.savefig(options.outname+".png",dpi=900)
  else:
    #plt.savefig(infile.split(".")[0]+".png",dpi=900,aspect=0.5)
    plt.savefig(infile.split(".")[0]+".png",dpi=900)
  return 0

def main_addshift_ffts():
  print("Entering options.combine_ffts DUMB ADDSHIFT branch\n")
  FFslist = []
  Xlist = []
  for i in range(0,len(args)):
    #X, Y = import_corrfn(args[i])
    X, Y = import_corrfn(args[i])
    X_ev, Y, W = process_corrfn(X, Y)
    X_eV_temp, FFs_temp, W_FFs_temp, maxf = calculate_fft(X, Y, args[i])
    # plot individual FFT
    if not options.no_plot:
      filename = args[i].split(".")[0]+"_fft.csv"
      drawgraph(X_eV_temp, FFs_temp, filename, maxf)
      if (options.window_scaling and options.windowfunc):
        windowgraph(X_eV_temp, W_FFs_temp)

    # remove negative frequencies
    X_eVlist.append([])
    FFslist.append([])
    for j in range(0,len(X_eV_temp)):
      if (float(X_eV_temp[j]) < 0.0):
        continue
      else:
        X_eVlist[i].append(X_eV_temp[j])
        FFslist[i].append(X_eV_temp[j])
  Y = combine_FFT_output(X_eVlist, FFslist, np.zeros(len(FFslist)))
  if options.outname:
    filename = options.outname+"_combinedfft.csv"
  else: filename = args[0].split(".")[0]+"_combinedfft.csv"
  f = open(filename, 'w')
  f.write("X(eV), Yr, Yi\n")
  for j in range(0,len(X_eVlist[0])):
    f.write( str(X_eVlist[0][j])+","+str(Y[j].real)+","+str(Y[j].imag)+"\n")
  f.close();del f
  if options.outname:
    filename = options.outname+"_combinedfft.png"
  else: filename = args[0].split(".")[0]+"_combinedfft.png"
  print((Y[0], Y[1], Y[2]))
  drawgraph(np.array(X_eVlist[0]), np.array(Y), filename, maxf)

def main_mulshift_corrfns(E, tds):
  print("Entering options.floquet mulshift branch\n")
  X_eV = []
  FFs = []
  W_FFs = []
  Ylist = []
  X = []
  Y = []
  maxi = 0
  maxf = 0
  if options.no_fft:
    X_eV, FFs = import_from_csv(datafile)

  else: 
    for datafile in args:
      X, Y = import_corrfn(datafile)
      Ylist.append(Y)
    Y = combine_corrfns(X, Ylist, options.floquet)
    filename = args[0].split(".")[0]+"_combined.csv"
    if options.outname:  filename = str(options.outname)+"_combined.csv"
    f = open(filename, 'w')
    f.write("X(as), Yr, Yi\n")
    for j in range(0,len(X)):
      f.write( str(X[j])+","+str(Y[j].real)+","+str(Y[j].imag)+"\n")
    f.close();del f
    X_eV, Y, W = process_corrfn(X, Y)
    X_eV, FFs, W_FFs, maxf = calculate_fft(X_eV, Y, filename)
  #print("max abs: "+str( max( abs(FFs[1:]))))
  # Scale to arbitrary units
  if options.yunits:
    FFs, maxi, maxf = scale_to_arbitrary_units(X_eV, FFs)
  if not options.no_plot:
    print("before drawgraph: (maxi, maxf, FFs[maxi]): "+str((maxi,maxf, FFs[maxi])))
    drawgraph(X_eV, FFs, filename, maxf)
    if (options.window_scaling and options.windowfunc):
      windowgraph(X_eV, W_FFs)
  if options.tcoutfile and options.popcheck:
    PopCheck(X_eV, FFs, E, tds)

def main_regular(E, tds):
  #for datafile in args:
  X_eV = []
  FFs = [] 
  W_FFs = []
  maxf = None
  maxi = None
  if options.no_fft:
    #X_eV, FFs = import_from_csv(datafile)
    X_eV, FFs = import_from_csv(args[0])
  else:
    X, Y = import_corrfn(args[0])
    X_eV, Y, W = process_corrfn(X, Y)
    X_eV, FFs, W_FFs, maxf = calculate_fft(X_eV, Y, args[0])
  #print("max abs: "+str( max( abs(FFs[1:]))))
  # Scale to arbitrary units
  if options.yunits:
    FFs, maxi, maxf = scale_to_arbitrary_units(X_eV, FFs)
  if options.xshift:
    X_eV = np.array(X_eV) + float(options.xshift)
  if not options.no_plot:
    #drawgraph(X_eV, FFs, datafile, maxf)
    #print("before drawgraph: (maxi, maxf, FFs[maxi]): "+str((maxi,maxf, FFs[maxi])))
    drawgraph(X_eV, FFs, args[0], maxf)
    if (options.window_scaling and options.windowfunc):
      windowgraph(X_eV, W_FFs)
  #print("X_eV: "+str(X_eV[0:10]))
  if options.tcoutfile and options.popcheck:
    PopCheck(X_eV, FFs, E, tds)


# "main" below

E = []
tds = []
if options.tcoutfile:
  E = getE(options.tcoutfile)
  tds = get_tds(options.tcoutfile)
  #print("len(E): "+str(len(E))+", len(tds): "+str(len(tds))+"\n")
  #print(E)
  #print(tds)

# This is the dumb way to do floquet
# Adds FFT outputs together, shifted.
#if options.combine_ffts:
#  main_addshift_ffts()

# This is the right way to do floquet
# Multiples CorrFns by e^(-i shift), sums them together, then FFTs the whole thing.
if options.floquet != "No":
  main_mulshift_corrfns(E, tds)

else: # Regular ol' FFT
  main_regular(E, tds)


 
  
































