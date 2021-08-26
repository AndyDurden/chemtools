
# includes tddft
class cis:
  def __init__(self, filename):
    self.filename = filename


  def get_gs_dipole(self,filename):
    f = open(filename, 'r')
    l = f.readline()
    while l:
      if l.split()[:2] != ['DIPOLE', 'MOMENT:']: l=f.readline()
      else:
        x = float(l.split()[2].lstrip("{").rstrip(","))
        y = float(l.split()[3].rstrip(","))
        z = float(l.split()[4].rstrip("}"))
        return ( (x,y,z) )

  def get_state_dipoles(self,filename):
    dipoles = [self.get_gs_dipole(filename)]
    f = open(filename, 'r')
    l = f.readline()
    while l:
      if l.split() != ["Unrelaxed", "excited", "state", "dipole", "moments:"]: l = f.readline()
      else:
        f.readline();f.readline()
        f.readline();l=f.readline()
        while len(l.split()) == 5:
          ls = map(float, l.split())
          dipoles.append( (ls[1], ls[2], ls[3]) )
          l=f.readline()
    f.close()
    return dipoles


  def get_tds(self,filename):
    tds = []
    f = open(filename,'r')
    l = f.readline()
    gslines = []
    eslines = []
    while l:
      if l.split() != ["Transition","dipole", "moments:"]: l = f.readline()
      else:
        print("found gs tds")
        f.readline();f.readline() # headers not used
        f.readline();l=f.readline() # loads first row
        i=1
        while len(l.split()) == 5:
          gslines.append(l)
          l=f.readline()

    f.close()
    g = open(filename, 'r')
    l = g.readline()
    while l:
      if l.split() != ["Transition", "dipole", "moments", "between","excited","states:"]: l=g.readline()
      else:
        print("found es tds")
        g.readline();g.readline() # headers not used
        g.readline();l=g.readline() # loads first row
        i=1
        while len(l.split()) == 7:
          eslines.append(l)
          l=g.readline()

    N = int(gslines[-1].split()[0]) +1
    for i in range(0,N):
      tds.append([])
      for j in range(0,N):
        tds[i].append( (0.,0.,0.) )

    for line in gslines:
      i = int(line.split()[0])
      tds[0][i] = tds[i][0] = tuple(map(float, line.split()[1:4]))

    print("eslines: ")
    print(eslines)

    for line in eslines:
      i = int(line.split()[0])
      j = int(line.split()[2])
      tds[i][j] = tds[j][i] = tuple(map(float, line.split()[3:6]))

    return tds


  def get_E(self,filename):
    E = [0.0]
    f = open(filename,'r')
    l = f.readline()
    while l:
      if l.split()!=['Root', 'Total', 'Energy', '(a.u.)', 'Ex.', 'Energy', '(a.u.)', 'Ex.', 'Energy', '(eV)', 'Ex.', 'Energy', '(nm)']: l=f.readline()
      else:
        print("found Es")
        #l = l.split() #headers
        f.readline()
        break
    #rows.append([l[2], "0.0"]) # dont care about ground state
    l=f.readline()
    while l!="\n":
      E.append(float(l.split()[2]))
      l = f.readline()
    f.close()
    return E



class casci:
  def __init__(self, filename):
    self.filename = filename

  def get_state_dipoles(self,filename):
    dipoles = []
    f = open(filename, 'r')
    l = f.readline()
    while l:
      if l.split() != ["Singlet", "state", "dipole", "moments:"]: l = f.readline()
      else:
        f.readline();f.readline()
        f.readline();l=f.readline()
        while len(l.split()) == 5:
          ls = map(float, l.split())
          dipoles.append( (ls[1], ls[2], ls[3]) )
    f.close()
    return dipoles

  def get_tds(self,filename):
    tds = []
    f = open(filename,'r')
    l = f.readline()
    lines = []
    while l:
      if l.split() != ["Singlet","state", "electronic", "transitions:"]: l = f.readline()
      else:
        print("found tds")
        f.readline();f.readline() # headers not used
        f.readline();l=f.readline() # loads first row
        i=1
        while len(l.split()) == 8:
          lines.append(l)

    print(lines)
    print(lines[-1].split())
    N = int(lines[-1].split()[2])
    for i in range(0,N):
      tds.append([])
      for j in range(0,N):
        tds[i].append( (0.,0.,0.) )
    for line in lines:
      i = int(line.split()[0])
      j = int(line.split()[2])
      tds[i][j] = tds[j][i] = tuple(map(float, line.split()[3:6]))

    f.close()
    return tds

  def get_E(self,filename):
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




import math
import numpy
norm = numpy.linalg.norm
def make_GMH(filename, fileclass):
  gmh = []
  d = fileclass.get_state_dipoles(filename)
  td = fileclass.get_tds(filename)
  E = fileclass.get_E(filename)
  print("td ("+str(len(td))+"): ");print(td)
  print("E ("+str(len(E))+"): ");print(E)
  print("d ("+str(len(d))+"): ");print(d)
  print("lens: "+str( (len(d), len(td), len(E))))
  for i in range(0,len(td)):
    gmh.append([])
    for j in range(0,len(td)):
      #gmhij =  (E[j] - E[i])*td[i][j] / math.sqrt( (d[i]-d[j])**2 + 4*(td[i][j]**2) )
      gmhij =  (E[j] - E[i])*norm(td[i][j]) / math.sqrt( norm(vecdiff(d[i],d[j]))**2 + 4*(norm(td[i][j])**2) )
      gmh[i].append( gmhij )
  return gmh

def vecdiff(a,b):
  out = []
  for i in range(0,len(a)):
    out.append(a[i]-b[i])
  return tuple(out)


def matrix2csv(matrix,outname):
  f = open(outname,"w")
  for row in matrix:
    for e in row:
      f.write(str(e)+",")
    f.write("\n")
  f.close()






import sys
fileclass = cis(sys.argv[1])
gmh = make_GMH(sys.argv[1], fileclass)
#print(gmh)

matrix2csv(gmh, sys.argv[1]+"_gmh.csv")




