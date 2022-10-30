#!/bin/python
import numpy as np



#############################
# Author: Andy Durden ^_^
# Email: andys.durden@gmail.com
# Usage: python dihedral.py coords.xyz atom_i1 atom_i2 atom_i3 atom_i4
#
# Feed me an xyz file and the indices of 4 atoms, I'll spit out the dihedral
# 
#############################



import numpy as np


# just some basic vector operations for working on xyz files
# ONLY FOR 3-VECTORS
class vec:
  def __init__(self, x):
    if isinstance(x, vec):
      self.x = x.x
      return None
    elif len(x) != 3:
      print("ONLY 3-VECTORS!!")
      return 1
    self.x=np.array(x)

  # Vector Length
  def len(self):
    return np.sqrt(self.x[0]**2+self.x[1]**2+self.x[2]**2)
  def __len__(self):
    return np.sqrt(self.x[0]**2+self.x[1]**2+self.x[2]**2)
 
  def __getitem__(self, i):
    return self.x[i]

  # THIS IS CROSS PRODUCT, NOT DOT OR COMPONENT-WISE MULTIPLICATION
  # OUT = SELF * OTHER,   SELF ON LEFT OTHER ON RIGHT.
  def __mul__(self,other):
    # If the other instance isnt a vector, assume a scalar
    # and return the scalar product
    if not isinstance(other, vec):
      return vec(np.array((self.x[0]*other ,self.x[1]*other , self.x[2]*other )))
    else:
      return vec(np.array(( self.x[1]*other.x[2]-self.x[2]*other.x[1] ,
               self.x[2]*other.x[0]-self.x[0]*other.x[2] ,
               self.x[0]*other.x[1]-self.x[1]*other.x[0] )))

  # Scalar product, float * vec
  # should NOT EVER trigger in othervec * selfvec ( tested, it doesnt)
  def __rmul__(self, other):
    return vec(np.array((other*self.x[0] ,other*self.x[1] , other*self.x[2] )))
  
  def __eq__(self,other):
    if type(other) != type(self): return False
    elif other.__class__.__name__ != 'vec': return False
    if np.allclose(self.x,other.x): return True
    else: return False

  # Componentwise Subtraction
  # Note: when changing ref frames a-b goes FROM b to a

  def __sub__(self,other):
    return vec(( self.x[0]-other.x[0], 
             self.x[1]-other.x[1], 
             self.x[2]-other.x[2]))
  def __neg__(self):
    return vec([0,0,0])-self

  def __add__(self,other):
    return vec(( self.x[0]+other.x[0],
             self.x[1]+other.x[1],
             self.x[2]+other.x[2]))
  def __radd__(self,other): # its commutative
    return vec(( self.x[0]+other.x[0],
             self.x[1]+other.x[1],
             self.x[2]+other.x[2]))
  # Returns unit vector in direction of vec
  def norm(self):
    if np.allclose(self.x, np.array([0,0,0])):
      return self
    return vec(( self.x[0]/self.len(),
             self.x[1]/self.len(),
             self.x[2]/self.len()))

  def colinear(self,other):
    if self.norm()==other.norm() or self.norm()== -(other.norm()): return True
    else: return False

  def dot(self,other):
    return( self.x[0]*other.x[0]+self.x[1]*other.x[1]+self.x[2]*other.x[2])

def dist(p1,p2):
  return (vec(p2)-vec(p1)).len()

def angle(p1,p2,p3):
  a = vec(p2)-vec(p1)
  b = vec(p3)-vec(p2)
  return np.atan2( (a*b).len() , a.dot(b) )

# Returns dihedral angle (in radians) between four points
def dihedral(p1,p2,p3,p4):
  b1 = vec(p2)-vec(p1) # from p1 to p2
  b2 = vec(p3)-vec(p2)
  b3 = vec(p4)-vec(p3)
  n1 = (b1*b2).norm() # norm of plane that b1 and b2 are on
  n2 = (b2*b3).norm()
  # New axes: n1, b2.norm(), and n1*b2.norm()
  # Calculate the coordinates of n2 on these axes
  x = n1.dot(n2)
  y = (n1*b2.norm()).dot(n2)
  # NOTE: arctan2( 0.0, -1) = pi
  #       arctan2(-0.0, -1) = -pi
  # this is IEEE compliant
  # but ruins dihedrals under transformations
  # so lets force all zeros to be 'positive zeros'
  if x == 0: x=+0.0
  if y == 0: y=+0.0
  return np.arctan2(y,x) # Better behavior than arccos

def centroid(veclist):
  # cannot get sum(veclist) to work. runs into an int+vec. why???
  z = (0.,0.,0.)
  vecsum = vec(z)
  for veci in veclist:
    vecsum = vecsum + veci
  return vecsum*(1./float(len(veclist)))

import sys
# into our weird dictionary construct
def readxyz(filename):
  f = open(filename,'r')
  out = []
  f.readline();f.readline()
  currline = f.readline()
  while currline != "":
    out.append(currline.split())
    currline = f.readline()
  refx = []
  i = 0
  while i < len(out):
    #print(type(x[i][1]))
    refx.append({'refnum':i, 'z':out[i][0], 'x':np.array([float(out[i][1]),float(out[i][2]),float(out[i][3])])})
    i=i+1
  return refx


def writexyz(xyzdict, filename):
  atoms = map(lambda x: x['z'], xyzdict)
  coords = map(lambda x: x['x'], xyzdict)
  writexyz_(atoms, coords, filename)

def writexyz_(atoms, coords, filename):
  f = open(filename,'w')
  f.write(str(len(atoms))+'\n\n')
  for atom, coord in zip(atoms, coords):
      f.write(('{:>3s}'+'{:25.17f}'*3+'\n').format(atom, coord[0], coord[1], coord[2]))
  f.close()



def dihedral_main():
  refx = readxyz(sys.argv[2])
  a = int(sys.argv[2])-1
  b = int(sys.argv[3])-1
  c = int(sys.argv[4])-1
  d = int(sys.argv[5])-1
  d_rad = dihedral( refx[a]['x'], refx[b]['x'], refx[c]['x'], refx[d]['x'] ) 
  print( d_rad* 180.0/np.pi)
  return ((d_rad*180.0/np.pi))

# distance between two centroids
def centroid_main():
  # yeah im hardcoding the atom numbers so i dont have to pass 12 arguments what about it
  refx = readxyz(sys.argv[2])
  atoms1 = [1,2,3,4,5,6]
  atoms2 = [40,42,43,44,46]
  vec1 = list(map( lambda x: vec(refx[x]['x']), atoms1 ))
  vec2 = list(map( lambda x: vec(refx[x]['x']), atoms2 ))
  c1 = centroid(vec1)
  c2 = centroid(vec2)
  print(dist(c1,c2))
  return dist(c1,c2)

def dist_main():
  refx = readxyz(sys.argv[2])
  atom1 = 24
  atom2 = 53
  vec1 = vec(refx[atom1]['x'])
  vec2 = vec(refx[atom2]['x'])
  print(dist(vec1,vec2))
  return dist(vec1,vec2)


# k replaces K

import gcutil
def substitute_fragment(base_xyz, fragment_xyz, I, J, K, k, l, m, dist_Jk, angle_Jkl, dihedral_IJkl, dihedral_Jklm, base_atoms2remove=[]):
#def substitute_fragment(base_xyz, fragment_xyz, base_atoms2remove=[]):
  base_str, frag_str = "", ""
  with f as open(base_xyz, 'r'):
    base_str = f.read()
  with f as open(fragment_xyz, 'r'):
    frag_str = f.read()
  basel, fragl = base_str.split("\n"), frag_str.split("\n")
  while len(basel[-1]) < 4: basel.pop() # remove empty lines
  while len(fragl[-1]) < 4: fragl.pop()
  # Make sure our reference atoms are the last 3 in base
  Iline, Jline, Kline = basel[I+1], basel[J+1], basel[K+1]
  basel.remove(Iline); basel.remove(Jline); basel.remove(Kline)
  basel.append(Iline); basel.append(Jline); basel.append(Kline)
  I = len(basel)-3
  J, K = I+1, I+2
  # Make sure our reference atoms are the first 3 in frag
  kline, lline, mline = fragl[k+1], fragl[l+1], fragl[m+1]
  fragl.remove(kline); fragl.remove(lline); fragl.remove(mline)
  fragl.insert(0, mline); fragl.insert(0, lline); fragl.insert(0, kline)
  k, l, m = 0, 1, 2
  # Back to strings, then zmat
  base_str, frag_str = "\n".join(basel), "\n".join(fragl)
  base_zstr, frag_zstr = xyz_str2zmat(base_str), xyz_str2zmat(frag_zstr)
  basezl, fragzl = base_zstr.split("\n"), frag_zstr.split("\n")
  # Substitute K for k
  lastbase = basezl[-1].split()
  firstfrag = fragzl[0].split()
  lastbase[0] = firstfrag[0]
  if int(lastbase[1]) != J:
    print("WARNING: lastbase[1] != J ("+str(lastbase[1])+" != "+str(J)")")
  lastbase[2] = str(dist_Jk)
  basezl[-1] = "  ".join(lastbase)
  if len(fragzl) > 1:
    # Add l
    lline = fragzl[1].split()
    lline[1] = str(K)
    lline.extend([str(J), str(angle_Jkl), str(I) , str(dihedral_IJkl)])
    basezl.append( "  ".join(lline)  )
    if len(fragzl) >2:
      # Add m
      mline = fragzl[2].split()
      mline[1] = int(mline[1])+K-1
      mline[3] = int(mline[3])+K-1
      mline.extend([str(J), str(dihedral_Jklm)])
      basezl.append( "  ".join(mline) )
      if len(fragzl) >3:
        # Add any other lines
        for n in range(3,len(fragzl)):
          nline = fragzl[n].split(0
          nline[1] = int(nline[1])+K-1
          nline[3] = int(nline[3])+K-1
          nline[5] = int(nline[5])+K-1
	  basezl.append( "  ".join(nline) )
  basez_out = "\n".join(basezl)
  zoutf = open("test_out.zmat", 'w')
  zoutf.write(basez_out)
  zoutf.close()
  xyzout = zmat2xyz(basez_out)
  xyzoutf = open("test_out.xyz", 'w')
  xyzoutf.write(xyzout)
  xyzoutf.close()
  print("OUTPUT:\nZMAT:")
  print(basez_out)
  print("\nXYZ:")
  print(xyzout)
  

  
  
   
def xyz_str2zmat(xyz_str):
  xyzarr, atomnames = gcutil.readxyz(xyz_str)
  distmat = gcutil.distance_matrix(xyzarr)
  return gcutil.write_zmat(xyzarr, distmat, atomnames, False, False, False)
  

def xyz2zmat(xyzfile):
  xyzf = open(xyzfile, 'r')
  xyz_str = xyzf.read()
  xyzf.close()
  xyzarr, atomnames = gcutil.readxyz(xyz_str)
  distmat = gcutil.distance_matrix(xyzarr)
  return gcutil.write_zmat(xyzarr, distmat, atomnames, False, False, False)

def zmat2xyz(zmat_str):
  atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist = gcutil.readzmat(zmat_str)
  return gcutil.write_xyz(atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist)



# real main

if sys.argv[1] == "dihedral":
  dihedral_main()
elif sys.argv[1] == "centroid":
  centroid_main()
elif sys.argv[1] == "dist":
  dist_main()
elif sys.argv[1] == "zmat":
  substitute_fragment("temp.xyz","temp.xyz")
elif sys.argv[1] == "xyz":
  zmat2xyz(sys.argv[2])
















