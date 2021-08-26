#!/bin/python
<<<<<<< HEAD
import numpy as np




# just some basic vector operations for working on xyz files

=======

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
>>>>>>> 789a18831bf8891012a8c882da4793eeca6ede16
# ONLY FOR 3-VECTORS
class vec:
  def __init__(self, x):
    if len(x) != 3:
      print("ONLY 3-VECTORS!!")
      return 1
    self.x=np.array(x)

  # Vector Length
  def len(self):
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
# 
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



import sys
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


refx = readxyz(sys.argv[1])
a = int(sys.argv[2])-1
b = int(sys.argv[3])-1
c = int(sys.argv[4])-1
d = int(sys.argv[5])-1
#for i in [a,b,c,d]: print(refx[i])


d_rad = dihedral( refx[a]['x'], refx[b]['x'], refx[c]['x'], refx[d]['x'] ) 
print( d_rad* 180.0/np.pi)






