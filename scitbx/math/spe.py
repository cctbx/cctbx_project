from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from scitbx.stdlib import math
from six.moves import range

"""
This is a python implementation of the
Stochastic Proximity Embedding Algorithm
(J.Chem. Information and Modelling, 2011,51,2852-2859.

The algorithm implemented below is the
non-pivotted one. At the moment, we don't have any convergence criteria implemented
and run the routine for a fixed number of cycles only.

The input needed is a similarity maxtrix that looks like this

[ a0,a1,a2,a3,.....,an  ]

ai = [ (1,di1), (4,di4), ..., (j,dij) ]

the object ai is a list of tuples. The tuple has two items
- j
- dij
These two value convey the the message: object i has a distance of dij to object j.

The algorithm takes these distance constraint and embeds them in an object in N-dimensional space.

One could use this for very rough conformation generation from distance constraints from a molecule,
or for more classic multidimensional scaling type problems.

"""

class classic_spe_engine(object):
  def __init__(self, dmat,l=.850,max_cycle=10000):
    self.dmat = dmat
    self.l = l
    self.max_cycle = max_cycle
    self.dl = self.l/max_cycle
    self.eps=1e-12

  def embed(self,n_dimensions,n_points):
    x = []
    for ii in range(n_points):
      x.append( flex.random_double(n_dimensions)*100 )

    l = float(self.l)
    for mm in range(self.max_cycle):
      atom_order = flex.sort_permutation( flex.random_double(len(x))  )
      strain = 0.0
      for ii in atom_order:
        n_contacts = len(self.dmat[ii])
        jj_index = flex.sort_permutation( flex.random_double( n_contacts ) )[0]
        jj_info = self.dmat[ii][jj_index]
        jj = jj_info[0]
        td = jj_info[1]
        xi = x[ii]
        xj = x[jj]
        cd = math.sqrt( flex.sum( (xi-xj)*(xi-xj) ) )
        new_xi = xi + l*0.5*(td-cd)/(cd+self.eps)*(xi-xj)
        new_xj = xj + l*0.5*(td-cd)/(cd+self.eps)*(xj-xi)
        strain += abs(cd-td)
        x[ii] = new_xi
        x[jj] = new_xj
      l = l-self.dl
    return x,strain/len(x)


def tst():
  M=10
  xy = []
  for ii in range(M):
    r=10.0
    phi = ii*(math.pi*2/M)
    x = r*math.cos(phi)
    y = r*math.sin(phi)
    xy.append( flex.double([x,y]) )
  # build distance matrix
  dmat = []
  for ii in range(M):
    tmp = []
    for jj in range(M):

      x1=xy[ii][0]
      x2=xy[jj][0]

      y1=xy[ii][1]
      y2=xy[jj][1]

      dd = math.sqrt( (x1-x2)**2.0 +(y1-y2)**2.0  )
      if jj != ii:
          tmp.append( (jj,dd) )
    dmat.append( tmp )
  spee=classic_spe_engine( dmat, l=1.0,max_cycle=5000 )
  x,s = spee.embed(2,M)
  assert s<1e-4
  print("OK")

if __name__ == "__main__":
  tst()
