from __future__ import absolute_import, division, print_function
from six.moves import range
from scitbx import lbfgs
from libtbx import adopt_init_args

from cctbx.array_family import flex

class NoCrystalSystem(Exception): pass

class s_minimizer:

  def __init__(self, orient, constraint='triclinic',
               min_iterations=25, max_calls=1000):
    self.constraint=constraint
    adopt_init_args(self, locals())
    self.n = 9
    self.x = flex.double(orient.direct_matrix())
    self.minimizer = lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs.termination_parameters(
        traditional_convergence_test=00000,
        min_iterations=min_iterations,
        max_calls=max_calls))
    del self.g

  def compute_functional_and_gradients(self):
    self.internals()
    return self.f,self.g

  def internals(self):
    O = self.newOrientation()
    A = O.A
    B = O.B
    C = O.C
    D = O.D
    E = O.E
    F = O.F
    Aij = self.x
    dA_dAij = (2.*Aij[0],2.*Aij[1],2.*Aij[2],0,0,0,0,0,0)
    dB_dAij = (0,0,0,2.*Aij[3],2.*Aij[4],2.*Aij[5],0,0,0)
    dC_dAij = (0,0,0,0,0,0,2.*Aij[6],2.*Aij[7],2.*Aij[8])
    dD_dAij = (0,0,0,Aij[6],Aij[7],Aij[8],Aij[3],Aij[4],Aij[5])
    dE_dAij = (Aij[6],Aij[7],Aij[8],0,0,0,Aij[0],Aij[1],Aij[2])
    dF_dAij = (Aij[3],Aij[4],Aij[5],Aij[0],Aij[1],Aij[2],0,0,0)
    self.g = flex.double()

    if self.constraint == 'monoclinic':
      self.f = D*D + F*F
      for x in range(self.n):
        self.g.append(2.*D*dD_dAij[x] + 2.*F*dF_dAij[x])

    elif self.constraint == 'orthorhombic':
      self.f = D*D + E*E + F*F
      for x in range(self.n):
        self.g.append(2.*D*dD_dAij[x] + 2.*E*dE_dAij[x] + 2.*F*dF_dAij[x])

    elif self.constraint == 'tetragonal':
      self.f = D*D + E*E + F*F + (A-B)*(A-B)
      for x in range(self.n):
        self.g.append(2.*D*dD_dAij[x] + 2.*E*dE_dAij[x] + 2.*F*dF_dAij[x] +
                      2.*(A-B)*(dA_dAij[x]-dB_dAij[x]))

    elif self.constraint == 'cubic':
      self.f = D*D + E*E + F*F + (A-B)*(A-B) + (A-C)*(A-C)
      for x in range(self.n):
        self.g.append(2.*D*dD_dAij[x] + 2.*E*dE_dAij[x] + 2.*F*dF_dAij[x] +
                      2.*(A-B)*(dA_dAij[x]-dB_dAij[x]) +
                      2.*(A-C)*(dA_dAij[x]-dC_dAij[x]))

    elif self.constraint == 'rhombohedral':
      self.f = (A-B)*(A-B) + D*D + E*E + (F+A/2.)*(F+A/2.)
      for x in range(self.n):
        self.g.append(2.*(A-B)*(dA_dAij[x]-dB_dAij[x]) +
                      2.*D*dD_dAij[x] + 2.*E*dE_dAij[x] +
                      2.*(F+A/2.)*(dF_dAij[x]+dA_dAij[x]/2.) )

    elif self.constraint == 'hexagonal':
      self.f = D*D + E*E + (A-B)*(A-B) + (F+A/2.)*(F+A/2.)
      for x in range(self.n):
        self.g.append(2.*D*dD_dAij[x] + 2.*E*dE_dAij[x] +
                      2.*(A-B)*(dA_dAij[x]-dB_dAij[x])  +
                      2.*(F+A/2.)*(dF_dAij[x]+dA_dAij[x]/2.))

    elif self.constraint == 'triclinic':
      self.f = 0.0
      for x in range(self.n):
        self.g.append(0.0)

    else:
      raise NoCrystalSystem

  def newOrientation(self):
    #trick to instantiate Orientation given self.x, the direct space matrix
    from rstbx.dps_core import Orientation
    from cctbx.crystal_orientation import basis_type
    return Orientation(tuple(self.x),basis_type.direct)
