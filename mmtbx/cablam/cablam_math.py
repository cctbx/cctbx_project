# (jEdit options) :folding=explicit:collapseFolds=1:
#The cablam_math module provides functions to calculate the backbone measures
#  used by the cablam system.
#The linked_residue.measures = {} dictionary will be empty unless functions from
# this module are called.
#One general note is that many of these calculations require atoms from many
#  different residues.  The first part of most functions deals with finding the
#  necessary atoms.  If a required atom (or residue) is missing, the function
#  will not add an entry to the measures dictionary in the current residue.  An
#  attempt to access the missing value will give a KeyError exception that must
#  be caught.
#Yes, the residue numbering is idiosyncratic: res2 is usually the res of
#  interest because 2 is its zero-indexed position in the CA pseudodihedral
#  calculations that define cablam-space.
#2012-03-02, This module has taken over the functions formerly in
#  cjw_vectormath.py  Also, results of these calculations will be stored in
#  linked_residue.measures instead of .results

#Note that all the geometry measures will currently break on alternate backbone
#  trace.  This will have to be fixed.  Plan: default to 'first_alt' in all
#  cases, develop special-use program to deal with alternates more completely.

import math
from cctbx import geometry_restraints

#{{{ vector math functions, formerly in cjw_vectormath.py
#rolling cjw_vectormath into cablam_math is a potential source of bugs

#Returns a vector connecting point p1 to point p2
def vectorize(p1, p2):
  v = [ p2[0]-p1[0] , p2[1]-p1[1] , p2[2]-p1[2] ]
  return v

#Returns the scalar length of a vector
def veclen(v):
  return math.sqrt( v[0]**2 + v[1]**2 + v[2]**2 )

#Dot product of two vectors
def dot(v1, v2):
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

#Cross product of two vectors
def cross(v1, v2):
  x = v1[1]*v2[2] - v1[2]*v2[1]
  y = v1[2]*v2[0] - v1[0]*v2[2]
  z = v1[0]*v2[1] - v1[1]*v2[0]
  return [x,y,z]

#Finds the line from a1 to a2, drops a perpendicular to it from b1, and returns
#  the point of intersection.
def perptersect(a1, a2, b1):
  #Find the slope of line A in each direction, A is in vector notation
  A = [a2[0]-a1[0], a2[1]-a1[1], a2[2]-a1[2]]
  #Solve the parametric equations (dot of perpendiculars=0). . .
  t = (A[0]*(b1[0]-a1[0]) + A[1]*(b1[1]-a1[1]) + A[2]*(b1[2]-a1[2])) / ((A[0]**2)+(A[1]**2)+(A[2]**2))
  # . . . and use the result to find the new point b2 on the line
  b2 = [a1[0]+A[0]*t, a1[1]+A[1]*t, a1[2]+A[2]*t]
  return b2

#Returns the perpendicular distance from point b1 to the a1-a2 line
def perpdist(a1, a2, b1):
  b2 = perptersect(a1, a2, b1)
  distance = veclen(vectorize(b1,b2))
  return distance

#}}}

#{{{ CA pseudodihedral/angle calculator
#Adds 'CA_d_in','CA_d_out', 'CA_a_in', 'CA_a', and 'CA_a_out' to
#  residue.measures={} for each residue in protein where protein is a dictionary
#  of cablam_res classes
#-------------------------------------------------------------------------------
def CApseudos(protein, dodihedrals = True, doangles = True):
  for resid2 in protein:
    res2 = protein[resid2]  #residue n
    gotall = False
    if res2.prevres:
      res1 = res2.prevres  # n-1
      if res1.prevres:
        res0 = res1.prevres  # n-2
        if res2.nextres:
          res3 = res2.nextres  # n+1
          if res3.nextres:
            res4 = res3.nextres  # n+2
            gotall = True
    if not gotall:
      continue

    CA_0 = res0.getatomxyz('CA')
    CA_1 = res1.getatomxyz('CA')
    CA_2 = res2.getatomxyz('CA')
    CA_3 = res3.getatomxyz('CA')
    CA_4 = res4.getatomxyz('CA')
    if None in [CA_0,CA_1,CA_2,CA_3,CA_4]:
      continue

    if dodihedrals:
      d_in = geometry_restraints.dihedral(sites=[CA_0,CA_1,CA_2,CA_3],
        angle_ideal=-40, weight=1)
      d_out = geometry_restraints.dihedral(sites=[CA_1,CA_2,CA_3,CA_4],
        angle_ideal=-40, weight=1)

      res2.measures['CA_d_in']  = d_in.angle_model
      res2.measures['CA_d_out'] = d_out.angle_model

    if doangles:
      a_in  = geometry_restraints.angle(sites=[CA_0,CA_1,CA_2],
        angle_ideal=120, weight=1)
      a     = geometry_restraints.angle(sites=[CA_1,CA_2,CA_3],
        angle_ideal=120, weight=1)
      a_out = geometry_restraints.angle(sites=[CA_2,CA_3,CA_4],
        angle_ideal=120, weight=1)
      res2.measures['CA_a_in']  = a_in.angle_model
      res2.measures['CA_a']     = a.angle_model
      res2.measures['CA_a_out'] = a_out.angle_model
#-------------------------------------------------------------------------------
#}}}

#{{{ cablam measures calculator (CA_d_in, CA_d_out, CA_a)
#Adds 'CA_d_in', 'CA_d_out', and 'CA_a' only to residue.measures={} for each
#  residue in portein where protein is a dictionary of cablem_res classes.
#This function is a condensed version of CApseudos() above, returning only the
#  measures currently in use by cablam annotation.
def cablam_measures(protein):
  for resid2 in protein:
    res2 = protein[resid2]  #residue n
    gotall = False
    if res2.prevres:
      res1 = res2.prevres  # n-1
      if res1.prevres:
        res0 = res1.prevres  # n-2
        if res2.nextres:
          res3 = res2.nextres  # n+1
          if res3.nextres:
            res4 = res3.nextres  # n+2
            gotall = True
    if not gotall:
      continue

    CA_0 = res0.getatomxyz('CA')
    CA_1 = res1.getatomxyz('CA')
    CA_2 = res2.getatomxyz('CA')
    CA_3 = res3.getatomxyz('CA')
    CA_4 = res4.getatomxyz('CA')
    if None in [CA_0,CA_1,CA_2,CA_3,CA_4]:
      #getatomxyz returns either an xyz list (tuple?) or None
      continue

    d_in = geometry_restraints.dihedral(sites=[CA_0,CA_1,CA_2,CA_3],
      angle_ideal=-40, weight=1)
    d_out = geometry_restraints.dihedral(sites=[CA_1,CA_2,CA_3,CA_4],
      angle_ideal=-40, weight=1)
    a     = geometry_restraints.angle(sites=[CA_1,CA_2,CA_3],
      angle_ideal=120, weight=1)

    res2.measures['CA_d_in']  = d_in.angle_model
    res2.measures['CA_d_out'] = d_out.angle_model
    res2.measures['CA_a']     = a.angle_model
#}}}

#{{{ CO pseudodihedral calculator
#Adds 'CO_d_in' and 'CA_O_out' dihedrals to residue.measures={} for each residue
#  in protein where protein is a dictionary of cablam_res classes
#-------------------------------------------------------------------------------
def COpseudodihedrals(protein):
  for resid2 in protein:
    res2 = protein[resid2]        #residue n
    gotall = False
    if res2.prevres:
      res1 = res2.prevres  # n-1
      if res2.nextres:
        res3 = res2.nextres  # n+1
        if res3.nextres:
          res4 = res3.nextres  # n+2
          gotall = True
    if not gotall:
      continue

    CA_1, O_1 = res1.getatomxyz('CA'),res1.getatomxyz('O')
    CA_2, O_2 = res2.getatomxyz('CA'),res2.getatomxyz('O')
    CA_3, O_3 = res3.getatomxyz('CA'),res3.getatomxyz('O')
    CA_4      = res4.getatomxyz('CA')
    if None in [CA_1,CA_2,CA_3,CA_4,O_1,O_2,O_3]:
      continue

    pseudoC_1 = perptersect(CA_1,CA_2,O_1)
    pseudoC_2 = perptersect(CA_2,CA_3,O_2)
    pseudoC_3 = perptersect(CA_3,CA_4,O_3)

    d_in = geometry_restraints.dihedral(sites=[O_1, pseudoC_1, pseudoC_2, O_2],
      angle_ideal=-40, weight=1)
    d_out = geometry_restraints.dihedral(sites=[O_2, pseudoC_2, pseudoC_3, O_3],
      angle_ideal=-40, weight=1)

    res2.measures['CO_d_in']  = d_in.angle_model
    res2.measures['CO_d_out'] = d_out.angle_model
#-------------------------------------------------------------------------------
#}}}

#{{{ Ramachandran calculator
#Adds 'phi' 'psi' 'phi+1' and 'psi-1' dehedrals  to residue.measures={} for each
#  residue in protein where protein is a dictionary of cablam_res classes
#-------------------------------------------------------------------------------
def phipsi(protein):
  for resid2 in protein:
    res2 = protein[resid2]        #residue n
    gotall = False
    if res2.prevres:
      res1 = res2.prevres  # n-1
      if res2.nextres:
        res3 = res2.nextres  # n+1
        gotall = True
    if not gotall:
      continue

    CO_1     = res1.getatomxyz('C')
    N_2      = res2.getatomxyz('N')
    CA_2,C_2 = res2.getatomxyz('CA'),res2.getatomxyz('C')
    N_3      = res3.getatomxyz('N')
    if None in [CO_1,N_2,CA_2,C_2,N_3]:
      continue

    phi = geometry_restraints.dihedral(sites=[CO_1, N_2, CA_2, C_2],
      angle_ideal=-40, weight=1)
    psi = geometry_restraints.dihedral(sites=[N_2, CA_2, C_2, N_3],
      angle_ideal=-40, weight=1)

    res2.measures['phi'] = phi.angle_model
    res2.measures['psi'] = psi.angle_model

    res1.measures['phi+1'] = phi.angle_model
    res3.measures['psi-1'] = psi.angle_model
#-------------------------------------------------------------------------------
#}}}

#{{{ tau angle calculator
#Adds 'tau' angle (N-CA-C) to residue.measures={} for each residue in protein
#  where protein is a dictionary of cablam_res classes
#-------------------------------------------------------------------------------
def taucalc(protein):
  for resid2 in protein:
    res2 = protein[resid2]

    N  = res2.getatomxyz('N')
    CA = res2.getatomxyz('CA')
    C  = res2.getatomxyz('C')
    if None in [N,CA,C]:
      continue

    tau = geometry_restraints.angle(sites=[N,CA,C], angle_ideal=120, weight=1)
    res2.measures['tau']  = tau.angle_model
#-------------------------------------------------------------------------------
#}}}

#{{{ omega angle calculator
#Adds 'omega' dihedral (peptide plane dihedral) to residue.measures={} for each
#  residue in protein where protein is a dictionary of cablam_res classes
#-------------------------------------------------------------------------------
def omegacalc(protein):
  for resid2 in protein:
    res2 = protein[resid2]
    gotall = False
    if res2.nextres:
      res3 = res2.nextres #n+1
      gotall = True
    if not gotall:
      continue

    CA  = res2.getatomxyz('CA')
    C   = res2.getatomxyz('C')
    N_2 = res3.getatomxyz('N')
    CA_2= res3.getatomxyz('CA')
    if None in [CA, C, N_2, CA_2]:
      #print res2.resnum, CA, C, N_2, CA_2
      continue

    omega = geometry_restraints.dihedral(sites=[CA,C,N_2,CA_2],
      angle_ideal=-40, weight=1)

    res2.measures['omega'] = omega.angle_model
#-------------------------------------------------------------------------------
