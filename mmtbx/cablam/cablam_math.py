# (jEdit options) :folding=explicit:collapseFolds=1:
#The cablam_math module provides functions to calculate the backbone measures
#  used by the cablam system.
#The ResConnect.results= {} dictionary will be empty unless functions from this
#  module are called.
#One general note is that many of these calculations require atoms from many
#  different residues.  The first part of most functions deals with finding the
#  necessary atoms.  If a required atom (or residue) is missing, the function
#  will not add an entry to the results dictionary in the current residue.  An
#  attempt to access the missing value will give a KeyError exception that must
#  be caught.
#Yes, the residue numbering is idiosyncratic

import cjw_vectormath
from cctbx import geometry_restraints

#{{{ CA pseudodihedral/angle calculator
#Adds 'CA_d_in','CA_d_out', 'CA_a_in', 'CA_a', and 'CA_a_out' to
#  residue.results={} for each residue in protein where protein is a dictionary
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

    try:
      CA_0 = res0.atomxyz['']['CA']
      CA_1 = res1.atomxyz['']['CA']
      CA_2 = res2.atomxyz['']['CA']
      CA_3 = res3.atomxyz['']['CA']
      CA_4 = res4.atomxyz['']['CA']
    except KeyError:
      continue

    if dodihedrals:
      d_in = geometry_restraints.dihedral(sites=[CA_0,CA_1,CA_2,CA_3],
        angle_ideal=-40, weight=1)
      d_out = geometry_restraints.dihedral(sites=[CA_1,CA_2,CA_3,CA_4],
        angle_ideal=-40, weight=1)

      res2.results['CA_d_in']  = d_in.angle_model
      res2.results['CA_d_out'] = d_out.angle_model

    if doangles:
      a_in  = geometry_restraints.angle(sites=[CA_0,CA_1,CA_2],
        angle_ideal=120, weight=1)
      a     = geometry_restraints.angle(sites=[CA_1,CA_2,CA_3],
        angle_ideal=120, weight=1)
      a_out = geometry_restraints.angle(sites=[CA_2,CA_3,CA_4],
        angle_ideal=120, weight=1)
      res2.results['CA_a_in']  = a_in.angle_model
      res2.results['CA_a']     = a.angle_model
      res2.results['CA_a_out'] = a_out.angle_model
#-------------------------------------------------------------------------------
#}}}

#{{{ CO pseudodihedral calculator
#Adds 'CO_d_in' and 'CA_O_out' to residue.results={} for each residue in
#  protein where protein is a dictionary of cablam_res classes
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

    try:
      CA_1, O_1 = res1.atomxyz['']['CA'],res1.atomxyz['']['O']
      CA_2, O_2 = res2.atomxyz['']['CA'],res2.atomxyz['']['O']
      CA_3, O_3 = res3.atomxyz['']['CA'],res3.atomxyz['']['O']
      CA_4      = res4.atomxyz['']['CA']
    except KeyError:
      continue

    pseudoC_1 = cjw_vectormath.perptersect(CA_1,CA_2,O_1)
    pseudoC_2 = cjw_vectormath.perptersect(CA_2,CA_3,O_2)
    pseudoC_3 = cjw_vectormath.perptersect(CA_3,CA_4,O_3)

    d_in = geometry_restraints.dihedral(sites=[O_1, pseudoC_1, pseudoC_2, O_2],
      angle_ideal=-40, weight=1)
    d_out = geometry_restraints.dihedral(sites=[O_2, pseudoC_2, pseudoC_3, O_3],
      angle_ideal=-40, weight=1)

    res2.results['CO_d_in']  = d_in.angle_model
    res2.results['CO_d_out'] = d_out.angle_model
#-------------------------------------------------------------------------------
#}}}

#{{{ Ramachandran calculator
#Adds 'phi' 'psi' 'phi+1' and 'psi-1' to residue.results={} for each residue in
#  protein where protein is a dictionary of cablam_res classes
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

    try:
      CO_1     = res1.atomxyz['']['C']
      N_2      = res2.atomxyz['']['N']
      CA_2,C_2 = res2.atomxyz['']['CA'],res2.atomxyz['']['C']
      N_3      = res3.atomxyz['']['N']
    except KeyError:
      continue

    phi = geometry_restraints.dihedral(sites=[CO_1, N_2, CA_2, C_2],
      angle_ideal=-40, weight=1)
    psi = geometry_restraints.dihedral(sites=[N_2, CA_2, C_2, N_3],
      angle_ideal=-40, weight=1)

    res2.results['phi'] = phi.angle_model
    res2.results['psi'] = psi.angle_model

    res1.results['phi+1'] = phi.angle_model
    res3.results['psi-1'] = psi.angle_model
#-------------------------------------------------------------------------------
#}}}

#{{{ tau angle calculator
#Adds 'tau' to residue.results={} for each residue in protein where protein is
#  a dictionary of cablam_res classes
#-------------------------------------------------------------------------------
def taucalc(protein):
  for resid2 in protein:
    res2 = protein[resid2]

    try:
      N  = res2.atomxyz['']['N']
      CA = res2.atomxyz['']['CA']
      C  = res2.atomxyz['']['C']
    except KeyError:
      continue

    tau = geometry_restraints.angle(sites=[N,CA,C], angle_ideal=120, weight=1)
    res2.results['tau']  = tau.angle_model
#-------------------------------------------------------------------------------
#}}}
