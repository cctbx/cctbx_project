from __future__ import absolute_import, division, print_function
import math

import iotbx.pdb
from scitbx import matrix
import six
from six.moves import range

def construct_xyz(ba, bv,
                  aa, av,
                  da, dv,
                  period=3,
                  ):
  """Construct a list of cartesian coordinates based on Z-matrix information

  Args:
      ba (hierarchy atom): 1st atom
      bv (float): bond length to ba
      aa (hierarchy atom): 2nd atom
      av (float): angle value to ba-aa
      da (hierarchy atom): 3rd atom
      dv (float): torsion value to da-aa-ba
      period (int, optional): periodicity of torsion

  Returns:
      matrix: List of cartesian coordinates of lenght period


           bv
      ?--------ba
                 \
                  \
                av \
                    \
                     \
                     aa---------da
  """
  assert ba is not None
  assert aa is not None
  assert da is not None
  rn = matrix.col(ba.xyz)
  rca = matrix.col(aa.xyz)
  rc = matrix.col(da.xyz)
  rcca = rc -rca

  e0 = (rn - rca).normalize()
  e1 = (rcca - (rcca.dot(e0))*e0).normalize()
  e2 = e0.cross(e1)

  pi = math.pi
  alpha = math.radians(av)
  phi = math.radians(dv)

  rh_list = []
  for n in range(0, period):
    rh = rn + bv * (math.sin(alpha)*(math.cos(phi + n*2*pi/period)*e1 +
                                     math.sin(phi + n*2*pi/period)*e2) -
                    math.cos(alpha)*e0)
    rh_list.append(rh)
  return rh_list

def get_hierarchy_atom(name, element, xyz, occ=1., b=20.):
  atom = iotbx.pdb.hierarchy.atom()
  atom.name = name
  atom.element = element #"H"
  atom.xyz = xyz # rh3[i]
  atom.occ = occ # n.occ
  atom.b = b #n.b
  atom.segid = ' '*4
  return atom

def get_hierarchy_h_atom(name, xyz, heavy_atom, proton_element='H'):
  return get_hierarchy_atom(name,
                            proton_element,
                            xyz,
                            heavy_atom.occ,
                            heavy_atom.b,
                            )

def is_perdeuterated(ag):
  protons = {}
  for atom in ag.atoms():
    if atom.element_is_hydrogen():
      protons.setdefault(atom.element, 0)
      protons[atom.element]+=1
  if len(protons) in [0,2]: return False
  if len(protons)==1:
    if 'D' in protons:
      return True
    else:
      return False
  assert 0

def get_proton_info(ag):
  proton_name=proton_element='H'
  if is_perdeuterated(ag):
    proton_name=proton_element='D'
  return proton_element, proton_name

def generate_atom_group_atom_names(rg, names, return_Nones=False, verbose=True):
  '''
  Generate all alt. loc. groups of names
  '''
  atom_groups = rg.atom_groups()
  atom_altlocs = {}
  for ag in atom_groups:
    for atom in ag.atoms():
      atom_altlocs.setdefault(atom.parent().altloc, [])
      atom_altlocs[atom.parent().altloc].append(atom)
  if len(atom_altlocs)>1 and '' in atom_altlocs:
    for key in atom_altlocs:
      if key=='': continue
      for atom in atom_altlocs['']:
        atom_altlocs[key].append(atom)
    del atom_altlocs['']
  for key, value in six.iteritems(atom_altlocs):
    atoms=[]
    for name in names:
      for atom in value:
        if atom.name.strip()==name.strip():
          atoms.append(atom)
          break
      else:
        if return_Nones:
          atoms.append(None)
        else:
          if verbose:
            print('not all atoms found. missing %s from %s' % (name, names))
          break
    if len(atoms)!=len(names):
      yield None, None
    else:
      yield atoms[0].parent(), atoms

if __name__ == '__main__':
  ba = get_hierarchy_atom('BA', 'B', (0,0,0))
  aa = get_hierarchy_atom('AA', 'C', (0,1,0))
  da = get_hierarchy_atom('DA', 'O', (1,1,0))
  print(ba.format_atom_record())
  print(aa.format_atom_record())
  print(da.format_atom_record())
  rc = construct_xyz(ba, 1, aa, 90, da, 0)
  print(rc)
