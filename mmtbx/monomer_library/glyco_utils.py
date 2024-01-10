from __future__ import absolute_import, division, print_function
from string import digits
from cctbx import geometry_restraints
from libtbx.utils import Sorry
from six.moves import range

beta_1_4 = """
data_link_BETA1-4
loop_
_chem_link_bond.link_id
_chem_link_bond.atom_1_comp_id
_chem_link_bond.atom_id_1
_chem_link_bond.atom_2_comp_id
_chem_link_bond.atom_id_2
_chem_link_bond.type
_chem_link_bond.value_dist
_chem_link_bond.value_dist_esd
 BETA1-4  1  O4  2  C1  single  1.439  0.020

loop_
_chem_link_angle.link_id
_chem_link_angle.atom_1_comp_id
_chem_link_angle.atom_id_1
_chem_link_angle.atom_2_comp_id
_chem_link_angle.atom_id_2
_chem_link_angle.atom_3_comp_id
_chem_link_angle.atom_id_3
_chem_link_angle.value_angle
_chem_link_angle.value_angle_esd
 BETA1-4  1  C4  1  O4  2  C1  108.700  3.000
 BETA1-4  1  O4  2  C1  2  O5  112.300  3.000
 BETA1-4  1  O4  2  C1  2  C2  109.470  3.000
 BETA1-4  1  O4  2  C1  2  H1  109.470  3.000

loop_
_chem_link_chir.link_id
_chem_link_chir.atom_centre_comp_id
_chem_link_chir.atom_id_centre
_chem_link_chir.atom_1_comp_id
_chem_link_chir.atom_id_1
_chem_link_chir.atom_2_comp_id
_chem_link_chir.atom_id_2
_chem_link_chir.atom_3_comp_id
_chem_link_chir.atom_id_3
_chem_link_chir.volume_sign
 BETA1-4  2  C1  1  O4  2  O5  2  C2  positiv

"""

from mmtbx.monomer_library import linking_utils
from mmtbx.monomer_library import glyco_chiral_values

# atoms
anomeric_carbon   = "2 C1"
ring_oxygen       = "2 O5"
link_oxygen       = "1 O4"
link_carbon       = "1 C4"
anomeric_hydrogen = "2 H1"
ring_carbon       = "2 C2"

def get_chiral_sign(code):
  return glyco_chiral_values.volumes.get(code, None)

def get_alpha_beta(code, fake=True): # the fake chiral is needed to apply the
                                     # correct link from monomer library
                                     # needed because FUC is alpha but has a
                                     # chiral volume of beta
  if not fake:
    return glyco_chiral_values.alpha_beta.get(code, None)
  cs = get_chiral_sign(code)
  if cs is None: return None
  elif cs < 0: return "ALPHA"
  else: return "BETA"

atom_types = ["anomeric_carbon",
              "ring_oxygen",
              "ring_carbon",
              "link_oxygen",
              "link_carbon",
              "anomeric_hydrogen",
              ]


class glyco_link_class:
  def __init__(self,
               anomeric_carbon,
               ring_oxygen=None,
               ring_carbon=None,
               link_oxygen=None,
               link_carbon=None,
               anomeric_hydrogen=None,
               link_phi_carbon=None,
               ):
    self.anomeric_carbon=anomeric_carbon
    self.ring_oxygen=ring_oxygen
    self.ring_carbon=ring_carbon
    self.link_oxygen=link_oxygen
    self.link_carbon=link_carbon
    self.anomeric_hydrogen=anomeric_hydrogen
    self.anomeric_carbon_linking=None
    self.link_phi_carbon=link_phi_carbon

  def __repr__(self):
    outl = "\nGlycosidic atoms\n"
    for attr in atom_types:
      try: outl += "  %-20s : %s" % (attr, getattr(self, attr).quote())
      except Exception: outl += "  %-20s : ???" % (attr)
      if attr=="anomeric_carbon":
        outl += " linking : %s" % self.anomeric_carbon_linking
      outl += "\n"
    return outl

  def is_correct(self, verbose=False):
    if (self.anomeric_carbon is None or
        self.link_oxygen is None or
        self.ring_oxygen is None or
        #ring_carbon is None or
        self.link_carbon is None
        ):
      if verbose:
        attrs = [
          'anomeric_carbon',
          'link_oxygen',
          'ring_oxygen',
          'link_carbon',
          ]
        for attr in attrs:
          atom = getattr(self, attr)
          if atom: atom = atom.quote()
          print('  %-15s : %s' % (attr, atom))
      return False
    return True

  def get_chiral_i_seqs(self, return_names=False):
    i_seqs = []
    for atom in [self.anomeric_carbon,
                 self.link_oxygen,
                 self.ring_oxygen,
                 self.ring_carbon,
                 ]:
      if atom is None: return None
      if return_names:
        i_seqs.append(getattr(atom, "name"))
      else:
        i_seqs.append(getattr(atom, "i_seq"))
    return i_seqs

  def get_isomer(self, verbose=False):
    isomer = get_alpha_beta(self.anomeric_carbon.parent().resname)
    if verbose: print('1 isomer',isomer)
    if isomer is None: isomer = "?"
    if verbose: print('2 isomer',isomer)
    if self.anomeric_carbon.name.strip()[-1] in digits:
      isomer += self.anomeric_carbon.name.strip()[-1]
      if verbose: print('3 isomer',isomer)
    else:
      isomer += " %s " % self.anomeric_carbon.name.strip()
      if verbose: print('4 isomer',isomer)
    if self.link_oxygen.name.strip()[-1] in digits:
      isomer += "-%s" % self.link_oxygen.name.strip()[-1]
      if verbose: print('5 isomer',isomer)
    else:
      isomer += "- %s " % self.link_oxygen.name.strip()
      if verbose: print('6 isomer',isomer)
    if verbose: print('get_isomer',isomer)
    return isomer

  def get_code(self):
    names = []
    resnames = []
    for attr in atom_types:
      atom = getattr(self, attr, None)
      if atom is None: continue
      names.append(atom.name)
      resnames.append(atom.parent().resname)
    return "%s_%s" % (resnames[0], resnames[3])
    #if names == [' C1 ', ' O5 ', ' C2 ', ' O4 ', ' C4 ']: return ""

  def as_cif(self):
    import iotbx.cif
    isomer = self.get_isomer()
    code = self.get_code()
    code = "%s-%s" % (isomer, code)
    code = isomer
    cif_block = iotbx.cif.model.block()
    loop = iotbx.cif.model.loop(header=(
      "_chem_link_bond.link_id",
      "_chem_link_bond.atom_1_comp_id",
      "_chem_link_bond.atom_id_1",
      "_chem_link_bond.atom_2_comp_id",
      "_chem_link_bond.atom_id_2",
      "_chem_link_bond.type",
      "_chem_link_bond.value_dist",
      "_chem_link_bond.value_dist_esd", # neutron????
      ))
    loop.add_row((code,
                  "1",
                  self.anomeric_carbon.name.strip(),
                  "2",
                  self.link_oxygen.name.strip(),
                  "single",
                  "1.439", # need to check this
                  "0.020",
                  ))
    cif_block.add_loop(loop)
    loop = iotbx.cif.model.loop(header=(
      "_chem_link_angle.link_id",
      "_chem_link_angle.atom_1_comp_id",
      "_chem_link_angle.atom_id_1",
      "_chem_link_angle.atom_2_comp_id",
      "_chem_link_angle.atom_id_2",
      "_chem_link_angle.atom_3_comp_id",
      "_chem_link_angle.atom_id_3",
      "_chem_link_angle.value_angle",
      "_chem_link_angle.value_angle_esd",
      ))
    for (id1, a1, id2, a2, id3, a3), angle, esd in [
      [['1', self.link_carbon, '1', self.link_oxygen,     '2', self.anomeric_carbon],   108.7,  3.],
      [['1', self.link_oxygen, '2', self.anomeric_carbon, '2', self.ring_oxygen],       112.3,  3.],
      [['1', self.link_oxygen, '2', self.anomeric_carbon, '2', self.ring_carbon],       109.47, 3.],
      [['1', self.link_oxygen, '2', self.anomeric_carbon, '2', self.anomeric_hydrogen], 109.47, 3.],
      ]:
      if a1 is None or a2 is None or a3 is None: continue
      loop.add_row((code,
                    id1,
                    a1.name.strip(),
                    id2,
                    a2.name.strip(),
                    id3,
                    a3.name.strip(),
                    "%0.1f" % angle,
                    "%0.1f" % esd,
                    ))
    cif_block.add_loop(loop)
    loop = iotbx.cif.model.loop(header=(
      "_chem_link_chir.link_id",
      "_chem_link_chir.atom_centre_comp_id",
      "_chem_link_chir.atom_id_centre",
      "_chem_link_chir.atom_1_comp_id",
      "_chem_link_chir.atom_id_1",
      "_chem_link_chir.atom_2_comp_id",
      "_chem_link_chir.atom_id_2",
      "_chem_link_chir.atom_3_comp_id",
      "_chem_link_chir.atom_id_3",
      "_chem_link_chir.volume_sign",
      ))
    value = get_chiral_sign(self.anomeric_carbon.parent().resname)
    if value is not None and value>0:
      value = "positiv"
    else:
      value = "negativ"
    names = self.get_chiral_i_seqs(return_names=True)
    loop.add_row((code,
                  "2",
                  names[0].strip(),
                  "1",
                  names[1].strip(),
                  "2",
                  names[2].strip(),
                  "2",
                  names[3].strip(),
                  value,
                  ))
    cif_block.add_loop(loop)
    return cif_block

def get_distance2(a1, a2):
  d2 = 0
  for i in range(3):
    d2 += (a1.xyz[i]-a2.xyz[i])**2
  return d2

def generate_atoms_from_atom_groups(atom_group1, atom_group2):
  for atom in atom_group1.atoms(): yield atom
  for atom in atom_group2.atoms(): yield atom

def get_anomeric_carbon(atom_group1, atom_group2, bonds, verbose=False):
  for i, atom in enumerate(generate_atoms_from_atom_groups(atom_group1,
                                                           atom_group2)
                                                           ):
    if atom.element.strip() not in ["C"]: continue
    oxygens = []
    residues = []
    for ba in bonds.get(atom.i_seq, []):
      if ba.element.strip() in ["O"]:
        oxygens.append(ba)
        residues.append(ba.parent())
    if len(oxygens)==2:
      if residues[0].id_str() != residues[1].id_str():
        return atom, True
      else:
        # Eli changed this but not sure why
        continue
        #return atom, False
##         raise Sorry("""
##         Trying to find the anomeric carbon but found a carbon
##         linked to two oxygens.
##           anomeric carbon %s
##           linked oxygens  %s
##                           %s
##         The anomeric carbons should link to another residue.
##         """ % (atom.quote(),
##                oxygens[0].quote(),
##                oxygens[1].quote(),
##                )
##                )
  return None

def get_any_linking_carbon(atom_group1, atom_group2, bonds, verbose=False):
  for i, atom in enumerate(generate_atoms_from_atom_groups(atom_group1,
                                                           atom_group2)
                                                           ):
    if atom.element.strip() not in ["C"]: continue
    oxygens = []
    residues = []
    linking = False
    for ba in bonds.get(atom.i_seq, []):
      if ba.element.strip() in ["O"]:
        if atom.parent().id_str() != ba.parent().id_str():
          linking=True
        oxygens.append(ba)
        residues.append(ba.parent())
    if len(oxygens)==2: # anomeric carbon not linkable
      outl = 'Anomeric carbon "%s" has two oxygens. Need to remove non-ring O.' % (
        atom.quote(),
        )
      raise Sorry(outl)
    if linking:
      return atom, True

  return None

def get_C1_carbon(atom_group1,
                  atom_group2,
                  distance_cutoff=2.,
                  ):
  distance_cutoff *= distance_cutoff
  c1s = []
  for i, atom in enumerate(generate_atoms_from_atom_groups(atom_group1,
                                                           atom_group2)
                                                           ):
    if atom.name.strip()[:2] =="C1": c1s.append(atom)
    #Fix for atoms with A or B at end of name
  if not c1s:
    assert 0
    return None
  for c1 in c1s:
    oxygens = []
    for i, atom in enumerate(generate_atoms_from_atom_groups(atom_group1,
                                                             atom_group2)
                                                             ):
      if atom.element.strip() in ["O"]:
        d2 = get_distance2(c1, atom)
        if d2<distance_cutoff: # need from outside
          oxygens.append(atom)
    if len(oxygens)==2:
      break
  else:
    outl = ""
    for atom in c1s:
      outl += "\n\t\t%s" % atom.quote()
    raise Sorry("""Trying to find the linking carbons but could not find
        a suitable candidate.
%s
        Check carbohydrate geometry.
                """ % outl
               )
  if oxygens[0].parent().id_str()!=oxygens[1].parent().id_str():
    return c1, True
  else:
    return c1, False
##     raise Sorry("""
##         Trying to find the anomeric carbon but found a carbon
##         linked to two oxygens.
##           anomeric carbon %s
##           linked oxygens  %s
##                           %s
##         """ % (c1.quote(),
##                oxygens[0].quote(),
##                oxygens[1].quote(),
##                )
##                )
  return None

def get_ring_oxygen(anomeric_carbon, bonds, element='O'):
  for ba in bonds.get(anomeric_carbon.i_seq, []):
    if ba.element.strip() not in [element]: continue
    # check in same atom group
    if ba.parent().id_str() == anomeric_carbon.parent().id_str():
      return ba
    # check in same residue group
    if ba.parent().parent().id_str() == anomeric_carbon.parent().parent().id_str():
      return ba

def get_ring_oxygen_substitute(anomeric_carbon, bonds):
  return get_ring_oxygen(anomeric_carbon, bonds, element='C')

def get_ring_carbon(anomeric_carbon, bonds):
  for ba in bonds.get(anomeric_carbon.i_seq, []):
    if ba.element.strip() not in ["C"]: continue
    if ba.parent().id_str() == anomeric_carbon.parent().id_str():
      return ba
    if ba.parent().parent().id_str() == anomeric_carbon.parent().parent().id_str():
      return ba

def get_anomeric_hydrogen(anomeric_carbon, bonds):
  for ba in bonds.get(anomeric_carbon.i_seq, []):
    if ba.element.strip() not in ["H"]: continue
    if ba.parent().id_str() == anomeric_carbon.parent().id_str():
      return ba
    if ba.parent().parent().id_str() == anomeric_carbon.parent().parent().id_str():
      return ba

def get_link_oxygen(anomeric_carbon, bonds, verbose=False):
  if verbose:
    print(anomeric_carbon.quote())
    print(bonds.get(anomeric_carbon.i_seq))
  for ba in bonds.get(anomeric_carbon.i_seq, []):
    if verbose: print(ba.quote())
    if ba.element.strip() not in ["O"]: continue
    if ba.parent().id_str() != anomeric_carbon.parent().id_str():
      return ba
    if ba.parent().parent().id_str() != anomeric_carbon.parent().parent().id_str():
      return ba

def get_link_oxygen_on_distance(anomeric_carbon, atom_group1, atom_group2):
  link_group = None
  for atom in atom_group2.atoms():
    if atom.quote()==anomeric_carbon.quote():
      link_group = atom_group1
      break
  if link_group is None:
    for atom in atom_group1.atoms():
      if atom.quote()==anomeric_carbon.quote():
        link_group = atom_group2
        break
  if link_group is None: assert 0
  for atom in link_group.atoms():
    if atom.element.strip()!="O": continue
    d2 = get_distance2(atom, anomeric_carbon)
    if d2<5.:
      return atom
  return None

def get_link_carbon(anomeric_carbon, link_oxygen, bonds):
  for ba in bonds.get(link_oxygen.i_seq, []):
    if ba.element.strip() not in ["C"]: continue
    if ba.i_seq==anomeric_carbon.i_seq: continue
    if ba.parent().id_str() != anomeric_carbon.parent().id_str():
      return ba

def get_link_carbon_on_distance(anomeric_carbon, atom_group1, atom_group2):
  link_group = None
  for atom in atom_group2.atoms():
    if atom.quote()==anomeric_carbon.quote():
      link_group = atom_group1
      break
  if link_group is None:
    for atom in atom_group1.atoms():
      if atom.quote()==anomeric_carbon.quote():
        link_group = atom_group2
        break
  if link_group is None: assert 0
  for atom in link_group.atoms():
    if atom.element.strip()!="O": continue
    d2 = get_distance2(atom, anomeric_carbon)
    if d2<5.:
      return atom
  return None

def get_link_phi_carbon(link_carbon, bonds):
  phi_carbon = None
  if not link_carbon: return None
  for ba in bonds.get(link_carbon.i_seq, []):
    if ba.element.strip() not in ["C"]: continue
    if ba.i_seq==link_carbon.i_seq: continue
    phi_carbon = ba
  return phi_carbon

def get_glyco_link_atoms(atom_group1,
                         atom_group2,
                         link_carbon_dist=2.0,
                         verbose=False,
                         ):
  # maybe should be restraints based?
  bonds = linking_utils.get_bonded_from_atom_groups(atom_group1,
                                                    atom_group2,
                                                    link_carbon_dist,
    )
  rc = get_anomeric_carbon(atom_group1,
                           atom_group2,
                           bonds,
                           verbose=verbose)
  if rc is None:
    rc = get_any_linking_carbon(atom_group1,
                                atom_group2,
                                bonds,
                                verbose=verbose)
  if rc is None:
    rc = get_C1_carbon(atom_group1,
                       atom_group2,
                       distance_cutoff=link_carbon_dist,
      )
  if rc is None: return None

  anomeric_carbon, linking_carbon = rc
  if anomeric_carbon is None:
    assert 0
    return None
  if verbose: print('anomeric_carbon',anomeric_carbon.quote())
  ring_oxygen = get_ring_oxygen(anomeric_carbon, bonds)
  if verbose:
    try: print('ring_oxygen',ring_oxygen.quote())
    except AttributeError: print('ring_oxygen',ring_oxygen)
  if ring_oxygen is None:
    ring_oxygen = get_ring_oxygen_substitute(anomeric_carbon, bonds)
  ring_carbon = get_ring_carbon(anomeric_carbon, bonds)
  if verbose: print('ring_carbon',ring_carbon.quote())
  anomeric_hydrogen = get_anomeric_hydrogen(anomeric_carbon, bonds)
  if verbose: print('anomeric_hydrogen',anomeric_hydrogen)
  link_oxygen = get_link_oxygen(anomeric_carbon, bonds, verbose=verbose)
  if link_oxygen is None:
    link_oxygen = get_link_oxygen_on_distance(anomeric_carbon,
                                              atom_group1,
                                              atom_group2)
  if link_oxygen is None:
    return None
  if verbose: print('link_oxygen',link_oxygen.quote())
  link_carbon = get_link_carbon(anomeric_carbon, link_oxygen, bonds)
  if link_carbon is None and link_carbon_dist:
    link_carbon = get_link_carbon_on_distance(anomeric_carbon,
                                              atom_group1,
                                              atom_group2,
      )
  if verbose:
    try: print('link_carbon',link_carbon.quote())
    except Exception: print()

  link_phi_carbon = get_link_phi_carbon(link_carbon, bonds)
  if verbose:
    try: print('link_phi_carbon',link_phi_carbon.quote(), link_phi_carbon.name, ba.element.strip())
    except Exception: print()

  gla = glyco_link_class(anomeric_carbon,
                         ring_oxygen,
                         ring_carbon,
                         link_oxygen,
                         link_carbon,
                         anomeric_hydrogen,
                         link_phi_carbon,
                         )
  gla.anomeric_carbon_linking = linking_carbon
  return gla

def apply_glyco_link_using_proxies_and_atoms(atom_group1,
                                             atom_group2,
                                             bond_params_table,
                                             bond_asu_table,
                                             geometry_proxy_registries,
                                             rt_mx_ji,
                                             link_carbon_dist=2.0,
                                             origin_id=None,
                                             ):
  origin_ids = geometry_restraints.linking_class.linking_class()
  def _add_bond(i_seqs,
                bond_params_table,
                bond_asu_table,
                value,
                esd,
                rt_mx_ji,
                origin_id,
                ):
    proxy = geometry_restraints.bond_simple_proxy(
      i_seqs=i_seqs,
      distance_ideal=value,
      weight=1/esd**2,
      origin_id=origin_id)
    bond_params_table.update(i_seq=i_seqs[0],
                             j_seq=i_seqs[1],
                             params=proxy)
    bond_asu_table.add_pair(
      i_seq=i_seqs[0],
      j_seq=i_seqs[1],
      rt_mx_ji=rt_mx_ji,
      )
  #
  def _add_angle(i_seqs, geometry_proxy_registries, value, esd, origin_id):
    proxy = geometry_restraints.angle_proxy(
      i_seqs=i_seqs,
      angle_ideal=value,
      weight=1/esd**2,
      origin_id=origin_id)
    geometry_proxy_registries.angle.add_if_not_duplicated(proxy=proxy)
  #
  def _add_chiral(i_seqs, geometry_proxy_registries, value, esd, origin_id, both_signs=False):
    proxy = geometry_restraints.chirality_proxy(
      i_seqs=i_seqs,
      volume_ideal=value,
      both_signs=both_signs,
      weight=1/esd**2,
      origin_id=origin_id,
      )
    geometry_proxy_registries.chirality.add_if_not_duplicated(proxy=proxy)

  def atom_group_output(atom_group):
    outl = ""
    for atom in atom_group.atoms():
      outl += "%s%s\n" % (' '*10, atom.quote())
    return outl

  ########
  from mmtbx.monomer_library import glyco_utils
  gla = glyco_utils.get_glyco_link_atoms(atom_group1,
                                         atom_group2,
                                         link_carbon_dist=link_carbon_dist,
    )
  # checks
  if gla and not gla.is_correct():
    gla = glyco_utils.get_glyco_link_atoms(atom_group2,
                                           atom_group1,
                                           link_carbon_dist=link_carbon_dist,
      )
  if gla and not gla.is_correct():
    raise Sorry("""Failed to get all the atoms needed for glycosidic bond
between

      group 1
%s
      group 2
%s
    """ % (atom_group_output(atom_group1), atom_group_output(atom_group2)
    )
    )
  if gla is None:
    raise Sorry("""Unspecified problem with carbohydrate groups. Could be that
the linking oxygen is on the linking residue instead of the docking residue.

    group 1
%s
    group 2
%s
    """ % (atom_group_output(atom_group1), atom_group_output(atom_group2)
           )
           )
  if not gla.anomeric_carbon_linking:
    raise Sorry("""The linking carbohydrate unit has the oxygen attached to the
anomeric carbon.

  Consider replacing oxygen %s
  with an oxygen linked to  %s in the same residue
    %s""" % (gla.link_oxygen.quote(),
             gla.link_carbon.quote(),
             gla)
             )
  origin_id = origin_ids.get('link_%s' % gla.get_isomer(), None)
  if not origin_id: origin_id=origin_ids['glycosidic custom']
  i_seqs = [gla.anomeric_carbon.i_seq, gla.link_oxygen.i_seq]
  bond_i_seqs = i_seqs
  # bonds
  _add_bond(i_seqs, bond_params_table, bond_asu_table, 1.439, 0.02, rt_mx_ji, origin_id)
  # angles
  for i_atoms, value, esd in [
      [[gla.link_carbon, gla.link_oxygen,     gla.anomeric_carbon],   108.7,  3.],
      [[gla.link_oxygen, gla.anomeric_carbon, gla.ring_oxygen],       112.3,  3.],
      [[gla.link_oxygen, gla.anomeric_carbon, gla.ring_carbon],       109.47, 3.],
      [[gla.link_oxygen, gla.anomeric_carbon, gla.anomeric_hydrogen], 109.47, 3.],
    ]:
    if None in i_atoms: continue
    i_seqs = [atom.i_seq for atom in i_atoms]
    _add_angle(i_seqs, geometry_proxy_registries, value, esd, origin_id)
  # chiral
  i_seqs = gla.get_chiral_i_seqs()
  if i_seqs is None:
    raise Sorry("""Unable to determine the linking chiral atoms for atom groups
    group 1
%s
    group 2
%s
    """ % (atom_group_output(atom_group1), atom_group_output(atom_group2)
    )
    )
  value = get_chiral_sign(gla.anomeric_carbon.parent().resname)
  if value:
    esd = 0.02
    _add_chiral(i_seqs, geometry_proxy_registries, value, esd, origin_id)

  return gla.get_isomer(), gla.as_cif(), bond_i_seqs

def is_standard_glyco_link(key, link):
  if key.find("ALPHA")>-1 or key.find("BETA")>-1:
    return True
  return False
