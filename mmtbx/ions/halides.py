
"""
Functions specific to identifying halide ions.  Very preliminary - needs
much more sophisticated analysis.
"""

from __future__ import absolute_import, division, print_function
from mmtbx.ions import server
from scitbx.matrix import col, distance_from_plane
from libtbx import Auto

chloride_params_str = """
max_distance_to_amide_n = 3.5
  .type = float
  .input_size = 80
  .help = Max distance to amide N atom
max_distance_to_cation = 3.5
  .type = float
  .input_size = 80
min_distance_to_anion = 3.5
  .type = float
  .input_size = 80
min_distance_to_other_sites = 2.1
  .type = float
  .input_size = 80
max_distance_to_hydroxyl = 3.5
  .type = float
delta_amide_h_angle = 20
  .type = float
  .input_size = 80
  .help = Allowed deviation (in degrees) from ideal N-H-X angle
delta_planar_angle = 10
  .type = float
  .input_size = 80
max_deviation_from_plane = 0.8
  .type = float
"""

def is_negatively_charged_oxygen(atom_name, resname):
  """
  Determine whether the oxygen atom of interest is either negatively charged
  (usually a carboxyl group or sulfate/phosphate), or has a lone pair (and
  no hydrogen atom) that would similarly repel anions.

  Parameters
  -----------
  atom_name : str
  resname : str

  Returns
  -------
  bool
  """
  if ((atom_name in ["OD1","OD2","OE1","OE2"]) and
      (resname in ["GLU","ASP","GLN","ASN"])):
    return True
  elif ((atom_name == "O") and (not resname in ["HOH","WAT"])):
    return True # sort of - the lone pair acts this way
  elif ((len(atom_name) == 3) and (atom_name[0:2] in ["O1","O2","O3"]) and
        (atom_name[2] in ["A","B","G"])):
    return True
  elif (resname in ["SO4","PO4"]):
    return True
  return False

# XXX distance cutoff may be too generous, but 0.5 is too strict
def _is_coplanar_with_sidechain(atom, residue, distance_cutoff=0.75):
  """
  Given an isolated atom and an interacting residue with one or more amine
  groups, determine whether the atom is approximately coplanar with the terminus
  of the sidechain (and thus interacting with the amine hydrogen(s) along
  approximately the same axis as the N-H bond).

  Parameters
  ----------
  atom : iotbx.pdb.hierarchy.atom
  residue : iotbx.pdb.hierarchy.residue
  distance_cutoff : float, optional

  Returns
  -------
  bool
  """
  sidechain_sites = []
  resname = residue.resname
  for other in residue.atoms():
    name = other.name.strip()
    if (resname == "ARG") and (name in ["NH1","NH2","NE"]):
      sidechain_sites.append(other.xyz)
    elif (resname == "GLN") and (name in ["OE1","NE2","CD"]):
      sidechain_sites.append(other.xyz)
    elif (resname == "ASN") and (name in ["OD1","ND2","CG"]):
      sidechain_sites.append(other.xyz)
  if (len(sidechain_sites) != 3) : # XXX probably shouldn't happen
    return False
  D = distance_from_plane(atom.xyz, sidechain_sites)
  #print atom.id_str(), D
  return (D <= distance_cutoff)

def is_favorable_halide_environment(
    i_seq,
    contacts,
    pdb_atoms,
    sites_frac,
    connectivity,
    unit_cell,
    params,
    assume_hydrogens_all_missing=Auto):
  """
  Detects if an atom's site exists in a favorable environment for a halide
  ion. This includes coordinating by a positively charged sidechain or backbone
  as well as an absense of negatively charged coordinating groups.

  Parameters
  ----------
  i_seq : int
  contacts : list of mmtbx.ions.environment.atom_contact
  pdb_atoms : iotbx.pdb.hierarchy.af_shared_atom
  sites_frac : tuple of float, float, float
  connectivity : scitbx.array_family.shared.stl_set_unsigned
  unit_cell : uctbx.unit_cell
  params : libtbx.phil.scope_extract
  assume_hydrogens_all_missing : bool, optional

  Returns
  -------
  bool
  """
  if (assume_hydrogens_all_missing in [None, Auto]):
    elements = pdb_atoms.extract_element()
    assume_hydrogens_all_missing = not ("H" in elements or "D" in elements)
  atom = pdb_atoms[i_seq]
  binds_amide_hydrogen = False
  near_cation = False
  near_lys = False
  near_hydroxyl = False
  xyz = col(atom.xyz)
  min_distance_to_cation = max(unit_cell.parameters()[0:3])
  min_distance_to_hydroxyl = min_distance_to_cation
  for contact in contacts :
    # to analyze local geometry, we use the target site mapped to be in the
    # same ASU as the interacting site
    def get_site(k_seq):
      return unit_cell.orthogonalize(
        site_frac = (contact.rt_mx * sites_frac[k_seq]))
    other = contact.atom
    resname = contact.resname()
    atom_name = contact.atom_name()
    element = contact.element
    distance = abs(contact)
    j_seq = other.i_seq
    # XXX need to figure out exactly what this should be - CL has a
    # fairly large radius though (1.67A according to ener_lib.cif)
    if (distance < params.min_distance_to_other_sites):
      return False
    if not element in ["C", "N", "H", "O", "S"]:
      charge = server.get_charge(element)
      if charge < 0 and distance <= params.min_distance_to_anion:
        # Nearby anion that is too close
        return False
      if charge > 0 and distance <= params.max_distance_to_cation:
        # Nearby cation
        near_cation = True
        if (distance < min_distance_to_cation):
          min_distance_to_cation = distance
    # Lysine sidechains (can't determine planarity)
    elif (atom_name in ["NZ"] and #, "NE", "NH1", "NH2"] and
          resname in ["LYS"] and
          distance <= params.max_distance_to_cation):
      near_lys = True
      if (distance < min_distance_to_cation):
        min_distance_to_cation = distance
    # sidechain amide groups, no hydrogens (except Arg)
    # XXX this would be more reliable if we also calculate the expected
    # hydrogen positions and use the vector method below
    elif (atom_name in ["NZ","NH1","NH2","ND2","NE2"] and
          resname in ["ARG","ASN","GLN"] and
          (assume_hydrogens_all_missing or resname == "ARG") and
          distance <= params.max_distance_to_cation):
      if (_is_coplanar_with_sidechain(atom, other.parent(),
            distance_cutoff = params.max_deviation_from_plane)):
        binds_amide_hydrogen = True
        if (resname == "ARG") and (distance < min_distance_to_cation):
          min_distance_to_cation = distance
    # hydroxyl groups - note that the orientation of the hydrogen is usually
    # arbitrary and we can't determine precise bonding
    elif ((atom_name in ["OG1", "OG2", "OH1"]) and
          (resname in ["SER", "THR", "TYR"]) and
          (distance <= params.max_distance_to_hydroxyl)):
      near_hydroxyl = True
      if (distance < min_distance_to_hydroxyl):
        min_distance_to_hydroxyl = distance
    # Backbone amide, explicit H
    elif atom_name in ["H"]:
      # TODO make this more general for any amide H?
      xyz_h = col(contact.site_cart)
      bonded_atoms = connectivity[j_seq]
      if (len(bonded_atoms) != 1):
        continue
      xyz_n = col(get_site(bonded_atoms[0]))
      vec_hn = xyz_h - xyz_n
      vec_hx = xyz_h - xyz
      angle = abs(vec_hn.angle(vec_hx, deg = True))
      # If Cl, H, and N line up, Cl binds the amide group
      if abs(angle - 180) <= params.delta_amide_h_angle:
        binds_amide_hydrogen = True
      else :
        pass #print "%s N-H-X angle: %s" % (atom.id_str(), angle)
    # Backbone amide, implicit H
    elif atom_name in ["N"] and assume_hydrogens_all_missing:
      xyz_n = col(contact.site_cart)
      bonded_atoms = connectivity[j_seq]
      ca_same = c_prev = None
      for k_seq in bonded_atoms :
        other2 = pdb_atoms[k_seq]
        if other2.name.strip().upper() in ["CA"]:
          ca_same = col(get_site(k_seq))
        elif other2.name.strip().upper() in ["C"]:
          c_prev = col(get_site(k_seq))
      if ca_same is not None and c_prev is not None:
        xyz_cca = (ca_same + c_prev) / 2
        vec_ncca = xyz_n - xyz_cca
        # 0.86 is the backbone N-H bond distance in geostd
        xyz_h = xyz_n + (vec_ncca.normalize() * 0.86)
        vec_nh = xyz_n - xyz_h
        vec_nx = xyz_n - xyz
        angle = abs(vec_nh.angle(vec_nx, deg = True))
        if abs(angle - 180) <= params.delta_amide_h_angle:
          binds_amide_hydrogen = True
    # sidechain NH2 groups, explicit H
    elif ((atom_name in ["HD1","HD2"] and resname in ["ASN"]) or
          (atom_name in ["HE1","HE2"] and resname in ["GLN"])):
          # XXX not doing this for Arg because it can't handle the bidentate
          # coordination
          #(atom_name in ["HH11","HH12","HH21","HH22"] and resname == "ARG")):
      bonded_atoms = connectivity[j_seq]
      assert (len(bonded_atoms) == 1)
      xyz_n = col(get_site(bonded_atoms[0]))
      xyz_h = col(contact.site_cart)
      vec_nh = xyz_n - xyz_h
      vec_xh = xyz - xyz_h
      angle = abs(vec_nh.angle(vec_xh, deg = True))
      if abs(angle - 180) <= params.delta_amide_h_angle:
        binds_amide_hydrogen = True
      else :
        pass #print "%s amide angle: %s" % (atom.id_str(), angle)
  # now check again for negatively charged sidechain (etc.) atoms (e.g.
  # carboxyl groups), but with some leeway if a cation is also nearby.
  # backbone carbonyl atoms are also excluded.
  for contact in contacts :
    if (contact.altloc() not in ["", "A"]):
      continue
    resname = contact.resname()
    atom_name = contact.atom_name()
    distance = abs(contact)
    if ((distance < 3.2) and
        (distance < (min_distance_to_cation + 0.2)) and
        is_negatively_charged_oxygen(atom_name, resname)):
      #print contact.id_str(), distance
      return False
  return (binds_amide_hydrogen or near_cation or near_lys)
