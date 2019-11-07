from __future__ import absolute_import, division, print_function
from scitbx.matrix import rotate_point_around_axis
from scitbx.math import dihedral_angle
from six.moves import cStringIO as StringIO

# flippable_sidechains:
#   key: flippable residue name
#   value: list of 5 atom names: first 4 define one angle,
#          1st, 2nd, 3rd and 5th define another torsion angle
#          6th, 7th - additional atoms to flip
flippable_sidechains = {
    "LEU": [" CA ", " CB ", " CG ", " CD1", " CD2"],
    "GLU": [" CB ", " CG ", " CD ", " OE1", " OE2"],
    "ASP": [" CA ", " CB ", " CG ", " OD1", " OD2"],
    "ASN": [" CA ", " CB ", " CG ", " OD1", " ND2"],
    "GLN": [" CB ", " CG ", " CD ", " OE1", " NE2"],
    "ARG": [" CD ", " NE ", " CZ ", " NH1", " NH2"],
    "VAL": [" C  ", " CA ", " CB ", " CG1", " CG2"],
    "PHE": [" CA ", " CB ", " CG ", " CD1", " CD2", " CE1", " CE2"],
    "HIS": [" CA ", " CB ", " CG ", " ND1", " CD2", " CE1", " NE2"],
    "TYR": [" CA ", " CB ", " CG ", " CD1", " CD2", " CE1", " CE2"],
}

def flip_residue(residue, mon_lib_srv=None):
  if residue.resname in flippable_sidechains:
    from mmtbx.utils import rotatable_bonds
    if mon_lib_srv is None:
      mon_lib_srv = mmtbx.monomer_library.server.server()
    null_log = StringIO()
    tardy_model = rotatable_bonds.tardy_model_one_residue(
        residue=residue,
        mon_lib_srv=mon_lib_srv,
        log=null_log)
    if tardy_model is None:
      return None
    clusters = tardy_model.tardy_tree.cluster_manager.clusters[1:]
    axes = tardy_model.tardy_tree.cluster_manager.hinge_edges[1:]
    assert len(axes) == len(clusters)
    residue_atoms = residue.atoms()
    ctr = 0
    ic = 0
    axis = None
    while ctr < len(axes):
      if residue_atoms[axes[ctr][0]].name.strip().upper() == flippable_sidechains[residue.resname][1].strip().upper() and \
          residue_atoms[axes[ctr][1]].name.strip().upper() == flippable_sidechains[residue.resname][2].strip().upper():
        axis = axes[ctr]
        ic = ctr
      ctr += 1
    atoms_to_rotate = []
    for ci in clusters[ic:]:
      atoms_to_rotate.extend(ci)
    if (axis is None):
      return None

    for atom in atoms_to_rotate:
      new_xyz = rotate_point_around_axis(
        axis_point_1=residue.atoms()[axis[0]].xyz,
        axis_point_2=residue.atoms()[axis[1]].xyz,
        point=residue.atoms()[atom].xyz,
        angle=180.0, deg=True)
      residue.atoms()[atom].xyz = new_xyz

def should_be_flipped(residue_1, residue_2):
  """ are these residues in similar flip orientation?"""
  # assert residue_1.resname == residue_2.resname, "%s %s" % (
  #     residue_1.id_str(), residue_2.id_str())
  if residue_1.resname != residue_2.resname:
    # e.g. 1u54, 3srv: PTR mutation of TYR residues
    return False
  if residue_1.resname in flippable_sidechains:
    tor_1_sites = []
    for aname in flippable_sidechains[residue_1.resname][:4]:
      a = residue_1.find_atom_by(name=aname)
      if a is None:
        return False
      else:
        tor_1_sites.append(a.xyz)
    tor_23_sites = []
    for aname in flippable_sidechains[residue_1.resname][:5]:
      a = residue_2.find_atom_by(name=aname)
      if a is None:
        return False
      else:
        tor_23_sites.append(a.xyz)
    tor1 = dihedral_angle(
      sites=tor_1_sites,
      deg=True)
    tor2 = dihedral_angle(
      sites=tor_23_sites[:4],
      deg=True)
    tor3 = dihedral_angle(
      sites=tor_23_sites[:3]+[tor_23_sites[4]],
      deg=True)
    if tor1 is None or tor2 is None or tor3 is None:
      return False
    return abs(tor1-tor2) > abs(tor1-tor3)+5
  return False
