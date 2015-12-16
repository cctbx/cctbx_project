from __future__ import division

from mmtbx.conformation_dependent_library import generate_protein_threes
from scitbx.matrix import rotate_point_around_axis
from mmtbx.validation import ramalyze #, RAMALYZE_FAVORED
import math
from cStringIO import StringIO

def get_phi_psi_atoms(hierarchy):
  phi_psi_atoms = []
  for three in generate_protein_threes(
        hierarchy=hierarchy,
        geometry=None):
    phi_atoms, psi_atoms = three.get_phi_psi_atoms()
    rama_key = three.get_ramalyze_key()
    # print "rama_key", rama_key
    phi_psi_atoms.append(([phi_atoms, psi_atoms],rama_key))
  return phi_psi_atoms

def get_dihedral_angle(atoms):
  from scitbx.math import dihedral_angle
  return dihedral_angle(
      sites = [x.xyz for x in atoms],
      deg=True)

def rama_score_evaluate(resType, value):
  # copy-paste from cctbx_project/mmtbx/validation/ramalyze.py
  # added 0.002 -> 0.0021, 0.0005 -> 0.00051, etc

  return ramalyze.ramalyze.evalScore(resType, value)

def pair_info(phi_psi_pair):
  return phi_psi_pair[0][2].id_str()

def list_rama_outliers_h(hierarchy, r):
  phi_psi_atoms = get_phi_psi_atoms(hierarchy)
  outp = list_rama_outliers(phi_psi_atoms, r)
  return outp


def list_rama_outliers(phi_psi_atoms, r):
  result = ""
  for phi_psi_pair, rama_key in phi_psi_atoms:
    rama_score = get_rama_score(phi_psi_pair, r, rama_key)
    if rama_evaluate(phi_psi_pair, r, rama_key) == ramalyze.RAMALYZE_OUTLIER:
      result += "  !!! OUTLIER %s, score=%f\n" % (pair_info(phi_psi_pair), rama_score)
    # print_rama_stats(phi_psi_atoms, r)
  return result


def get_rama_score(phi_psi_pair, r, rama_key):
  phi_psi_angles = get_pair_angles(phi_psi_pair)
  rama_score = r.evaluate(ramalyze.res_types[rama_key], phi_psi_angles)
  return rama_score

def rama_evaluate(phi_psi_pair, r, rama_key):
  score = get_rama_score(phi_psi_pair, r, rama_key)
  return rama_score_evaluate(rama_key, score)

def get_pair_angles(phi_psi_pair):
  phi_psi_angles = [0,0]
  phi_psi_angles[0] = get_dihedral_angle(phi_psi_pair[0])
  phi_psi_angles[1] = get_dihedral_angle(phi_psi_pair[1])
  return phi_psi_angles

def print_rama_stats(phi_psi_atoms, r):
  result = StringIO()
  for phi_psi_pair, rama_key in phi_psi_atoms:
    for i, atoms in enumerate(phi_psi_pair):
      for a in atoms:
        print >> result, a.id_str()
    rama_score = get_rama_score(phi_psi_pair, r, rama_key)
    print >> result, "rama score:", get_pair_angles(phi_psi_pair), rama_score,
    print >> result, rama_score_evaluate(rama_key, rama_score), rama_key
    print >> result, "="*20
  print >> result, "*"*80
  r = result.getvalue()
  return r

def get_rmsd(fixed_points, moving_points):
  rmsd = 0
  for fp, mp in zip(fixed_points, moving_points):
    rmsd += fp.distance(mp)**2
  return math.sqrt(rmsd)

def rotate_atoms_around_bond(
    moving_h, atom_axis_point_1, atom_axis_point_2, angle, degrees=True):
  # changes moving_h
  # print "in rotation, iseqs:", atom_axis_point_1.i_seq, atom_axis_point_2.i_seq
  #
  # find xyz based on i_seqs
  rotate_xyz1 = None
  rotate_xyz2 = None
  for a in moving_h.atoms():
    if a.i_seq == atom_axis_point_1.i_seq:
      rotate_xyz1 = a.xyz
    elif a.i_seq == atom_axis_point_2.i_seq:
      rotate_xyz2 = a.xyz
  # rotate stuff
  for a in moving_h.atoms():
    if a.i_seq > atom_axis_point_1.i_seq:
      new_xyz = rotate_point_around_axis(
          axis_point_1=rotate_xyz1,
          axis_point_2=rotate_xyz2,
          point=a.xyz,
          angle=angle,
          deg=degrees)
      # print "actually setting coordinates:", a.i_seq, a.xyz, "->", new_xyz
      a.set_xyz(new_xyz)

def find_nearest_non_outlier_region(phi_psi_pair, r, rama_key):
  def spiral(N, M):
      x,y = 0,0
      dx, dy = 0, -1
      for dumb in xrange(N*M):
          if abs(x) == abs(y) and [dx,dy] != [1,0] or x>0 and y == 1-x:
              dx, dy = -dy, dx            # corner, change direction
          if abs(x)>N/2 or abs(y)>M/2:    # non-square
              dx, dy = -dy, dx            # change direction
              x, y = -y+dx, x+dy          # jump
          yield x, y
          x, y = x+dx, y+dy
  # ==
  phi_psi_angles = get_pair_angles(phi_psi_pair)
  for dx,dy in spiral(360, 360):
    angles = [phi_psi_angles[0]+dx, phi_psi_angles[1]+dy]
    rama_score = r.evaluate(rama_key, angles)
    if rama_score_evaluate(rama_key, rama_score) == ramalyze.RAMALYZE_FAVORED:
      return angles
