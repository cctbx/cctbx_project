from __future__ import division
import math
from math import pi, sin, cos, atan2, sqrt

from mmtbx.conformation_dependent_library.cdl_utils import round_to_ten
from mmtbx.validation import mean_devs_PRO_phi_psi
from mmtbx.validation import mean_devs_VAL_THR_ILE_phi_psi
from mmtbx.validation import mean_devs_others_phi_psi
from mmtbx.validation import mean_devs_VAL_THR_ILE_pre_PRO_phi_psi
from mmtbx.validation import mean_devs_others_pre_PRO_phi_psi

# (phi, psi) -> (r in Angstrom, t in radians) of the mean C-beta deviation
tables = {
  'PRO'                 : mean_devs_PRO_phi_psi.mean_devs,
  'VAL_THR_ILE'         : mean_devs_VAL_THR_ILE_phi_psi.mean_devs,
  'others'              : mean_devs_others_phi_psi.mean_devs,
  'VAL_THR_ILE_pre_PRO' : mean_devs_VAL_THR_ILE_pre_PRO_phi_psi.mean_devs,
  'others_pre_PRO'      : mean_devs_others_pre_PRO_phi_psi.mean_devs,
  }

class radial_deviation:
  def __init__(self, r, t, input_in_radians=False, input_in_cartesian=False):
    if input_in_cartesian:
      # r, t are x, y
      r, t = sqrt(r*r + t*t), atan2(t, r)
      input_in_radians=True
    self.r = r
    self.t = t
    if not input_in_radians: self.t *= pi/180

  def __repr__(self):
    #return u'%0.2f \u2220 %6.1f' % (self.r, self.t/pi*180)
    return "'%4.1f/_%6.1f'" % (self.r*100, self.t/pi*180)

  def __add__(self, other):
    r = sqrt(self.r**2 + other.r**2 + 2*self.r*other.r*cos(other.t-self.t))
    numer = other.r * sin(other.t-self.t)
    denom = self.r + other.r * cos(other.t-self.t)
    t = self.t + atan2(numer, denom)
    t = atan2(sin(t), cos(t)) # wrap to (-pi, pi] to match dihedral_angle
    new = radial_deviation(0,0)
    new.r = r
    new.t = t
    return new

  def __sub__(self, other):
    new = radial_deviation(0,0)
    new.r = other.r
    new.t = other.t + pi
    return self.__add__(new)

  def __truediv__(self, other):
    if type(other)==type(1):
      self.r/=other
      return self
    else:
      assert 0

  def __cmp__(self, other):
    if self.r<other.r: return -1
    return 1

  def __eq__(self, other):
    return self.r == other.r

  def __ne__(self, other):
    return self.r != other.r

  def __lt__(self, other):
    return self.r < other.r

  def __le__(self, other):
    return self.r <= other.r

  def __gt__ (self, other):
    return self.r > other.r

  def __ge__(self, other):
    return self.r >= other.r

def get_table_resname(resname, use_parent_for_modified_residues=False):
  """Residue name used to choose the phi/psi table: the name itself for
  standard amino acids and, if use_parent_for_modified_residues, the parent
  amino acid from iotbx.pdb.modified_aa_names for other residues. None if
  the residue isn't corrected, which includes the D-amino acids."""
  from iotbx.pdb import common_residue_names_get_class
  if common_residue_names_get_class(resname) == 'common_amino_acid':
    return resname
  if not use_parent_for_modified_residues: return None
  from iotbx.pdb.aa_utils import get_aa_parent
  return get_aa_parent(resname)

def get_correction_table(table_resname, next_resname):
  """The table for a residue: PRO whatever follows, otherwise VAL/THR/ILE or
  the others, each split by whether the next residue is PRO."""
  if table_resname == 'PRO': return tables['PRO']
  if table_resname in ['VAL', 'THR', 'ILE']: name = 'VAL_THR_ILE'
  else: name = 'others'
  if next_resname == 'PRO': name += '_pre_PRO'
  return tables[name]

def get_correction_vector(resname,
                          next_resname,
                          phi_psi,
                          use_parent_for_modified_residues=False):
  """Mean deviation vector (x, y) in Angstrom for the residue's (phi, psi)
  bin, in the (deviation, dihedral_NABB) polar plane, or None if the residue
  isn't corrected or its bin has no entry."""
  resname = get_table_resname(resname,
    use_parent_for_modified_residues=use_parent_for_modified_residues)
  if resname is None: return None
  table = get_correction_table(resname, next_resname)
  entry = table.get((round_to_ten(phi_psi[0]), round_to_ten(phi_psi[1])), None)
  if entry is None: return None
  return (entry[0]*cos(entry[1]), entry[0]*sin(entry[1]))

# the fitted strength is pulled toward prior_strength as if prior_weight more
# residues agreed with it
prior_strength = 0.5
prior_weight = 100

def apply_correction(deviation, dihedral, correction, strength):
  """(deviation, dihedral in degrees) of v - strength*c, where v is the
  deviation vector and c the correction (x, y) in Angstrom."""
  t = math.radians(dihedral)
  x = deviation*cos(t) - strength*correction[0]
  y = deviation*sin(t) - strength*correction[1]
  return sqrt(x*x + y*y), math.degrees(atan2(y, x))

def choose_strength(deviations, dihedrals, corrections, user_strength=None):
  """Correction strength for one model from its residues that have a
  correction. fitted is the least-squares strength sum(v.c)/sum(|c|**2),
  pulled the same with the prior; applied is pulled limited to 0..1, or 0 if
  that would raise the mean deviation, or user_strength if given."""
  from libtbx import group_args
  n = len(corrections)
  svc = scc = 0.
  for d, t, c in zip(deviations, dihedrals, corrections):
    t = math.radians(t)
    svc += d*cos(t)*c[0] + d*sin(t)*c[1]
    scc += c[0]*c[0] + c[1]*c[1]
  fitted = pulled = None
  if n and scc:
    fitted = svc/scc
    m = scc/n
    pulled = (svc + prior_weight*prior_strength*m)/(scc + prior_weight*m)
  rc = group_args(n=n,
                  fitted=fitted,
                  pulled=pulled,
                  applied=0.,
                  user_set=(user_strength is not None),
                  check_failed=False,
                  mean_uncorrected=None,
                  mean_corrected=None)
  if user_strength is not None:
    rc.applied = user_strength
  elif pulled is not None:
    rc.applied = min(1., max(0., pulled))
  if not n: return rc
  rc.mean_uncorrected = sum(deviations)/n
  rc.mean_corrected = sum(apply_correction(d, t, c, rc.applied)[0]
    for d, t, c in zip(deviations, dihedrals, corrections))/n
  if not rc.user_set and rc.mean_corrected > rc.mean_uncorrected:
    rc.check_failed = True
    rc.applied = 0.
  return rc

def get_phi_psi_correction(result,
                           residue,
                           phi_psi,
                           next_resname=None,
                           display_phi_psi_correction=False,
                           use_parent_for_modified_residues=False,
                           verbose=False):
  correction = get_correction_vector(residue.resname, next_resname, phi_psi,
    use_parent_for_modified_residues=use_parent_for_modified_residues)
  if correction:
    correction = radial_deviation(*correction, input_in_cartesian=True)
    current = radial_deviation(result.deviation, result.dihedral)
    rc = current-correction
    start = (current.r>=0.25)
    finish = (rc.r>=0.25)
    if verbose:
      print(f'{current=} corrected with {correction=} ~> {rc=}')
    if display_phi_psi_correction and (start or finish):
      show_phi_psi_correction(residue, phi_psi, current, correction, rc)
    dihedral = (rc.t/pi*180 + 180) % 360 - 180
    return rc.r, dihedral, start, finish
  else:
    return None

def show_phi_psi_correction(residue, phi_psi, current, correction, rc, units=100):
  import matplotlib.pyplot as plt

  xs = [current.t, correction.t, rc.t]
  ys = [current.r, correction.r, rc.r]
  labels = ['default', 'correction', 'result']

  for i, (x, y) in enumerate(zip(xs, ys)):
      y*=units
      plt.polar(x, y, 'ro')
      plt.text(x, y,
               '%s' % (labels[i]),
               )
  key = (round_to_ten(phi_psi[0]), round_to_ten(phi_psi[1]))
  plt.title('%s (%5.1f, %5.1f) (%4.0f, %4.0f)\n%s minus %s equals %s' %(
    residue.id_str(),
    phi_psi[0],
    phi_psi[1],
    key[0],
    key[1],
    current,
    correction,
    rc,
   )
  )
  plt.show()

if __name__ == '__main__':
  # tests
  v1 = radial_deviation(1,0)
  v2 = radial_deviation(-1,0)
  print(v1,v2)
  print(v1-v2)
  v1 = radial_deviation(1,1)
  print(v1)
  v1 = radial_deviation(1,1, input_in_radians=True)
  print(v1)
  v1 = radial_deviation(1,1, input_in_cartesian=True)
  print(v1)
