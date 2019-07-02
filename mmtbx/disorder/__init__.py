
"""
Tools for analyzing disorder (or heterogeneity in general) in ensembles and
and other multi-conformer models or collections of conformations.  For building
and refinement of such models, see mmtbx.building.alternate_conformations or
mmtbx.refinement.ensemble_refinement.
"""

from __future__ import absolute_import, division, print_function
import math
import six

def set_ensemble_b_factors_to_xyz_displacement(pdb_hierarchy,
    include_hydrogens=False,
    include_waters=False,
    use_c_alpha_values=False,
    method="rmsf",
    selection=None,
    substitute_b_value=-1.0,
    logarithmic=False,
    log=None):
  """
  Given an ensemble (multi-MODEL PDB hierarchy), calculate the deviation
  between copies of each atom (defined here as either the root-mean-square
  fluctuation, or the radius of the minimum covering sphere) and set the
  isotropic B-factors to this value.
  """
  if (log is None) : log = null_out()
  assert (method in ["rmsf", "mcs"])
  from scitbx.math import minimum_covering_sphere
  from scitbx.array_family import flex
  pdb_atoms = pdb_hierarchy.atoms()
  pdb_atoms.reset_i_seq()
  xyz_by_atom = {}
  def get_key(atom):
    labels = atom.fetch_labels()
    return (labels.chain_id, labels.resid(), labels.altloc, atom.name)
  def get_c_alpha(atom):
    if (atom.name.strip() == "CA") and (atom.element.strip() == "C"):
      return atom
    for other in atom.parent().atoms():
      if (other.name.strip() == "CA") and (other.element.strip() == "C"):
        return other
    return None
  for model in pdb_hierarchy.models():
    for atom in model.atoms():
      if (selection is not None):
        if (not selection[atom.i_seq]) : continue
      elif (not include_hydrogens) and (atom.element.strip() in ["H","D"]):
        continue
      elif (not include_waters) and (atom.parent().resname in ["HOH"]):
        continue
      if (use_c_alpha_values) and (atom.name.strip() != "CA"):
        continue
      atom_key = get_key(atom)
      if (atom_key in xyz_by_atom):
        xyz_by_atom[atom_key].append(atom.xyz)
      else :
        xyz_by_atom[atom_key] = flex.vec3_double([atom.xyz])
  dev_by_atom = {}
  for atom_key, xyz in six.iteritems(xyz_by_atom):
    if (method == "mcs"):
      mcs = minimum_covering_sphere(points=xyz, epsilon=0.1)
      radius = mcs.radius()
      if (logarithmic):
        radius = math.log(radius + 1.0)
      dev_by_atom[atom_key] = radius
    else :
      mean_array = flex.vec3_double(xyz.size(), xyz.mean())
      rmsf = xyz.rms_difference(mean_array)
      dev_by_atom[atom_key] = rmsf
  # NOTE it seems flex.double handles list generators, not sure about funky python3 dict.values() though
  all_dev = flex.double(list(dev_by_atom.values()))
  if (method == "mcs"):
    print("Distribution of sphere radii:", file=log)
  else :
    print("Distribution of root-mean-square fluctuation values:", file=log)
  flex.histogram(all_dev, n_slots=20).show(f=log, prefix="  ",
    format_cutoffs="%.2f")
  for model in pdb_hierarchy.models():
    for atom in model.atoms():
      if (use_c_alpha_values):
        c_alpha = get_c_alpha(atom)
        if (c_alpha is None):
          atom.b = substitute_b_value
        else :
          atom_key = get_key(c_alpha)
          atom.b = dev_by_atom.get(atom_key, substitute_b_value)
      else :
        atom_key = get_key(atom)
        atom.b = dev_by_atom.get(atom_key, substitute_b_value)
