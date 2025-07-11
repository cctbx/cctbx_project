"""Distance between two PDB files: per atom, residue, chain, model and overall"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.model_model_distances

from scitbx.array_family import flex
import sys
import iotbx.pdb
from libtbx.utils import Sorry
import mmtbx.utils
from six.moves import zip

legend = """phenix.model_model_distances:
  Given two PDB files output distances per atom, residue, chain, model and overall.
  It is assumed (and asserted in the code!) that the amount and order of atoms
  in both files are identical.

How to run:
  phenix.model_model_distances model_1.pdb model_2.pdb

Feedback:
  PAfonine@lbl.gov
  phenixbb@phenix-online.org"""

def compute_overall(h1,h2,log):
  x1 = h1.extract_xray_structure()
  x2 = h2.extract_xray_structure()
  d = x1.distances(other=x2)
  d= d.min_max_mean().as_tuple()
  print("Overall  min/max/mean: %8.3f %8.3f %8.3f"%d, file=log)
  print(file=log)

def compute_overall_backbone(h1,h2,log):
  selection_cache1 = h1.atom_selection_cache()
  isel1 = selection_cache1.iselection("name ca or name n or name o or name c")
  selection_cache2 = h2.atom_selection_cache()
  isel2 = selection_cache2.iselection("name ca or name n or name o or name c")
  x1 = h1.select(isel1).extract_xray_structure()
  x2 = h2.select(isel2).extract_xray_structure()
  d = x1.distances(other=x2)
  d= d.min_max_mean().as_tuple()
  print("Backbone min/max/mean: %8.3f %8.3f %8.3f"%d, file=log)
  print(file=log)


def compute_per_model(h1,h2,log):
  if(len(h1.models())==1): return
  print("Per model (min/max/mean):", file=log)
  for m1,m2 in zip(h1.models(), h2.models()):
    r1 = m1.atoms().extract_xyz()
    r2 = m2.atoms().extract_xyz()
    d = flex.sqrt((r1 - r2).dot()).min_max_mean().as_tuple()
    print(m1.id, ": %-8.3f %-8.3f %-8.3f"%d, file=log)
  print(file=log)

def compute_per_chain(h1,h2,log):
  cs1 = list(h1.chains())
  cs2 = list(h2.chains())
  if(len(cs1)==1): return
  print("Per chain (min/max/mean):", file=log)
  for c1, c2 in zip(cs1, cs2):
    label = c1.id
    r1 = c1.atoms().extract_xyz()
    r2 = c2.atoms().extract_xyz()
    d = flex.sqrt((r1 - r2).dot()).min_max_mean().as_tuple()
    print(label, ": %-8.3f %-8.3f %-8.3f"%d, file=log)
  print(file=log)

def compute_per_residue(h1,h2,log):
  rgs1 = list(h1.residue_groups())
  rgs2 = list(h2.residue_groups())
  if(len(rgs1)==1): return
  print("Per residue (min/max/mean):", file=log)
  for rg1, rg2 in zip(rgs1, rgs2):
    label = "%10s"%"/".join([
      rg1.parent().id.strip(),
      rg1.resid().strip(),
      "_".join(list(rg1.unique_resnames()))])
    r1 = rg1.atoms().extract_xyz()
    r2 = rg2.atoms().extract_xyz()
    d = flex.sqrt((r1 - r2).dot()).min_max_mean().as_tuple()
    print(label, ": %-8.3f %-8.3f %-8.3f"%d, file=log)
  print(file=log)

def compute_per_atom(h1,h2,log):
  as1 = list(h1.atoms())
  as2 = list(h2.atoms())
  if(len(as1)==1): return
  print("Per atom:", file=log)
  for a1, a2 in zip(as1, as2):
    r1 = flex.vec3_double([a1.xyz])
    r2 = flex.vec3_double([a2.xyz])
    d = flex.sqrt((r1 - r2).dot())
    print(a1.format_atom_record()[:30], ": %-8.3f"%d[0], file=log)

def run(args, log=sys.stdout):
  print("-"*79, file=log)
  print(legend, file=log)
  print("-"*79, file=log)
  inputs = mmtbx.utils.process_command_line_args(args = args)
  file_names = inputs.pdb_file_names
  if(len(file_names) != 2): raise Sorry("Two PDB files has to given.")
  pi1 = iotbx.pdb.input(file_name = file_names[0])
  pi2 = iotbx.pdb.input(file_name = file_names[1])
  if pi1.crystal_symmetry_from_cryst1() is not None:
    if(not pi1.crystal_symmetry_from_cryst1().is_similar_symmetry(
           pi2.crystal_symmetry_from_cryst1())):
      raise Sorry("CRYST1 records must be identical.")
  h1 = pi1.construct_hierarchy()
  h2 = pi2.construct_hierarchy()
  if(not h1.is_similar_hierarchy(h2)):
    raise Sorry("Input PDB files have different content or atom order.")
  #
  print(file=log)
  compute_overall(    h1=h1, h2=h2, log=log)
  compute_overall_backbone(h1=h1, h2=h2, log=log)
  compute_per_model(  h1=h1, h2=h2, log=log)
  compute_per_chain(  h1=h1, h2=h2, log=log)
  compute_per_residue(h1=h1, h2=h2, log=log)
  compute_per_atom(   h1=h1, h2=h2, log=log)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])

