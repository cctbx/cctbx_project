
from __future__ import division
import libtbx.load_env
from cStringIO import StringIO
import os

def exercise () :
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/unsorted.pdb",
    test=os.path.isfile)
  if (pdb_file is None) :
    print "phenix_regression not found, skipping"
    return
  if (os.path.exists("unsorted_sorted.pdb")) :
    os.remove("unsorted_sorted.pdb")
  from mmtbx.command_line import sort_hetatms
  from iotbx import file_reader
  out = StringIO()
  sort_hetatms.run(
    args=[pdb_file, "--verbose"],
    out=out)
  assert ("""Residue group pdb=" O   HOH    10 " added to chain A (distance = 2.814, symop = -x-1,y-1/2,-z)""" in out.getvalue())
  pdb_in = file_reader.any_file("unsorted_sorted.pdb")
  hierarchy = pdb_in.file_object.construct_hierarchy()
  chains = hierarchy.models()[0].chains()
  assert (len(chains) == 4)
  rgsA = chains[-2].residue_groups()
  rgsB = chains[-1].residue_groups()
  assert (len(rgsA) == 3) and (len(rgsB) == 4)
  sort_hetatms.run(
    args=[pdb_file, "preserve_chain_id=True", "renumber=False"],
    out=StringIO())
  pdb_in = file_reader.any_file("unsorted_sorted.pdb")
  hierarchy = pdb_in.file_object.construct_hierarchy()
  chains = hierarchy.models()[0].chains()
  assert (len(chains) == 3)
  assert (chains[-1].id == "A")
  print "OK"

if (__name__ == "__main__") :
  exercise()
