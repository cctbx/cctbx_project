
from __future__ import division
from libtbx import easy_run
import os

def exercise () :
  if (os.path.isfile("ca_frag.mtz")) :
    os.remove("ca_frag.mtz")
  from mmtbx.regression import make_fake_anomalous_data
  from iotbx.file_reader import any_file
  mtz_file, pdb_file = make_fake_anomalous_data.run()
  pdb_in = any_file(pdb_file)
  hierarchy = pdb_in.file_object.construct_hierarchy()
  for chain in hierarchy.models()[0].chains() :
    for residue_group in chain.residue_groups() :
      for atom_group in residue_group.atom_groups() :
        if (atom_group.resname == "CA ") :
          atom_group.resname = "HOH"
          atom = atom_group.atoms()[0]
          atom.name = " O  "
          atom.element = " O"
          atom.charge = ""
          atom.segid = "CA"
          break
  f = open("ca_frag_hoh.pdb", "w")
  f.write(hierarchy.as_pdb_string(pdb_in.file_object.crystal_symmetry()))
  f.close()
  args = ["ca_frag_hoh.pdb", "ca_frag.mtz", "wavelength=1.12", "nproc=1"]
  result = easy_run.fully_buffered("mmtbx.water_screen %s" % " ".join(args)
    ).raise_if_errors()
  assert ("  Probable element: CA+2" in result.stdout_lines)
  print "OK"

if (__name__ == "__main__") :
  exercise()
