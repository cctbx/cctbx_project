"""Example of how to shake (randomize) a model"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
from cctbx.array_family import flex
import mmtbx.model
import sys

def run(args):
  assert len(args) == 1
  # Read file into pdb_input class
  inp = iotbx.pdb.input(file_name=args[0])

  # create a model manager
  model = mmtbx.model.manager(
      model_input = inp)

  # get number of atoms in the input model
  n_atoms = model.get_number_of_atoms()

  # extract atom coordinates
  old_sites_cart = model.get_sites_cart()
  # generate random additions
  random_addition = flex.vec3_double(
    flex.random_double(size=n_atoms*3)-0.5)
  # actually add them to old coordinates
  new_xyz = old_sites_cart + random_addition

  # Update coordinates in model manager
  model.set_sites_cart(sites_cart=new_xyz)

  # get xray structure
  xrs = model.get_xray_structure()

  # reset B-factors (min=1, max=20)
  # generate array of new B-factors
  new_b = flex.random_double(size=n_atoms, factor=19) + 1
  # set them in xray structure
  xrs.set_b_iso(values=new_b)
  # update model manager with this xray structure
  model.set_xray_structure(xrs)
  # output result in PDB format to the screen
  print(model.as_pdb_or_mmcif_string())
  print("END")

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
