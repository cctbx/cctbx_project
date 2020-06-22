from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import os, sys
from libtbx.utils import Sorry
from mmtbx.model import manager as model_manager

def exercise(file_name, out = sys.stdout):
  if not os.path.isfile(file_name):
    raise Sorry("Missing the file: %s" %(file_name)+"\n"+"Please run with phenix.python read_write_mrc.py my_mrcfile.mrc")

  print ("Reading from %s" %(file_name))
  from iotbx.map_manager import map_manager

  m = map_manager(file_name)

  # make a little model
  sites_cart = flex.vec3_double( ((8, 10, 12), (14, 15, 16)))
  model = model_manager.from_sites_cart(
         atom_name = ' CA ',
         resname = 'ALA',
         chain_id = 'A',
         b_iso = 30.,
         occ = 1.,
         scatterer = 'C',
         sites_cart = sites_cart,
         crystal_symmetry = m.crystal_symmetry())

  # make a map_model_manager with lots of maps and model and ncs
  from iotbx.map_model_manager import map_model_manager

  from mmtbx.ncs.ncs import ncs
  ncs_object=ncs()
  ncs_object.set_unit_ncs()
  mam = map_model_manager(
          map_manager =  m,
          ncs_object =  ncs_object,
          map_manager_1 =  m.deep_copy(),
          map_manager_2 =  m.deep_copy(),
          map_manager_list =  [m.deep_copy(),m.deep_copy()],
          map_manager_id_list = ["extra_1","extra_2"],
          model     = model.deep_copy(),)

  r_model=mam.as_r_model()

  print ("OK")

if __name__ == "__main__":
  args = sys.argv[1:]
  if not args:
    import libtbx.load_env
    file_name = libtbx.env.under_dist(
      module_name = "iotbx",
      path = "ccp4_map/tst_input.map")
    args = [file_name]
  exercise(file_name = args[0])
