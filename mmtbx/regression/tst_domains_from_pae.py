from __future__ import division, print_function
import sys, os, time

from libtbx import group_args

import libtbx.load_env
data_dir = libtbx.env.under_dist(
  module_name="mmtbx",
  path="regression",
  test=os.path.isdir)

from mmtbx.domains_from_pae import get_domain_selections_from_pae_matrix

pae_file=os.path.join(data_dir,'pae.json')
model_file=os.path.join(data_dir, 'pdbs','pae_model.pdb')

from iotbx.data_manager import DataManager
dm = DataManager()
distance_model = dm.get_model(model_file)
distance_model.add_crystal_symmetry_if_necessary()

def tst_01(log = sys.stdout):

    args = group_args(
      group_args_type = 'parameters',
      pae_file = pae_file,
      library = 'networkx',
      pae_power = 2.0,
      pae_cutoff = 5.0,
      resolution = 1.0,
      select_range = False)

    selections = get_domain_selections_from_pae_matrix(pae_file = args.pae_file,
        library = args.library,
        pae_power = args.pae_power, pae_cutoff = args.pae_cutoff,
         graph_resolution =  args.resolution,)
    print("Selections:")
    for s in selections:
       print(s)
    assert selections == [
      "(resseq 0:113) or (resseq 184:187)",
      "(resseq 114:182) or (resseq 188:291)",
      "(resseq 183:183)",
      "(resseq 292:308)"
   ]

def tst_02(log = sys.stdout):

    args = group_args(
      group_args_type = 'parameters',
      pae_file = pae_file,
      library = 'networkx',
      pae_power = 2.0,
      pae_cutoff = 5.0,
      resolution = 1.0,
      weight_by_ca_ca_distance = 1.0,
      distance_power = 1.0,
      distance_model = distance_model,
      select_range = False)



    selections = get_domain_selections_from_pae_matrix(pae_file = args.pae_file,
        library=args.library,
        pae_power = args.pae_power, pae_cutoff = args.pae_cutoff,
         graph_resolution =  args.resolution,
         weight_by_ca_ca_distance = args.weight_by_ca_ca_distance,
         distance_power = args.distance_power,
         distance_model = args.distance_model)
    print("Selections:")
    for s in selections:
       print(s)
    assert selections == [
      "(resseq 0:1) or (resseq 22:113) or (resseq 184:187)",
      "(resseq 2:21)",
      "(resseq 114:183) or (resseq 188:291)",
      "(resseq 292:308)"
   ]


if __name__ == "__main__":

  t0 = time.time()
  tst_01()
  print ("Time 01: ", time.time()-t0)
  t1 = time.time()
  tst_02()
  print ("Time 02: ", time.time()-t1)
  print ("OK")

