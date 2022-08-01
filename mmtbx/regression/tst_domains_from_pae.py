from __future__ import division, print_function
import sys, os, time

from libtbx import group_args
import networkx as nx
from packaging import version

import libtbx.load_env
data_dir = libtbx.env.under_dist(
  module_name="mmtbx",
  path="regression",
  test=os.path.isdir)

from mmtbx.domains_from_pae import get_domain_selections_from_pae_matrix

pae_file=os.path.join(data_dir,'pae.json')
model_file=os.path.join(data_dir, 'pdbs','pae_model.pdb')
pae_file_v3=os.path.join(data_dir,'AF_json_v3.json')
model_file_v3=os.path.join(data_dir, 'pdbs','AF_json_v3.pdb')

from iotbx.data_manager import DataManager
dm = DataManager()
distance_model = dm.get_model(model_file)
distance_model.add_crystal_symmetry_if_necessary()
distance_model_v3 = dm.get_model(model_file_v3)
distance_model_v3.add_crystal_symmetry_if_necessary()

def tst_01(log = sys.stdout):

    if version.parse(nx.__version__) < version.parse('2.6.2'):
      pae_power = 2.0
      pae_cutoff = 5.0
      resolution = 1.0
    else:
      pae_power = 1.0
      pae_cutoff = 5.0
      resolution = 0.5

    args = group_args(
      group_args_type = 'parameters',
      pae_file = pae_file,
      library = 'networkx',
      pae_power = pae_power,
      pae_cutoff = pae_cutoff,
      resolution = resolution,
      select_range = False)

    selections = get_domain_selections_from_pae_matrix(pae_file = args.pae_file,
        library = args.library,
        pae_power = args.pae_power, pae_cutoff = args.pae_cutoff,
         graph_resolution =  args.resolution,)
    if version.parse(nx.__version__) < version.parse('2.6.2'):
      assert selections == [
         "(resseq 0:113) or (resseq 184:187)",
         "(resseq 114:182) or (resseq 188:291)",
         "(resseq 183:183)",
         "(resseq 292:308)"
       ], selections
    else:
      assert selections == ['(resseq 0:117) or (resseq 181:181) or (resseq 183:187)', '(resseq 118:180) or (resseq 182:182) or (resseq 188:308)'], selections


def tst_02(log = sys.stdout):

    if version.parse(nx.__version__) < version.parse('2.6.2'):
      pae_power = 2.0
      pae_cutoff = 5.0
      resolution = 1.0
    else:
      pae_power = 1.0
      pae_cutoff = 5.0
      resolution = 0.5

    args = group_args(
      group_args_type = 'parameters',
      pae_file = pae_file,
      library = 'networkx',
      pae_power = pae_power,
      pae_cutoff = pae_cutoff,
      resolution = resolution,
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
    if version.parse(nx.__version__) < version.parse('2.6.2'):
        assert selections == [
          "(resseq 0:1) or (resseq 22:113) or (resseq 184:187)",
          "(resseq 2:21)",
          "(resseq 114:183) or (resseq 188:291)",
          "(resseq 292:308)"
       ]
    else:
        assert selections == ['(resseq 0:24)', '(resseq 25:62) or (resseq 66:96) or (resseq 100:111)', '(resseq 63:65)', '(resseq 97:99)', '(resseq 112:117)', '(resseq 118:135) or (resseq 142:181) or (resseq 190:199)', '(resseq 136:141)', '(resseq 182:189)', '(resseq 200:269) or (resseq 274:288)', '(resseq 270:273)', '(resseq 289:308)'], selections

def tst_03(log = sys.stdout):

    if version.parse(nx.__version__) < version.parse('2.6.2'):
      pae_power = 2.0
      pae_cutoff = 5.0
      resolution = 1.0
    else:
      pae_power = 1.0
      pae_cutoff = 5.0
      resolution = 0.5

    args = group_args(
      group_args_type = 'parameters',
      pae_file = pae_file_v3,
      library = 'networkx',
      pae_power = pae_power,
      pae_cutoff = pae_cutoff,
      resolution = resolution,
      weight_by_ca_ca_distance = 1.0,
      distance_power = 1.0,
      distance_model = distance_model_v3,
      select_range = False)



    selections = get_domain_selections_from_pae_matrix(pae_file = args.pae_file,
        library=args.library,
        pae_power = args.pae_power, pae_cutoff = args.pae_cutoff,
         graph_resolution =  args.resolution,
         weight_by_ca_ca_distance = args.weight_by_ca_ca_distance,
         distance_power = args.distance_power,
         distance_model = args.distance_model)
    if version.parse(nx.__version__) < version.parse('2.6.2'):
        assert selections == ['(resseq 0:8)', '(resseq 9:16)',
   '(resseq 17:85) or (resseq 95:226) or (resseq 229:253) or (resseq 258:330)',
    '(resseq 86:93)', '(resseq 94:94)', '(resseq 227:228)',
       '(resseq 254:257)', '(resseq 331:334)'], selections
    else:
        assert selections == ['(resseq 0:5)', '(resseq 6:14)', '(resseq 15:35) or (resseq 44:59) or (resseq 67:83) or (resseq 97:111) or (resseq 121:135) or (resseq 145:156) or (resseq 172:177)', '(resseq 36:43) or (resseq 60:66)', '(resseq 84:86)', '(resseq 87:92)', '(resseq 93:96)', '(resseq 112:120) or (resseq 136:144) or (resseq 157:171) or (resseq 179:180)', '(resseq 178:178) or (resseq 181:181) or (resseq 199:208)', '(resseq 182:192) or (resseq 209:214)', '(resseq 193:198) or (resseq 215:223) or (resseq 235:249) or (resseq 266:285) or (resseq 293:309) or (resseq 318:332)', '(resseq 224:226)', '(resseq 227:234)', '(resseq 250:254)', '(resseq 255:259)', '(resseq 260:265)', '(resseq 286:292) or (resseq 310:317)', '(resseq 333:334)'], selections

if __name__ == "__main__":

  t0 = time.time()
  tst_01()
  print ("Time 01: ", time.time()-t0)
  t1 = time.time()
  tst_02()
  print ("Time 02: ", time.time()-t1)
  tst_03()
  print ("Time 03: ", time.time()-t1)

  print ("OK")

