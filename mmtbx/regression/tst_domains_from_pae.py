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

def tst_01(log = sys.stdout):

    args = group_args(
      group_args_type = 'parameters',
      pae_file = pae_file,
      library = 'networkx',
      pae_power = 1.0,
      pae_cutoff = 5.0,
      resolution = 1.0,
      select_range = False)

    selections = get_domain_selections_from_pae_matrix(pae_file = args.pae_file,
        pae_power = args.pae_power, pae_cutoff = args.pae_cutoff,
         resolution =  args.resolution,)
    print("Selections:")
    for s in selections:
       print(s)
    assert selections == [
   "(resseq 0:391)",
   "(resseq 392:401)",
   "(resseq 402:721)"]

if __name__ == "__main__":

  t0 = time.time()
  tst_01()
  print ("Time:", time.time()-t0)
  print ("OK")

