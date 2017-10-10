 # -*- coding: utf-8; py-indent-offset: 2 -*-

from __future__ import division
from __future__ import print_function

import os
import sys

from libtbx.easy_pickle import load
from libtbx.easy_run import fully_buffered
from mmtbx.regression.make_fake_anomalous_data import generate_zinc_inputs
from mmtbx.ions.svm import ion_class

def exercise():
  wavelength = 1.025
  mtz_file, pdb_file = generate_zinc_inputs(anonymize=True)
  args = [
    "input.pdb.file_name=" + pdb_file,
    "input.xray_data.file_name=" + mtz_file,
    "wavelength={}".format(wavelength),
    ]
  fully_buffered(
    "phenix.python -m mmtbx.ions.svm.dump_sites " + " ".join(args),
    ).raise_if_errors()

  os.remove(pdb_file)
  os.remove(os.path.splitext(pdb_file)[0][:-4] + ".pdb")
  os.remove(mtz_file)
  # "zn_frag_hoh.pdb" => "zn_frag_fmodel.eff"
  os.remove(os.path.splitext(pdb_file)[0][:-4] + "_fmodel.eff")

  sites_path = os.path.splitext(pdb_file)[0] + "_sites.pkl"
  sites = load(sites_path)

  os.remove(sites_path)

  assert len(sites) == 7
  for chem_env, scatter_env in sites:
    assert chem_env is not None
    assert scatter_env is not None
    for name in chem_env.__slots__:
      if getattr(chem_env, name) is None:
        print("Error: chem_env.{} is not set".format(name))
        sys.exit()
    for name in scatter_env.__slots__:
      # f' is not set by phaser
      if name in ["fp"]:
        continue
      # Only check f'' for heavy metals
      if name != "fpp" or ion_class(chem_env) != "HOH":
        if getattr(scatter_env, name) is None:
          print("Error: scatter_env.{} is not set".format(name))
          sys.exit()

  print("OK")

if __name__ == "__main__":
  exercise()
