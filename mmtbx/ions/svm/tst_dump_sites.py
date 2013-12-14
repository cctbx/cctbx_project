 # -*- coding: utf-8; py-indent-offset: 2 -*-

from __future__ import division

import os
from pickle import load

# We must make sure to import sklearn before boost python
import sklearn.svm # import dependency

import libtbx
from mmtbx.regression.make_fake_anomalous_data import generate_zinc_inputs
from mmtbx.ions.svm.dump_sites import main

def exercise():
  wavelength = 1.025
  mtz_file, pdb_file = generate_zinc_inputs(anonymize = False)
  null_out = libtbx.utils.null_out()
  main(["skip_twin_detection=True", "use_phaser=False",
        "input.pdb.file_name=" + pdb_file,
        "input.xray_data.file_name=" + mtz_file,
        "wavelength={}".format(wavelength)], out = null_out)

  sites_path = os.path.splitext(pdb_file)[0] + "_sites.pkl"
  sites = load(open(sites_path))

  assert len(sites) == 7
  for chem_env, scatter_env in sites:
    assert chem_env is not None
    assert scatter_env is not None

  os.remove(pdb_file)
  os.remove(mtz_file)
  os.remove(sites_path)
  os.remove(os.path.splitext(pdb_file)[0] + "_fmodel.eff")

  print "OK"

if __name__ == "__main__":
  exercise()
