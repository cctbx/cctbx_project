 # -*- coding: utf-8; py-indent-offset: 2 -*-

from __future__ import division

import os
from pickle import loads, dumps
import sys

import libtbx
from mmtbx.command_line.water_screen import master_phil
from mmtbx import ions
import mmtbx.monomer_library.pdb_interpretation
from mmtbx import monomer_library
from mmtbx.regression.make_fake_anomalous_data import generate_zinc_inputs
import mmtbx.utils

def exercise():
  wavelength = 0
  mtz_file, pdb_file = generate_zinc_inputs(anonymize = False)
  null_out = libtbx.utils.null_out()
  null_out = sys.stdout

  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
    args = [pdb_file, mtz_file, "wavelength={}".format(wavelength),
            "use_phaser=False"],
    master_phil = master_phil,
    out = null_out,
    process_pdb_file = True,
    create_fmodel = True,
    prefer_anomalous = True
    )

  cmdline.xray_structure.set_inelastic_form_factors(
    photon = cmdline.params.wavelength,
    table = "sasaki"
    )

  cmdline.fmodel.update_xray_structure(
    cmdline.xray_structure,
    update_f_calc = True
    )

  manager = ions.create_manager(
    pdb_hierarchy = cmdline.pdb_hierarchy,
    fmodel = cmdline.fmodel,
    geometry_restraints_manager = cmdline.geometry,
    wavelength = cmdline.params.wavelength,
    params = cmdline.params,
    nproc = cmdline.params.nproc,
    log = null_out
    )

  os.remove(pdb_file)
  os.remove(mtz_file)

  manager = loads(dumps(manager))

  print "OK"

if __name__ == "__main__":
  exercise()
