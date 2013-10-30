 # -*- coding: utf-8; py-indent-offset: 2 -*-

from __future__ import division

import os
from pickle import loads, dumps

import libtbx
from mmtbx.command_line.water_screen import master_phil
from mmtbx import ions
from mmtbx.ions import environment
from mmtbx.regression.make_fake_anomalous_data import generate_zinc_inputs
import mmtbx.utils

def exercise():
  wavelength = 1.025
  mtz_file, pdb_file = generate_zinc_inputs(anonymize = False)
  null_out = libtbx.utils.null_out()

  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
    args = [pdb_file, mtz_file, "wavelength={}".format(wavelength),
            "use_phaser=False"],
    master_phil = master_phil,
    out = null_out,
    process_pdb_file = True,
    create_fmodel = True,
    prefer_anomalous = True
    )

  os.remove(pdb_file)
  os.remove(mtz_file)

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

  manager.validate_ions(
    out = null_out
    )

  fo_map = manager.get_map("mFo")

  for atom_props in manager.atoms_to_props.values():
    env = environment.Environment(
      atom_props.i_seq, atom_props.nearby_atoms, manager, fo_map = fo_map)
    new_env = loads(dumps(env))
    for attr in dir(env):
      if attr == "atom":
        # The two won't be directly comparable, but we will trust atom_labels is
        # tested fully in its own module
        assert env.atom.id_str() == new_env.atom.id_str()
      elif not attr.startswith("_"):
        assert getattr(env, attr) == getattr(new_env, attr)

  del fo_map

  print "OK"

if __name__ == "__main__":
  exercise()
