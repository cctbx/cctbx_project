 # -*- coding: utf-8; py-indent-offset: 2 -*-

from __future__ import division
from __future__ import print_function

import os
from pickle import loads, dumps
from types import MethodType

import libtbx
from mmtbx.command_line.water_screen import master_phil
import mmtbx.command_line
from mmtbx import ions
from mmtbx.ions import environment
from mmtbx.regression.make_fake_anomalous_data import generate_zinc_inputs
import mmtbx.utils

def exercise():
  wavelength = 1.025
  mtz_file, pdb_file = generate_zinc_inputs(anonymize = False)
  null_out = libtbx.utils.null_out()

  cmdline = mmtbx.command_line.load_model_and_data(
    args = [pdb_file, mtz_file, "wavelength={}".format(wavelength),
            "use_phaser=False"],
    master_phil = master_phil(),
    out = null_out,
    process_pdb_file = True,
    create_fmodel = True,
    prefer_anomalous = True
    )

  os.remove(pdb_file)
  os.remove(mtz_file)

  cmdline.xray_structure.set_inelastic_form_factors(
    photon = cmdline.params.input.wavelength,
    table = "sasaki"
    )

  cmdline.fmodel.update_xray_structure(
    cmdline.xray_structure,
    update_f_calc = True
    )

  manager = ions.identify.create_manager(
    pdb_hierarchy = cmdline.pdb_hierarchy,
    fmodel = cmdline.fmodel,
    geometry_restraints_manager = cmdline.geometry,
    wavelength = cmdline.params.input.wavelength,
    params = cmdline.params,
    nproc = cmdline.params.nproc,
    log = null_out
    )

  manager.validate_ions(
    out = null_out
    )

  fo_map = manager.get_map("mFo")
  fofc_map = manager.get_map("mFo-DFc")

  for atom_props in manager.atoms_to_props.values():
    chem_env = environment.ChemicalEnvironment(
      atom_props.i_seq,
      manager.find_nearby_atoms(atom_props.i_seq, far_distance_cutoff = 3.5),
      manager
      )
    new_chem_env = loads(dumps(chem_env))
    for attr in dir(chem_env):
      if attr == "atom":
        # The two won't be directly comparable, but we will trust atom_labels is
        # tested fully in its own module
        assert chem_env.atom.id_str() == new_chem_env.atom.id_str()
      elif not attr.startswith("_") and \
        not isinstance(getattr(chem_env, attr), MethodType):
        assert getattr(chem_env, attr) == getattr(new_chem_env, attr)

    scatter_env = environment.ScatteringEnvironment(
      atom_props.i_seq, manager, fo_map, fofc_map
      )
    new_scatter_env = loads(dumps(scatter_env))
    for attr in dir(scatter_env):
      if not attr.startswith("_") and \
        not isinstance(getattr(scatter_env, attr), MethodType):
        assert getattr(scatter_env, attr) == getattr(new_scatter_env, attr)

  del fo_map
  del fofc_map

  print("OK")

if __name__ == "__main__":
  exercise()
