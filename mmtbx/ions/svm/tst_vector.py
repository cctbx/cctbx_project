 # -*- coding: utf-8; py-indent-offset: 2 -*-

from __future__ import division

import os, time

import libtbx
from mmtbx.command_line.water_screen import master_phil
from mmtbx import ions
from mmtbx.ions.environment import ChemicalEnvironment, ScatteringEnvironment
from mmtbx.ions.svm import ion_class, ion_vector
from mmtbx.regression.make_fake_anomalous_data import generate_zinc_inputs
import mmtbx.command_line

import libtbx.load_env

def exercise():
  wavelength = 1.025
  mtz_file, pdb_file = generate_zinc_inputs(anonymize = False)
  null_out = libtbx.utils.null_out()

  cmdline = mmtbx.command_line.load_model_and_data(
    args = [pdb_file, mtz_file, "wavelength={}".format(wavelength),
            "use_phaser=False", "use_svm=True"],
    master_phil = master_phil(),
    out = null_out,
    process_pdb_file = True,
    create_fmodel = True,
    prefer_anomalous = True
    )

  os.remove(pdb_file)
  os.remove(mtz_file)
  os.remove(os.path.splitext(pdb_file)[0] + "_fmodel.eff")

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

  for atom_props in manager.atoms_to_props.values():
    i_seq = atom_props.i_seq
    chem_env = ChemicalEnvironment(
      i_seq,
      manager.find_nearby_atoms(i_seq, far_distance_cutoff = 3.5),
      manager
      )
    scatter_env = ScatteringEnvironment(
      i_seq, manager,
      fo_density = manager.get_map_gaussian_fit("mFo", i_seq),
      fofc_density = manager.get_map_gaussian_fit("mFo-DFc", i_seq),
      anom_density = manager.get_map_gaussian_fit("anom", i_seq),
      )
    vector = ion_vector(chem_env, scatter_env)
    resname = ion_class(chem_env)
    assert vector is not None
    assert resname != ""

  print "OK"


if __name__ == "__main__":
  keep_going=True
  try:
    import numpy as np
    import svm
    import svmutil
  except ImportError:
    print "Required third-party dependencies are missing, skipping test."
    keep_going=False
  if (libtbx.env.find_in_repositories(relative_path="chem_data") is None):
    print "Skipping tst_vector exercise(): chem_data directory not available"
  else:
    t0 = time.time()
    if(keep_going):
      exercise()
    print "Time: %6.2f"%(time.time()-t0)
