 # -*- coding: utf-8; py-indent-offset: 2 -*-
from __future__ import division

import os
import sys

import libtbx
from mmtbx.command_line.water_screen import master_phil
from mmtbx.ions.environment import ChemicalEnvironment, ScatteringEnvironment
from mmtbx import ions
from mmtbx.ions.svm import ion_class, predict_ion, get_classifier
from mmtbx.regression.make_fake_anomalous_data import generate_zinc_inputs, \
     generate_calcium_inputs
import mmtbx.command_line

def exercise () :
  if get_classifier() is None:
    print "Skipping {}".format(os.path.split(__file__)[1])
    return

  # XXX: Should CA's wavelength be 1.025, too?
  fns = [generate_calcium_inputs, generate_zinc_inputs]
  wavelengths = [1.025, 1.025]

  for fn, wavelength in zip(fns, wavelengths):
    mtz_file, pdb_file = fn(anonymize = True)
    null_out = libtbx.utils.null_out()

    cmdline = mmtbx.command_line.load_model_and_data(
      args = [pdb_file, mtz_file, "wavelength={}".format(wavelength),
              "use_phaser=True"],
      master_phil = master_phil(),
      out = null_out,
      process_pdb_file = True,
      create_fmodel = True,
      prefer_anomalous = True
      )

    os.remove(pdb_file)
    os.remove(mtz_file)
    os.remove(os.path.splitext(mtz_file)[0] + "_fmodel.eff")
    os.remove(os.path.splitext(mtz_file)[0] + ".pdb")

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
    fofc_map = manager.get_map("mFo-DFc")
    anom_map = manager.get_map("anom")

    for atom_props in manager.atoms_to_props.values():
      chem_env = ChemicalEnvironment(
        atom_props.i_seq,
        manager.find_nearby_atoms(atom_props.i_seq, far_distance_cutoff = 3.5),
        manager,
        )
      scatter_env = ScatteringEnvironment(
        atom_props.i_seq, manager, fo_map, fofc_map, anom_map
        )
      resname = ion_class(chem_env)
      prediction = predict_ion(chem_env, scatter_env,
                               elements = ["ZN", "MN", "CA"])
      assert resname != ""
      if resname != prediction[0][0]:
        print "Prediction did not match expected: {} != {}" \
          .format(resname, prediction[0][0])
        for element, prob in prediction:
          print "  {}: {}".format(element, prob)
        # sys.exit()

    del fo_map
    del fofc_map
    del anom_map

  print "OK"

if __name__ == "__main__":
  exercise()
