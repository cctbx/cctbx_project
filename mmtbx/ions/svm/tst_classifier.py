 # -*- coding: utf-8; py-indent-offset: 2 -*-
from __future__ import division

import os
import sys

import libtbx
from mmtbx.command_line.water_screen import master_phil
from mmtbx.ions.environment import ChemicalEnvironment, ScatteringEnvironment
from mmtbx import ions
from mmtbx.ions.svm import ion_class, predict_ion
from mmtbx.regression.make_fake_anomalous_data import generate_zinc_inputs, \
     generate_calcium_inputs
import mmtbx.command_line

def exercise () :
  fns = [generate_calcium_inputs, generate_zinc_inputs]
  wavelengths = [1.025, 1.54]

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
      prefer_anomalous = True,
      set_inelastic_form_factors = "sasaki",
      )

    os.remove(pdb_file)
    os.remove(mtz_file)
    os.remove(os.path.splitext(mtz_file)[0] + "_fmodel.eff")
    os.remove(os.path.splitext(mtz_file)[0] + ".pdb")

    manager = ions.create_manager(
      pdb_hierarchy = cmdline.pdb_hierarchy,
      fmodel = cmdline.fmodel,
      geometry_restraints_manager = cmdline.geometry,
      wavelength = cmdline.params.wavelength,
      params = cmdline.params,
      nproc = cmdline.params.nproc,
      log = null_out
      )

    # Build a list of properties of each water / ion site
    waters = []
    for chain in manager.pdb_hierarchy.only_model().chains():
      for residue_group in chain.residue_groups():
        atom_groups = residue_group.atom_groups()
        if (len(atom_groups) > 1) : # alt conf, skip
          continue
        for atom_group in atom_groups :
          # Check for non standard atoms in the residue
          # Or a label indicating the residue is a water
          resname = atom_group.resname.strip().upper()

          if (resname in ions.WATER_RES_NAMES) :
            atoms = atom_group.atoms()
            if (len(atoms) == 1) : # otherwise it probably has hydrogens, skip
              waters.append(atoms[0].i_seq)

    assert len(waters) > 0

    atom_props = [ions.AtomProperties(i_seq, manager) for i_seq in waters]

    fo_map = manager.get_map("mFo")
    fofc_map = manager.get_map("mFo-DFc")
    anom_map = manager.get_map("anom")

    for atom_prop in atom_props:
      chem_env = ChemicalEnvironment(
        atom_prop.i_seq,
        manager.find_nearby_atoms(atom_prop.i_seq, far_distance_cutoff = 3.5),
        manager,
        )
      scatter_env = ScatteringEnvironment(
        atom_prop.i_seq, manager, fo_map, fofc_map, anom_map
        )
      resname = ion_class(chem_env)
      assert resname != ""

      predictions = predict_ion(chem_env, scatter_env,
                                elements = ["HOH", "ZN", "CA"])
      if predictions is None:
        print "Could not load SVM classifier"
        print "Skipping {}".format(os.path.split(__file__)[1])
        return

      if resname != predictions[0][0]:
        print "Prediction ({}) did not match expected: {}" \
          .format(predictions[0][0], resname)
        for element, prob in predictions:
          print "  {}: {:.2f}".format(element, prob)
        sys.exit()

    del fo_map
    del fofc_map
    del anom_map

  print "OK"

if __name__ == "__main__":
  exercise()
