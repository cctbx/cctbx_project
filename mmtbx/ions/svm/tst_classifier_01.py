 # -*- coding: utf-8; py-indent-offset: 2 -*-

from __future__ import absolute_import, division, print_function
from mmtbx.ions.svm.dump_sites import master_phil
from mmtbx.ions.environment import ChemicalEnvironment, ScatteringEnvironment
from mmtbx import ions
from mmtbx.ions.identify import WATER_RES_NAMES, AtomProperties
from mmtbx.ions.svm import ion_class, predict_ion
from mmtbx.regression.make_fake_anomalous_data import generate_calcium_inputs
import mmtbx.command_line
import libtbx.load_env
import warnings
import os
import sys
import time

def exercise(prefix="ions_svm_classifier_01"):
  wavelength = 1.025
  mtz_file, pdb_file = generate_calcium_inputs(
      file_base=prefix,
      anonymize = True)
  null_out = libtbx.utils.null_out()

  cmdline = mmtbx.command_line.load_model_and_data(
    args = [pdb_file, mtz_file, "wavelength={}".format(wavelength),
            "use_phaser=True", "use_svm=True"],
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

  manager = ions.identify.create_manager(
    pdb_hierarchy = cmdline.pdb_hierarchy,
    fmodel = cmdline.fmodel,
    geometry_restraints_manager = cmdline.geometry,
    wavelength = cmdline.params.input.wavelength,
    params = cmdline.params,
    nproc = cmdline.params.nproc,
    log = null_out,
    manager_class = ions.svm.manager,
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

        if (resname in WATER_RES_NAMES):
          atoms = atom_group.atoms()
          if (len(atoms) == 1) : # otherwise it probably has hydrogens, skip
            waters.append(atoms[0].i_seq)

  assert len(waters) > 0

  atom_props = [AtomProperties(i_seq, manager) for i_seq in waters]

  for atom_prop in atom_props:
    i_seq = atom_prop.i_seq
    chem_env = ChemicalEnvironment(
      i_seq,
      manager.find_nearby_atoms(i_seq, far_distance_cutoff = 3.5),
      manager,
      )
    scatter_env = ScatteringEnvironment(
      i_seq, manager,
      fo_density = manager.get_map_gaussian_fit("mFo", i_seq),
      fofc_density = manager.get_map_gaussian_fit("mFo-DFc", i_seq),
      anom_density = manager.get_map_gaussian_fit("anom", i_seq),
      )
    resname = ion_class(chem_env)
    assert resname != ""

    predictions = predict_ion(chem_env, scatter_env,
                              elements = ["HOH", "ZN", "CA"])
    if predictions is None:
      print("Could not load SVM classifier")
      print("Skipping {}".format(os.path.split(__file__)[1]))
      return

    if resname != predictions[0][0]:
      print("Prediction ({}) did not match expected: {}" \
        .format(predictions[0][0], resname))
      for element, prob in predictions:
        print("  {}: {:.2f}".format(element, prob))
      sys.exit()

  print("OK")

if __name__ == "__main__":
  if (not libtbx.env.find_in_repositories("chem_data")):
    warnings.warn("chem_data not available, skipping this test")
  else :
    try :
      import svm
    except ImportError :
      warnings.warn("libsvm not available, skipping this test")
    else :
      t0 = time.time()
      exercise()
      print("Time: %6.2f"%(time.time()-t0))
      print("OK")
