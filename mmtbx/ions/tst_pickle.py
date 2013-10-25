 # -*- coding: utf-8; py-indent-offset: 2 -*-

from __future__ import division

import os
from pickle import loads, dumps

import libtbx
from mmtbx.command_line.water_screen import master_phil
from mmtbx import ions
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

  manager.validate_ions(
    out = null_out
    )

  fo_map = manager.get_map("mFo")

  for atom_props in manager.atoms_to_props.values():
    new_atom_props = loads(dumps(atom_props))

    # Test all of the values used by the machine learning module
    assert atom_props.atom.id_str() == \
      new_atom_props.atom.id_str()
    assert atom_props.atom.segid.strip().upper() == \
      new_atom_props.atom.segid.strip().upper()
    assert atom_props.resname == new_atom_props.resname
    assert atom_props.b_iso == new_atom_props.b_iso
    assert atom_props.atom.occ == new_atom_props.atom.occ
    assert atom_props.d_min == new_atom_props.d_min
    assert atom_props.wavelength == new_atom_props.wavelength
    assert len(atom_props.nearby_atoms) == len(new_atom_props.nearby_atoms)
    assert atom_props.residue_counts == new_atom_props.residue_counts
    assert atom_props.peak_2fofc == new_atom_props.peak_2fofc
    assert atom_props.peak_fofc == new_atom_props.peak_fofc
    assert atom_props.peak_anom == new_atom_props.peak_anom
    assert atom_props.estimated_weight == new_atom_props.estimated_weight
    assert atom_props.fp == new_atom_props.fp
    assert atom_props.fpp == new_atom_props.fpp

    geometry = ions.find_coordination_geometry(
      atom_props.nearby_atoms, minimizer_method = True)
    new_geometry = ions.find_coordination_geometry(
      new_atom_props.nearby_atoms, minimizer_method = True)

    assert geometry == new_geometry
    assert ions.fit_gaussian(manager, atom_props.atom.xyz, fo_map) == \
      ions.fit_gaussian(manager, new_atom_props.atom.xyz, fo_map)

  print "OK"

if __name__ == "__main__":
  exercise()
