# -*- coding: utf-8; py-indent-offset: 2 -*-
"""
A script to open a model and its data and dump out all information about its
ions sites to a pickle file.

Run this module with:

phenix.python -m mmtbx.ions.svm.dump_sites [args]
"""

from __future__ import division
from mmtbx import ions
from mmtbx.ions.environment import ChemicalEnvironment, ScatteringEnvironment
from mmtbx.ions.svm.utils import iterate_sites
import mmtbx.command_line
from iotbx.pdb import common_residue_names_water as WATER_RES_NAMES
from cctbx.eltbx import sasaki
from libtbx.str_utils import make_header
from libtbx import easy_pickle
import os
import sys

def master_phil () :
  return mmtbx.command_line.generate_master_phil_with_inputs(
    enable_automatic_twin_detection=True,
    phil_string="""
include scope mmtbx.ions.ion_master_phil
debug = True
  .type = bool
elements = Auto
  .type = str
nproc = Auto
  .type = int
""")

def main(args, out = sys.stdout):
  usage_string = """\
mmtbx.dump_sites model.pdb data.mtz [options ...]

Utility to dump information about the properties of water and ion sites in a
model. This properties include local environment, electron density maps, and
atomic properties.
"""
  cmdline = mmtbx.command_line.load_model_and_data(
    args = args,
    master_phil = master_phil(),
    out = out,
    process_pdb_file = True,
    create_fmodel = True,
    prefer_anomalous = True,
    set_wavelength_from_model_header=True,
    set_inelastic_form_factors="sasaki",
    usage_string=usage_string)
  params = cmdline.params

  make_header("Inspecting sites", out = out)

  manager = ions.create_manager(
    pdb_hierarchy = cmdline.pdb_hierarchy,
    fmodel = cmdline.fmodel,
    geometry_restraints_manager = cmdline.geometry,
    wavelength = params.input.wavelength,
    params = params,
    verbose = params.debug,
    nproc = params.nproc,
    log = out)

  manager.show_current_scattering_statistics(out = out)

  sites = dump_sites(manager)

  out_name = os.path.splitext(params.input.pdb.file_name[0])[0] + "_sites.pkl"
  print >> out, "Dumping to", out_name
  easy_pickle.dump(out_name, sites)

def dump_sites (manager):
  """
  Iterate over all the ions and waters built into the model and dump out
  information about their properties.
  """

  atoms = iterate_sites(
    manager.pdb_hierarchy,
    res_filter = ions.SUPPORTED + ions.WATER_RES_NAMES,
    split_sites = True)

  # Can't pickle entire AtomProperties because they include references to the
  # Atom object. Instead, gather what properties we want and store them in a
  # second list
  fo_map = manager.get_map("mFo")
  fofc_map = manager.get_map("mFo-DFc")
  anom_map = manager.get_map("anom")

  properties = \
    [(
      ChemicalEnvironment(
        atom.i_seq,
        manager.find_nearby_atoms(atom.i_seq, far_distance_cutoff = 3.5),
        manager),
      ScatteringEnvironment(
        atom.i_seq,
        manager,
        fo_map,
        fofc_map,
        anom_map),
      )
     for atom in atoms]

  del fo_map
  del fofc_map
  del anom_map

  return properties

if __name__ == "__main__":
  main(sys.argv[1:])
