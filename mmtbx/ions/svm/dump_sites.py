# -*- coding: utf-8; py-indent-offset: 2 -*-
"""
A script to open a model and its data and dump out all information about its
ions sites to a pickle file.

Run this module with:

phenix.python -m mmtbx.ions.svm.dump_sites [args]
"""

from __future__ import absolute_import, division, print_function

import os
import sys

from libtbx import easy_pickle
from libtbx.str_utils import make_header
from mmtbx import ions
from mmtbx.ions.environment import ChemicalEnvironment, ScatteringEnvironment
from mmtbx.ions.identify import WATER_RES_NAMES
from mmtbx.ions.svm.utils import iterate_sites
from mmtbx.command_line import load_model_and_data

def master_phil():
  from mmtbx.command_line import generate_master_phil_with_inputs
  return generate_master_phil_with_inputs(
    enable_automatic_twin_detection=True,
    enable_pdb_interpretation_params=True,
    enable_stop_for_unknowns=False,
    phil_string="""
include scope mmtbx.ions.identify.ion_master_phil
include scope mmtbx.ions.svm.svm_phil_str
debug = True
  .type = bool
elements = Auto
  .type = str
use_svm = False
  .type = bool
nproc = Auto
  .type = int
""")

def _main(args, out=sys.stdout):
  """
  Main entry point to this script.

  Parameters
  ----------
  args : list of str
      List of arguments, should not include the first argument with the
      executable name.
  out : file, optional
  """
  usage_string = """\
phenix.python -m mmtbx.ions.svm.dump_sites model.pdb data.mtz [options ...]

Utility to dump information about the properties of water and ion sites in a
model. This properties include local environment, electron density maps, and
atomic properties.
"""
  cmdline = load_model_and_data(
    args=args,
    master_phil=master_phil(),
    out=out,
    process_pdb_file=True,
    create_fmodel=True,
    prefer_anomalous=True,
    set_wavelength_from_model_header=True,
    set_inelastic_form_factors="sasaki",
    usage_string=usage_string,
    )

  params = cmdline.params
  params.use_svm = True

  make_header("Inspecting sites", out=out)

  manager = ions.identify.create_manager(
    pdb_hierarchy=cmdline.pdb_hierarchy,
    fmodel=cmdline.fmodel,
    geometry_restraints_manager=cmdline.geometry,
    wavelength=params.input.wavelength,
    params=params,
    verbose=params.debug,
    nproc=params.nproc,
    log=out,
    )

  manager.show_current_scattering_statistics(out=out)

  sites = dump_sites(manager)

  out_name = os.path.splitext(params.input.pdb.file_name[0])[0] + "_sites.pkl"
  print("Dumping to", out_name, file=out)
  easy_pickle.dump(out_name, sites)

def dump_sites(manager):
  """
  Iterate over all the ions and waters built into the model and dump out
  information about their properties.

  Parameters
  ----------
  manager : mmtbx.ions.identify.manager

  Returns
  -------
  list of tuple of mmtbx.ions.environment.ChemicalEnvironment, \
  mmtbx.ions.environment.ScatteringEnvironment
  """

  atoms = iterate_sites(
    manager.pdb_hierarchy,
    res_filter=ions.SUPPORTED + WATER_RES_NAMES,
    split_sites=True,
    )

  # Can't pickle entire AtomProperties because they include references to the
  # Atom object. Instead, gather what properties we want and store them in a
  # second list
  properties = []
  for atom in atoms:
    map_stats = manager.map_stats(atom.i_seq)
    fo_density = manager.get_map_gaussian_fit("mFo", atom.i_seq)
    chem_env = ChemicalEnvironment(
      atom.i_seq,
      manager.find_nearby_atoms(atom.i_seq, far_distance_cutoff=3.5),
      manager,
      )
    scatter_env = ScatteringEnvironment(
      atom.i_seq,
      manager,
      fo_density=fo_density,
      fofc_density=(map_stats.fofc, 0),
      anom_density=(map_stats.anom, 0),
      )
    properties.append((chem_env, scatter_env))

  return properties

if __name__ == "__main__":
  _main(sys.argv[1:])
