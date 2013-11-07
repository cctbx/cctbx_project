# -*- coding: utf-8; py-indent-offset: 2 -*-
"""
A script to open a model and its data and dump out all information about its
ions sites to a pickle file.
"""
from __future__ import division, print_function

import os
from cPickle import dump
import sys

from ion_utils import iterate_sites

from cctbx.eltbx import sasaki
from iotbx.file_reader import any_file
from libtbx.str_utils import make_header
from libtbx.utils import Usage
from libtbx.phil import parse
from mmtbx import ions
from mmtbx.ions.environment import ChemicalEnvironment, ScatteringEnvironment
from mmtbx.ions.parameters import server
from mmtbx.utils import cmdline_load_pdb_and_data

master_phil = parse("""
include scope mmtbx.utils.cmdline_input_phil_str
include scope mmtbx.ions.ion_master_phil
debug = True
  .type = bool
elements = Auto
  .type = str
wavelength = None
  .type = float
nproc = Auto
  .type = int
""", process_includes = True)

METAL_SERVER = server()

def _main(args, out = sys.stdout):
  if len(args) == 0 or "--help" in args:
    raise Usage("""\
mmtbx.dump_sites model.pdb data.mtz [options ...]

Utility to dump information about the properties of water and ion sites in a
model. This properties include local environment, electron density maps, and
atomic properties.

Full parameters:
{}
""".format(master_phil.as_str(prefix=" ")))

  cmdline = cmdline_load_pdb_and_data(
    args = args,
    master_phil = master_phil,
    out = out,
    process_pdb_file = True,
    create_fmodel = True,
    prefer_anomalous = True)
  params = cmdline.params

  if params.wavelength is None:
    pdb_in = any_file(params.input.pdb.file_name[0],
      force_type="pdb")
    wavelength = pdb_in.file_object.extract_wavelength()
    if wavelength is not None:
      print("", file = out)
      print("Using wavelength = {} from PDB header".format(wavelength),
            file = out)
      params.wavelength = wavelength

  if params.wavelength is not None:
    cmdline.xray_structure.set_inelastic_form_factors(
      photon = params.wavelength,
      table="sasaki")
    cmdline.fmodel.update_xray_structure(
      cmdline.xray_structure,
      update_f_calc = True)

  make_header("Inspecting sites", out = out)

  manager = ions.create_manager(
    pdb_hierarchy = cmdline.pdb_hierarchy,
    fmodel = cmdline.fmodel,
    geometry_restraints_manager = cmdline.geometry,
    wavelength = params.wavelength,
    params = params,
    verbose = params.debug,
    nproc = params.nproc,
    log = out)

  manager.show_current_scattering_statistics(out = out)

  sites = dump_sites(manager)

  out_name = os.path.splitext(params.input.pdb.file_name[0])[0] + "_sites.pkl"
  print("Dumping to", out_name)

  with open(out_name, "w") as f:
    dump(sites, f)

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

  properties = \
    [(
      ChemicalEnvironment(
        atom.i_seq,
        manager.find_nearby_atoms(atom.i_seq, far_distance_cutoff = 3.5),
        manager),
      ScatteringEnvironment(
        atom.i_seq,
        manager,
        fo_map)
      )
     for atom in atoms]

  del fo_map

  return properties

if __name__ == "__main__":
  _main(sys.argv[1:])
