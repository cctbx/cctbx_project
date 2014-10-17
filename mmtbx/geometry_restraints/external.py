
"""
Accessory module for interfacing phenix.refine (or similar programs) with
various external third-party software such as Amber, AFITT, DivCon, or
Rosetta.
"""

from __future__ import division
import libtbx.load_env
import os

external_energy_params_str = ""

amber_installed = False
if libtbx.env.has_module("amber_adaptbx") :
  build_dir = libtbx.env.under_build("amber_adaptbx")
  if (build_dir is not None) and (os.path.isdir(build_dir)) :
    amber_installed = True

if (amber_installed) :
  external_energy_params_str += """
    amber
      .help = Parameters for using Amber in refinement.  EXPERIMENTAL
      .expert_level = 3
    {
      include scope amber_adaptbx.master_phil_str
    }
"""

qb_installed = os.environ.get("QBHOME", False)
if (qb_installed) :
  external_energy_params_str += """
    qbio
      .short_caption = QuantumBio refinement plugin
    {
      include scope qbpy.qb_params.qblib_master_params
    }
"""

#afitt
afitt_installed=False
if os.environ.get("OE_EXE", None) :
  if os.path.exists(os.environ["OE_EXE"]):
    afitt_installed = True

if (afitt_installed) :
  external_energy_params_str += """
    afitt
      .help = Parameters for using AFITT ligand gradients.
      .expert_level = 3
    {
      include scope mmtbx.geometry_restraints.afitt.master_phil_str
    }
"""

#
rosetta_installed = False
if libtbx.env.has_module("rosetta_adaptbx") :
  build_dir = libtbx.env.under_build("rosetta_adaptbx")
  if (build_dir is not None) and (os.path.isdir(build_dir)) :
    rosetta_installed = True

if rosetta_installed :
  external_energy_params_str += """
    rosetta
      .help = Rosetta force field extension
      .expert_level = 3
    {
      use_rosetta_energy = False
        .type = bool
      include scope rosetta_adaptbx.master_phil
    }
"""
