
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
