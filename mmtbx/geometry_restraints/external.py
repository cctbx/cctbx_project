
"""
Accessory module for interfacing phenix.refine (or similar programs) with
various external third-party software such as AFITT, DivCon, Schrodinger.
"""

from __future__ import absolute_import, division, print_function
import libtbx.load_env
import os

external_energy_params_str = ""

# amber_installed = False
# if libtbx.env.has_module("amber_adaptbx"):
#   build_dir = libtbx.env.under_build("amber_adaptbx")
#   try: import sander
#   except ImportError as e: sander = False
#   if sander:
#   #if (build_dir is not None) and (os.path.isdir(build_dir)):
#     amber_installed = True

# if (amber_installed):
#   external_energy_params_str += """
#     amber
#       .help = Parameters for using Amber in refinement.
#       .expert_level = 3
#     {
#       include scope amber_adaptbx.master_phil_str
#     }
# """

qb_installed = os.environ.get("QBHOME", False)
if qb_installed: qb_installed = os.path.isdir(qb_installed)
if (qb_installed):
  external_energy_params_str += """
    qbio
      .short_caption = QuantumBio refinement plugin
    {
      include scope qbpy.qb_params.qblib_master_params
    }
"""

#afitt
afitt_installed=False
if os.environ.get("OE_EXE", None):
  if os.path.exists(os.environ["OE_EXE"]):
    afitt_installed = True

if (afitt_installed):
  external_energy_params_str += """
    afitt
      .help = Parameters for using AFITT ligand gradients.
      .expert_level = 3
    {
      include scope mmtbx.geometry_restraints.afitt.master_phil_str
    }
"""

# Schrodinger
def is_schrodinger_installed(env):
  """
  Check if Schrodinger is installed and interface requested. Schrodinger is
  installed if SCHRODINGER env variable is set to root directory and if
  PHENIX_SCHRODINGER env variable is set.

  Parameters
  ----------
  env: dict
     Environment dictionary

  Returns
  -------
  out : bool
    Boolean indicating if Schrodinger is installed and requested.
  """
  return (
     env.get('PHENIX_SCHRODINGER', False)
     and env.get("SCHRODINGER", False)
     and os.path.exists(env["SCHRODINGER"])
  )

schrodinger_installed = is_schrodinger_installed(os.environ)
if schrodinger_installed:
  from glob import glob
  paths = glob(os.path.join(
      os.environ["SCHRODINGER"], 'psp-v*', 'data', 'phenix'))
  if len(paths) == 1:
    import sys
    sys.path.append(paths[0])
  try:
    import phenix_schrodinger
    external_energy_params_str += """
    schrodinger
      .help = Parameters for using Schrodinger's force fields.
      .expert_level = 3
    {
      include scope phenix_schrodinger.master_phil_str
    }
"""
  except ImportError:
    schrodinger_installed = False

def is_orca_installed(env):
  return (env.get('PHENIX_ORCA', False) and os.path.exists(env['PHENIX_ORCA']))

def orca_action():
  outl = '''
    orca
      .help = Orca
    {
      include scope mmtbx.geometry_restraints.qm_manager.orca_master_phil_str
    }
  '''
  return outl

def is_any_quantum_package_installed(env):
  installed = []
  actions = []
  for key, (question, action) in {'orca' : (is_orca_installed, orca_action),
                                  }.items():
    if question(os.environ):
      rc = action()
      actions.append(rc)
      installed.append(key)
  if installed:
    outl = '''
  qi
    .help = QM
    .expert_level = 3
  {
    use_quantum_interface = False
      .type = bool
    selection = None
      .type = atom_selection
    charge = 0
      .type = int
    multiplicity = 1
      .type = int
    buffer = 0.
      .type = float
    refine_buffer_hydrogen_atoms = True
      .type = bool
'''
    for action in actions:
      outl += action
    outl += '}'
    return outl

any_quantum_package_installed = is_any_quantum_package_installed(os.environ)
external_energy_params_str += any_quantum_package_installed
