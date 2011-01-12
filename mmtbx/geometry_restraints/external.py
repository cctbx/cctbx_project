
import libtbx.phil
import libtbx.load_env
from libtbx.utils import Sorry
import sys

is_installed = False
if libtbx.env.has_module("rosetta_adaptbx") :
  is_installed = True

if (not is_installed) :
  external_energy_params = libtbx.phil.parse("")
else :
  external_energy_params = libtbx.phil.parse("""
external_energy {
  use_rosetta_energy = False
    .type = bool
    .style = hidden
  rosetta {
    include scope rosetta_adaptbx.scoring.master_phil
  }
}
""", process_includes=True)

def check_if_enabled () :
  rosetta_adaptbx = None
  is_enabled = False
  if libtbx.env.has_module("rosetta_adaptbx") :
    try :
      import rosetta_adaptbx
    except ImportError, e :
      print "Error attempting to import rosetta_adaptbx:"
      print e
      rosetta_adaptbx = None
    else :
      is_enabled = True
  if (not is_enabled) :
    raise Sorry("External energy functions are not available in this " +
      "build environment.  Please set use_external_energies=False.")
  return is_enabled

def get_rosetta_manager (pdb_hierarchy,
                         params,
                         log=sys.stdout) :
  import rosetta_adaptbx # import dependency
  from rosetta_adaptbx import scoring
  manager = scoring.manager(
    pdb_hierarchy=pdb_hierarchy,
    params=params.rosetta,
    log=log)
  return manager
