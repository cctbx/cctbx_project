# -*- coding: utf-8; py-indent-offset: 2 -*-
from __future__ import division
import errno
from libtbx import slots_getstate_setstate
from libtbx.command_line.easy_qsub import run as run_easy_qsub
from libtbx.command_line.easy_qsub import wait as qsub_wait
from mmtbx.ions.parameters import server as metal_server
from mmtbx.monomer_library.server import server as cif_server

import os.path

xtal_params_str = """
screen = None
  .type = str
  .help = ...
screen_number = None
  .type = int
  .help = ...
"""

def _have_ligand(ligand):
  ligand = ligand.strip().upper()
  assert len(ligand) <= 3
  mon_lib_entry = cif_server().get_comp_comp_id_direct(ligand)
  return mon_lib_entry is not None

class server (object) :
  def __init__ (self) :
    import iotbx.cif
    params_path = os.path.join(os.path.split(__file__)[0],
      "crystallization_screens.cif")
    assert (os.path.isfile(params_path))
    self._cif_model = iotbx.cif.reader(file_path=params_path).model()

  def get_condition (self, screen_name, condition_id) :
    screen_name = screen_name.lower().replace(" ", "_")
    keys = self._cif_model.keys()
    data = None
    if screen_name in keys:
      data = self._cif_model[screen_name]
    else :
      for other_key in keys :
        other_data = self._cif_model[other_key]
        for name in other_data["_xtal_screen.name"]:
          if name == screen_name:
            data = other_data
            break
    if (data is None) :
      raise RuntimeError("Screen '%s' not recognized!" % screen_name)
    official_name = data["_xtal_screen.name"]
    well_ids = data["_lib_screen.well_number"]
    _id = None
    for i_well, well_id in enumerate(data["_lib_screen.well_number"]) :
      if (condition_id == well_id) :
        _id = i_well
      elif(condition_id == data["_lib_screen.condition_number"][i_well]) :
        _id = i_well
      if (_id is not None) :
        # XXX sloppy
        kwds = dict([ ("screen_name_", official_name) ] + [
          (name, data["_lib_screen."+name[:-1]][_id])
            for name in solution.__slots__[1:] ])
        return solution(**kwds)
    raise RuntimeError("Condition '%s' not found in '%s'." % (condition_id,
      screen_name))

  def generate_restraints(self, phenix_source = None):
    ligands = set(lig.strip()
                  for data in self._cif_model.values()
                  for ligs in data["_lib_screen.ligands"]
                  for lig in ligs.split(","))
    s = metal_server()
    ligands.difference_update(
      resname.strip() for resname in s.params["_lib_charge.resname"])

    cmds = ["phenix.elbow --chem={} --opt".format(lig)
            for lig in ligands
            if not _have_ligand(lig)]

    if phenix_source is None:
      phenix_source = \
        os.path.join(os.environ["HOME"], "build", "setpaths_all.csh")

    cache = os.path.join(os.path.split(__file__)[0], "cache")

    try:
      os.makedirs(cache)
    except OSError as err:
      if err.errno != errno.EEXIST:
        raise err

    if not cmds:
      return

    job_id = run_easy_qsub(
      phenix_source = phenix_source,
      where = cache,
      commands = cmds,
      code = "gen_ligs",
      host_scratch_dir = "gen_ligs"
      )

    qsub_wait(task_name = "gen_ligs",
              job_id = job_id)

    print("Done! Find generated cif restraints in {}".format(cache))

class solution (slots_getstate_setstate) :
  """
  Container for information about a specific crystallization condition (as
  defined in crystallization_screens.cif).
  """
  __slots__ = [
    "screen_name_",
    "screen_number_",
    "well_number_",
    "tube_number_",
    "condition_number_",
    "ligands_",
    "pH_",
    "condition_"
  ]
  def __init__ (self, **kwds) :
    for name in self.__slots__ :
      setattr(self, name, kwds.get(name))

  def condition_name (self) :
    return "%s: %s" % (self.screen_name_, self.condition_number_)

  def pH (self) :
    if (self.pH_ == ".") : return None
    return float(self.pH_)

  def ligands (self) :
    if (self.ligands_ == ".") :
      return []
    else :
      return self.ligands_.split(",")
