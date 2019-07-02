# -*- coding: utf-8; py-indent-offset: 2 -*-
from __future__ import absolute_import, division, print_function
from libtbx import slots_getstate_setstate

import os.path

xtal_params_str = """
screen = None
  .type = str
  .help = ...
screen_number = None
  .type = int
  .help = ...
"""

class server(object):
  def __init__(self):
    import iotbx.cif
    params_path = os.path.join(os.path.split(__file__)[0],
      "crystallization_screens.cif")
    assert (os.path.isfile(params_path))
    self._cif_model = iotbx.cif.reader(file_path=params_path).model()

  def get_condition(self, screen_name, condition_id):
    """
    Gets the crystallization conditions used within a particular well.

    Parameters
    ----------
    screen_name : str
        The name of the crystal screen used. (i.e. Crystal Screen HT)
    condition_id : int
        The condition number within that screen. (i.e. D2 or 38)

    Returns
    -------
    mmtbx.ligands.xtal_screens.solution
        The information associated with that screen.
    """
    screen_name = screen_name.lower().replace(" ", "_")
    data = None
    # NOTE: I assumed cif_model is a dict
    if screen_name in self._cif_model:
      data = self._cif_model[screen_name]
    else :
      for other_key in self._cif_model :
        other_data = self._cif_model[other_key]
        for name in other_data["_xtal_screen.name"]:
          if name == screen_name:
            data = other_data
            break
    if (data is None):
      raise RuntimeError("Screen '%s' not recognized!" % screen_name)
    official_name = data["_xtal_screen.name"]
    well_ids = data["_lib_screen.well_number"]
    _id = None
    for i_well, well_id in enumerate(data["_lib_screen.well_number"]):
      if (condition_id == well_id):
        _id = i_well
      elif(condition_id == data["_lib_screen.condition_number"][i_well]):
        _id = i_well
      if (_id is not None):
        kwds = {}
        for name in solution.__slots__[1:]:
          kwds[name] = data["_lib_screen." + name[:-1]][_id]
        kwds["screen_name_"] = official_name
        return solution(**kwds)
    raise RuntimeError("Condition '%s' not found in '%s'." % (condition_id,
      screen_name))

class solution(slots_getstate_setstate):
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
  def __init__(self, **kwds):
    for name in self.__slots__ :
      setattr(self, name, kwds.get(name))

  def condition_name(self):
    return "%s: %s" % (self.screen_name_, self.condition_number_)

  def pH(self):
    if (self.pH_ == ".") : return None
    return float(self.pH_)

  def ligands(self):
    if (self.ligands_ == "."):
      return []
    else :
      return self.ligands_.split(",")
