
from __future__ import absolute_import, division, print_function
import os.path

from libtbx.utils import to_str

class watch_model(object):
  """
  Simple filesystem-only automatic loader for PDB files.  Once invoked with
  the specified file name, it will inspect the file mtime at regular intervals
  and reload if modification is detected.
  """
  def __init__(self, file_name, timeout=2000, model_imol=None):
    self.file_name = file_name
    self.model_imol = model_imol
    self.last_mtime = 0
    import gobject
    gobject.timeout_add(timeout, self.check_file)

  def check_file(self):
    import coot # import dependency
    if os.path.exists(self.file_name):
      file_mtime = os.path.getmtime(self.file_name)
      if file_mtime > self.last_mtime :
        if self.model_imol is not None :
          clear_and_update_molecule_from_file(self.model_imol,
            to_str(self.file_name))
        else :
          self.model_imol = read_pdb(to_str(self.file_name))
        self.last_mtime = file_mtime


class watch_model_and_maps(object):
  """
  Simple filesystem-only automatic loader for PDB and MTZ map coefficient
  files.  Identical to the watch_model class in most respects.
  """
  def __init__(self, file_prefix, timeout=2000, model_imol=None,
      map_imols=()):
    self.file_prefix = file_prefix
    self.pdb_file = self.file_prefix + ".pdb"
    self.map_file = self.file_prefix + ".mtz"
    self.model_imol = model_imol
    self.map_imols = map_imols
    self.last_mtime = 0
    import gobject
    gobject.timeout_add(timeout, self.check_files)

  def check_files(self):
    if os.path.exists(self.pdb_file) and os.path.exists(self.map_file):
      pdb_mtime = os.path.getmtime(self.pdb_file)
      map_mtime = os.path.getmtime(self.map_file)
      if (pdb_mtime > self.last_mtime) and (map_mtime > self.last_mtime):
        self.update_files()
        self.last_mtime = max(pdb_mtime, map_mtime)

  def update_files(self):
    import coot # import dependency
    if (self.model_imol is not None):
      clear_and_update_molecule_from_file(self.model_imol,
                                          to_str(self.pdb_file))
    else :
      self.model_imol = read_pdb(to_str(self.pdb_file))
    if (len(self.map_imols) > 0):
      set_colour_map_rotation_for_map(0)
      for imol in self._map_mols :
        close_molecule(imol)
    else :
      set_colour_map_rotation_for_map(10)
    map_imol = auto_read_make_and_draw_maps(to_str(self.map_file))
    if (isinstance(map_imol, int)):
      # XXX this may be dangerous, but auto_read_make_and_draw_maps only returns
      # the last imol
      self.map_imols = [ map_imol - 1, map_imol ]
    else :
      self.map_imols = map_imol
    set_imol_refinement_map(self.map_imols[0])
    set_scrollable_map(self.map_imols[0])
