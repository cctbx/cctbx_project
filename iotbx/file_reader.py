
import sys, os, re
from mmtbx.monomer_library import server
from iotbx.pdb import is_pdb_file
from iotbx.reflection_file_reader import any_reflection_file
from iotbx.reflection_file_utils import reflection_file_server
from libtbx.utils import Sorry

standard_file_types = ["hkl", "map", "pdb", "cif", "phil", "seq"]

standard_file_extensions = {
  'pdb'  : ["pdb", "ent"],
  'hkl'  : ["mtz", "hkl", "sca", "cns", "xplor"],
  'cif'  : ["cif"],
  'seq'  : ["fa", "faa", "seq", "pir", "dat"],
  'map'  : ["map", "ccp4"],
  'phil' : ["params", "eff", "def", "phil"]
}
compression_extensions = ["gz", "Z", "bz2", "zip"]

standard_file_descriptions = {
  'pdb'  : "Model",
  'hkl'  : "Reflections",
  'cif'  : "Restraints",
  'seq'  : "Sequence",
  'map'  : "Map",
  'phil' : "Parameters"
}

class any_file (object) :
  def __init__ (
      self,
      file_name,
      get_processed_file=False,
      valid_types=["pdb","hkl","cif","seq","map","phil"]) :

    if not os.path.isfile(file_name) :
      raise Sorry("%s is not a valid file.")

    self.file_name = file_name
    self.file_object = None
    self.file_type = None
    self.file_server = None
    self.file_size = os.path.getsize(file_name)
    self.get_processed_file = get_processed_file

    (file_base, file_ext) = os.path.splitext(file_name)
    try :
      for file_type in standard_file_types :
        if file_ext[1:] in standard_file_extensions[file_type] :
          read_method = getattr(self, "try_as_%s" % file_type)
          read_method()
          break
    except Exception, e :
      self.file_type = None

    if self.file_type is None :
      self.try_all_types()

  def try_as_pdb (self) :
    if self.file_type is not None : return False
    if is_pdb_file(self.file_name) :
      self.file_type = "pdb"

  def try_as_hkl (self) :
    if self.file_type is not None : return False
    hkl_file = any_reflection_file(self.file_name)
    self.file_server = reflection_file_server(
      crystal_symmetry=None,
      force_symmetry=True,
      reflection_files=[hkl_file],
      err=sys.stderr)
    self.file_type = "hkl"
    self.file_object = hkl_file

  def try_as_cif (self) :
    if self.file_type is not None : return False
    cif_object = server.read_cif(file_name=self.file_name)
    self.file_type = "cif"
    self.file_object = cif_object

  def try_as_phil (self) :
    if self.file_type is not None : return False
    phil_object = iotbx.phil.parse(file_name=self.file_name)
    self.file_type = "phil"
    self.file_object = phil_object

  def try_as_seq (self) :
    if self.file_type is not None : return False
    self.file_type = "seq"

  def try_as_map (self) :
    if self.file_type is not None : return False
    pass

  def try_all_types (self) :
    for filetype in standard_file_types :
      read_method = getattr(self, "try_as_%s" % filetype)
      try :
        read_method()
      except Exception, e :
        pass

  def file_info (self, show_file_size=True) :
    file_size_str = ""
    if show_file_size :
      file_size = self.file_size
      if file_size > 10000000 :
        file_size_str = " (%.1f MB)" % (self.file_size / 1000000.0)
      elif file_size > 1000000 :
        file_size_str = " (%.2f MB)" % (self.file_size / 1000000.0)
      elif file_size > 100000 :
        file_size_str = " (%d KB)" % (self.file_size / 1000.0)
      elif file_size > 1000 :
        file_size_str = " (%.1f KB)" % (self.file_size / 1000.0)
      else :
        file_size_str = " (%d B)" % self.file_size
    if self.file_type == None :
      return "Unknown file%s" % file_size_str
    else :
      return "%s%s" % (standard_file_descriptions[self.file_type],
        file_size_str)

#---end
