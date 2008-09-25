
import sys, os, re
from mmtbx.monomer_library import server
from iotbx.pdb import is_pdb_file
from iotbx.reflection_file_reader import any_reflection_file
from libtbx.utils import Sorry

standard_file_extensions = {
  'pdb'  : ["pdb", "ent"],
  'hkl'  : ["mtz", "hkl", "sca", "cns", "xplor"],
  'cif'  : ["cif"],
  'seq'  : ["fa", "seq", "pir", "dat"],
  'map'  : ["map", "ccp4"],
  'phil' : ["params", "eff", "def", "phil"]
}

standard_file_descriptions = {
  'pdb'  : "Model",
  'hkl'  : "Reflections",
  'cif'  : "Restraints",
  'seq'  : "Sequence",
  'map'  : "Map",
  'phil' : "Parameters"
}

def get_file_extension (file_name) :
  base_file = file_name

  if len(file_name) < 3 :
    return ""
  if file_name[-3:] == ".gz" :
    base_file = file_name[:-3]
  file_ext = re.sub(".*\.", "", base_file)
  return file_ext

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
    self.file_size = os.path.getsize(file_name)
    self.get_processed_file = get_processed_file

    file_ext = get_file_extension(file_name)
    try :
      if file_ext in standard_file_extensions['pdb'] :
        self.try_as_pdb()
      elif file_ext in standard_file_extensions['hkl'] :
        self.try_as_hkl()
      elif file_ext in standard_file_extensions['phil'] :
        self.try_as_phil()
      elif file_ext in standard_file_extensions['seq'] :
        self.try_as_seq()
      elif file_ext in standard_file_extensions['cif'] :
        self.try_as_cif()
    except Exception :
      pass

    if not self.file_type :
      self.try_all_types()

  def try_as_pdb (self) :
    if self.file_type != None : return False
    if is_pdb_file(self.file_name) :
      self.file_type = "pdb"
      if self.get_processed_file :
        pass

  def try_as_hkl (self) :
    if self.file_type != None : return False
    try :
      hkl_file = any_reflection_file(self.file_name)
      self.file_type = "hkl"
      self.file_object = hkl_file
    except Exception :
      pass

  def try_as_cif (self) :
    if self.file_type != None : return False
    try :
      cif_object = server.read_cif(file_name=self.file_name)
      self.file_type = "cif"
      self.file_object = cif_object
    except Exception :
      pass

  def try_as_phil (self) :
    if self.file_type != None : return False
    try :
      phil_object = iotbx.phil.parse(file_name=self.file_name)
      self.file_type = "phil"
      self.file_object = phil_object
    except Exception :
      pass

  def try_as_seq (self) :
    if self.file_type != None : return False
    pass

  def try_as_map (self) :
    if self.file_type != None : return False
    pass

  def try_all_types (self) :
    self.try_as_pdb()
    self.try_as_hkl()
    self.try_as_cif()
    self.try_as_phil()
    self.try_as_seq()
    self.try_as_map()

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
