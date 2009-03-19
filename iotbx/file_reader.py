
# this will try to guess file type based on extensions.  since this will
# frequently break, it will also try every other file type if necessary,
# stopping when it finds an appropriate format.

# TODO: map files

import sys, os, re, string
from mmtbx.monomer_library import server
from iotbx.phil import parse as parse_phil
from iotbx.pdb import is_pdb_file
from iotbx.reflection_file_reader import any_reflection_file
from iotbx.reflection_file_utils import reflection_file_server
from libtbx.utils import Sorry
import cPickle

standard_file_types = ["hkl", "map", "pdb", "cif", "phil", "seq", "xml", "pkl",
                       "txt"]

standard_file_extensions = {
  'pdb'  : ["pdb", "ent"],
  'hkl'  : ["mtz", "hkl", "sca", "cns", "xplor", "cv", "ref"],
  'cif'  : ["cif"],
  'seq'  : ["fa", "faa", "seq", "pir", "dat"],
  'map'  : ["map", "ccp4"],
  'phil' : ["params", "eff", "def", "phil"],
  'xml'  : ["xml"],
  'pkl'  : ["pickle", "pkl"],
  'txt'  : ["txt", "log", "html"],
}
compression_extensions = ["gz", "Z", "bz2", "zip"]

standard_file_descriptions = {
  'pdb'  : "Model",
  'hkl'  : "Reflections",
  'cif'  : "Restraints",
  'seq'  : "Sequence",
  'map'  : "Map",
  'phil' : "Parameters",
  'xml'  : "XML",
  'pkl'  : "Python pickle",
  'txt'  : "Text"
}

def any_file (file_name,
              get_processed_file=False,
              valid_types=["pdb","hkl","cif","pkl","seq","phil", "txt"],
              allow_directories=False) :
  if not os.path.exists(file_name) :
    raise Sorry("Couldn't find the file %s" % file_name)
  elif os.path.isdir(file_name) :
    if not allow_directories :
      raise Sorry("This application does not support folders as input.")
    else :
      return _dir_input(file_name)
  elif not os.path.isfile(file_name) :
    raise Sorry("%s is not a valid file.")
  else :
    return _any_file(file_name, get_processed_file, valid_types)

class _any_file (object) :
  def __init__ (self, file_name, get_processed_file, valid_types) :
    self.valid_types = valid_types
    self.file_name = file_name
    self.file_object = None
    self.file_type = None
    self.file_server = None
    self.file_size = os.path.getsize(file_name)
    self.get_processed_file = get_processed_file

    (file_base, file_ext) = os.path.splitext(file_name)
    for file_type in valid_types :
      if file_ext[1:] in standard_file_extensions[file_type] :
        try :
          read_method = getattr(self, "try_as_%s" % file_type)
          read_method()
        except Exception, e :
          print e
          self.file_type = None
          self.file_object = None
        else :
          break
    if self.file_type is None :
      self.try_all_types()

  def try_as_pdb (self) :
    if is_pdb_file(self.file_name) :
      self.file_type = "pdb"

  def try_as_hkl (self) :
    hkl_file = any_reflection_file(self.file_name)
    assert hkl_file.file_type() is not None
    self.file_server = reflection_file_server(
      crystal_symmetry=None,
      force_symmetry=True,
      reflection_files=[hkl_file],
      err=sys.stderr)
    self.file_type = "hkl"
    self.file_object = hkl_file

  def try_as_cif (self) :
    cif_object = server.read_cif(file_name=self.file_name)
    assert len(cif_object) != 0
    self.file_type = "cif"
    self.file_object = cif_object

  def try_as_phil (self) :
    phil_object = parse_phil(file_name=self.file_name, process_includes=True)
    self.file_type = "phil"
    self.file_object = phil_object

  def try_as_seq (self) :
    self.try_as_txt()
    assert len(self.file_object) != 0
    for _line in self.file_object.splitlines() :
      line = _line.rstrip()
      assert (len(line) == 0 or
              line[0] == ">" or
              (line[-1] == '*' and line[:-1].isalpha()) or
              line.isalpha())
    self.file_type = "seq"

  def try_as_map (self) :
    raise Sorry("Map input not currently supported.")

  def try_as_pkl (self) :
    pkl_object = cPickle.load(open(self.file_name, "rb"))
    self.file_type = "pkl"
    self.file_object = pkl_object

  def try_as_txt (self) :
    file_as_string = open(self.file_name).read()
    file_as_ascii = file_as_string.decode("ascii")
    self.file_type = "txt"
    self.file_object = file_as_string

  def try_all_types (self) :
    for filetype in self.valid_types :
      read_method = getattr(self, "try_as_%s" % filetype)
      try :
        read_method()
      except Exception, e :
        self.file_type = None
        self.file_object = None
        continue
      else :
        if self.file_type is not None :
          break

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

class _dir_input (object) :
  def __init__ (self, dir_name) :
    self.file_name = dir_name
    self.file_object = dircache.listdir(dir_name)
    self.file_server = None
    self.file_type = "dir"
    self.file_size = os.path.getsize(dir_name)

  def file_info (self, show_file_size=False) :
    return "Folder"

#---end
