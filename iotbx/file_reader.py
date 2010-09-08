
# this will try to guess file type based on extensions.  since this will
# frequently break, it will also try every other file type if necessary,
# stopping when it finds an appropriate format.

# XXX note that there is some cross-importing from mmtbx here, but it is done
# inline, not globally

import sys, os, re, string
from libtbx import smart_open
from libtbx.utils import Sorry
import cPickle

standard_file_types = ["hkl", "ccp4_map", "xplor_map", "pdb", "cif", "phil",
  "seq", "xml", "pkl", "txt"]

standard_file_extensions = {
  'pdb'  : ["pdb", "ent"],
  'hkl'  : ["mtz", "hkl", "sca", "cns", "xplor", "cv", "ref", "fobs"],
  'cif'  : ["cif"],
  'seq'  : ["fa", "faa", "seq", "pir", "dat"],
  'xplor_map' : ["xplor", "map"],
  'ccp4_map'  : ["ccp4", "map"],
  'phil' : ["params", "eff", "def", "phil"],
  'xml'  : ["xml"],
  'pkl'  : ["pickle", "pkl"],
  'txt'  : ["txt", "log", "html"],
  'mtz'  : ["mtz"],
}
compression_extensions = ["gz", "Z", "bz2", "zip"]

standard_file_descriptions = {
  'pdb'  : "Model",
  'hkl'  : "Reflections",
  'cif'  : "Restraints",
  'seq'  : "Sequence",
  'xplor_map'  : "XPLOR map",
  'ccp4_map' : "CCP4 map",
  'phil' : "Parameters",
  'xml'  : "XML",
  'pkl'  : "Python pickle",
  'txt'  : "Text",
  'mtz'  : "Reflections (MTZ)",
}

supported_file_types = ["pdb","hkl","cif","pkl","seq","phil", "txt",
  "xplor_map", "ccp4_map"]

class FormatError (Sorry) :
  pass

def guess_file_type (file_name, extensions=standard_file_extensions) :
  base, ext = os.path.splitext(file_name)
  for known_type, known_extensions in extensions.iteritems() :
    if ext[1:] in known_extensions :
      return known_type
  return None

def any_file (file_name,
              get_processed_file=False,
              valid_types=supported_file_types,
              allow_directories=False,
              force_type=None,
              input_class=None) :
  if not os.path.exists(file_name) :
    raise Sorry("Couldn't find the file %s" % file_name)
  elif os.path.isdir(file_name) :
    if not allow_directories :
      raise Sorry("This application does not support folders as input.")
    else :
      return directory_input(file_name)
  elif not os.path.isfile(file_name) :
    raise Sorry("%s is not a valid file.")
  else :
    if input_class is None :
      input_class = any_file_input
    return input_class(file_name=file_name,
      get_processed_file=get_processed_file,
      valid_types=valid_types,
      force_type=force_type)

class any_file_input (object) :
  __extensions = standard_file_extensions
  __descriptions = standard_file_descriptions

  def __init__ (self, file_name, get_processed_file, valid_types, force_type) :
    self.valid_types = valid_types
    self.file_name = file_name
    self.file_object = None
    self.file_type = None
    self.file_server = None
    self.file_description = None
    self._cached_file = None # XXX: used in phenix.file_reader
    self.file_size = os.path.getsize(file_name)
    self.get_processed_file = get_processed_file

    (file_base, file_ext) = os.path.splitext(file_name)
    if file_ext in [".gz"] : # XXX: does this work for anything other than pdb?
      (base2, ext2) = os.path.splitext(file_base)
      if ext2 != "" :
        file_ext = ext2
    if force_type is not None :
      read_method = getattr(self, "try_as_%s" % force_type, None)
      if read_method is None :
        raise Sorry("Couldn't force file type to '%s' - unrecognized format." %
                    force_type)
      else :
        read_method()
    else :
      for file_type in valid_types :
        if file_ext[1:] in self.__extensions[file_type] :
          read_method = getattr(self, "try_as_%s" % file_type)
          try :
            read_method()
          except KeyboardInterrupt :
            raise
          except FormatError, e :
            raise e
          except Exception, e :
            self.file_type = None
            self.file_object = None
          else :
            break
      if self.file_type is None :
        self.try_all_types()
    if self.file_type is not None :
      self.file_description = self.__descriptions[self.file_type]

  def try_as_pdb (self) :
    from iotbx.pdb import is_pdb_file
    if is_pdb_file(self.file_name) :
      from iotbx.pdb import input as pdb_input
      from scitbx.array_family import flex
      raw_records = flex.std_string()
      pdb_file = smart_open.for_reading(file_name=self.file_name)
      raw_records.extend(flex.split_lines(pdb_file.read()))
      structure = pdb_input(source_info=None, lines=raw_records)
      self.file_type = "pdb"
      self.file_object = structure

  def try_as_hkl (self) :
    from iotbx.reflection_file_reader import any_reflection_file
    from iotbx.reflection_file_utils import reflection_file_server
    try :
      hkl_file = any_reflection_file(self.file_name)
    except Exception, e :
      print e
      raise
    assert hkl_file.file_type() is not None
    self.file_server = reflection_file_server(
      crystal_symmetry=None,
      force_symmetry=True,
      reflection_files=[hkl_file],
      err=sys.stderr)
    self.file_type = "hkl"
    self.file_object = hkl_file

  def try_as_cif (self) :
    from mmtbx.monomer_library import server
    cif_object = server.read_cif(file_name=self.file_name)
    assert len(cif_object) != 0
    self.file_type = "cif"
    self.file_object = cif_object

  def try_as_phil (self) :
    from iotbx.phil import parse as parse_phil
    phil_object = parse_phil(file_name=self.file_name, process_includes=True)
    self.file_type = "phil"
    self.file_object = phil_object

  def try_as_seq (self) :
    self.try_as_txt()
    assert len(self.file_object) != 0
    for _line in self.file_object.splitlines() :
      assert not _line.startswith(" ")
      line = re.sub(" ", "", _line)
      assert (len(line) == 0 or
              line[0] == ">" or
              (line[-1] == '*' and line[:-1].isalpha()) or
              line.isalpha())
    self.file_type = "seq"

  def try_as_xplor_map (self) :
    import iotbx.xplor.map
    map_object = iotbx.xplor.map.reader(file_name=self.file_name)
    self.file_type = "xplor_map"
    self.file_object = map_object

  def try_as_ccp4_map (self) :
    import iotbx.ccp4_map
    map_object = iotbx.ccp4_map.map_reader(file_name=self.file_name)
    self.file_type = "ccp4_map"
    self.file_object = map_object

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
      except KeyboardInterrupt :
        raise
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
      return "%s%s" % (self.__descriptions[self.file_type],
        file_size_str)

  def assert_file_type (self, expected_type) :
    if (expected_type is None) :
      return None
    elif (self.file_type == expected_type) :
      return True
    else :
      raise Sorry(("Expected file type '%s' for %s, got '%s'.  This is " +
        "almost certainly a bug; please contact the developers.") %
        (str(self.file_name), expected_type, str(self.file_type)))

class directory_input (object) :
  def __init__ (self, dir_name) :
    self.file_name = dir_name
    self.file_object = dircache.listdir(dir_name)
    self.file_server = None
    self.file_type = "dir"
    self.file_size = os.path.getsize(dir_name)

  def file_info (self, show_file_size=False) :
    return "Folder"

#---end
