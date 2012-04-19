
# this will try to guess file type based on extensions.  since this will
# frequently break, it will also try every other file type if necessary,
# stopping when it finds an appropriate format.
#
# MTZ file handling is kludgy, but unfortunately there are circumstances
# where only an MTZ file will do, so it requires some extra code to work
# around the automatic behavior
#
# XXX note that there is some cross-importing from mmtbx here, but it is done
# inline, not globally

from libtbx import smart_open
from libtbx.utils import Sorry
import cPickle
import os
import re
import sys

standard_file_types = ["hkl", "ccp4_map", "xplor_map", "pdb", "cif", "phil",
  "hhr", "aln", "seq", "xml", "pkl", "txt",]

standard_file_extensions = {
  'pdb'  : ["pdb", "ent"],
  'hkl'  : ["mtz", "hkl", "sca", "cns", "xplor", "cv", "ref", "fobs"],
  'cif'  : ["cif"],
  'seq'  : ["fa", "faa", "seq", "pir", "dat", "fasta"],
  'xplor_map' : ["xplor", "map"],
  'ccp4_map'  : ["ccp4", "map"],
  'map'  : ["xplor", "map", "ccp4"],
  'phil' : ["params", "eff", "def", "phil", "param"],
  'xml'  : ["xml"],
  'pkl'  : ["pickle", "pkl"],
  'txt'  : ["txt", "log", "html", "geo"],
  'mtz'  : ["mtz"],
  'aln'  : ["aln", "ali", "clustal"],
  'hhr'  : ["hhr"],
  'img'  : ["img", "osc", "mccd", "cbf"],
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
  'aln'  : "Sequence alignment",
  'hhr'  : "HHpred alignment",
  'img'  : "Detector image",
}

supported_file_types = ["pdb","hkl","cif","pkl","seq","phil", "aln", "txt",
  "xplor_map", "ccp4_map"]

binary_types = ["hkl","ccp4_map","img","pkl"]
ascii_types = ["hkl","xplor_map","pdb","cif","phil","hhr", "aln", "seq", "xml",
"txt"]

def get_wildcard_string (format) :
  assert (format in standard_file_extensions)
  wildcards = [ "*.%s" % ext for ext in standard_file_extensions[format] ]
  wildcard_str = "%s file (%s)|%s" % (standard_file_descriptions[format],
    ", ".join(wildcards), ";".join(wildcards))
  return wildcard_str

def get_wildcard_strings (formats, include_any=True) :
  wildcards = [ get_wildcard_string(format) for format in formats ]
  if (include_any) :
    if (sys.platform == "win32") :
      wildcards.insert(0, "All files (*.*)|*.*")
    else :
      wildcards.append("All files (*.*)|*.*")
  wildcards_str = "|".join(wildcards)
  return wildcards_str

class FormatError (Sorry) :
  pass

def splitext (file_name) :
  base, ext = os.path.splitext(file_name)
  if (ext == ".gz") :
    base, ext = os.path.splitext(base)
  return base, ext

def guess_file_type (file_name, extensions=standard_file_extensions) :
  base, ext = splitext(file_name)
  if (ext == "") :
    return None
  if (ext == ".mtz") : # XXX gross
    return "hkl"
  for known_type, known_extensions in extensions.iteritems() :
    if ext[1:] in known_extensions :
      return known_type
  return None

def sort_by_file_type (file_names, sort_order=None) :
  if (sort_order is None) :
    sort_order = standard_file_types
  def _score_extension (ext) :
    for n, format in enumerate(sort_order) :
      extensions = standard_file_extensions.get(format, [])
      if (ext[1:] in extensions) :
        return len(sort_order) - n
    return 0
  def _sort (f1, f2) :
    base1, ext1 = splitext(f1)
    base2, ext2 = splitext(f2)
    s1 = _score_extension(ext1)
    s2 = _score_extension(ext2)
    if (s1 > s2) :
      return -1
    elif (s2 > s1) :
      return 1
    else :
      return 0
  file_names.sort(_sort)
  return file_names

def any_file (file_name,
              get_processed_file=False,
              valid_types=supported_file_types,
              allow_directories=False,
              force_type=None,
              input_class=None,
              raise_sorry_if_errors=False) :
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
      force_type=force_type,
      raise_sorry_if_errors=raise_sorry_if_errors)

def splitext (file_name) :
  (file_base, file_ext) = os.path.splitext(file_name)
  if file_ext in [".gz"] : # XXX: does this work for anything other than pdb?
    (base2, ext2) = os.path.splitext(file_base)
    if (ext2 in [".pdb", ".ent", ".cif"]) :
      file_ext = ext2
  return (file_base, file_ext)

class any_file_input (object) :
  __extensions__ = standard_file_extensions
  __descriptions__ = standard_file_descriptions

  def __init__ (self,
      file_name,
      get_processed_file,
      valid_types,
      force_type,
      raise_sorry_if_errors = False) :
    self.valid_types = valid_types
    self.file_name = file_name
    self.file_object = None
    self.file_type = None
    self.file_server = None
    self.file_description = None
    self._cached_file = None # XXX: used in phenix.file_reader
    self._tried_types = []
    self._errors = {}
    self.file_size = os.path.getsize(file_name)
    self.get_processed_file = get_processed_file

    (file_base, file_ext) = splitext(file_name)
    if (force_type not in [None, "None"]) :
      read_method = getattr(self, "try_as_%s" % force_type, None)
      if (read_method is None) :
        raise Sorry("Couldn't force file type to '%s' - unrecognized format." %
                    force_type)
      else :
        if (raise_sorry_if_errors) :
          try :
            read_method()
          except Exception, e :
            raise Sorry(str(e))
        else :
          read_method()
    else :
      for file_type in valid_types :
        if file_ext[1:] in self.__extensions__[file_type] :
          read_method = getattr(self, "try_as_%s" % file_type)
          self._tried_types.append(file_type)
          try :
            read_method()
          except KeyboardInterrupt :
            raise
          except FormatError, e :
            raise e
          except Exception, e :
            self._errors[file_type] = str(e)
            self.file_type = None
            self.file_object = None
          else :
            break
      if self.file_type is None :
        self.try_all_types()
    if self.file_type is not None :
      self.file_description = self.__descriptions__[self.file_type]

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
    else :
      raise ValueError("Can't parse this as a PDB file.")

  def try_as_hkl (self) :
    from iotbx.reflection_file_reader import any_reflection_file
    from iotbx.reflection_file_utils import reflection_file_server
    # XXX this is unfortunate, but unicode breaks Boost.Python extensions
    hkl_file = any_reflection_file(str(self.file_name))
    assert (hkl_file.file_type() is not None), "Not a valid reflections file."
    self.file_server = reflection_file_server(
      crystal_symmetry=None,
      force_symmetry=True,
      reflection_files=[hkl_file],
      err=sys.stderr)
    self.file_type = "hkl"
    self.file_object = hkl_file

  def try_as_cif (self) :
    import iotbx.cif
    from iotbx.reflection_file_reader import any_reflection_file
    from iotbx.reflection_file_utils import reflection_file_server
    cif_file = any_reflection_file(str(self.file_name))
    if cif_file.file_type() is not None:
      self.file_server = reflection_file_server(
        crystal_symmetry=None,
        force_symmetry=True,
        reflection_files=[cif_file],
        err=sys.stderr)
      self.file_object = cif_file
      self.file_type = "hkl"
    else:
      self.file_object = iotbx.cif.reader(file_path=str(self.file_name),
        strict=False)
      self.file_type = "cif"

  def try_as_phil (self) :
    from iotbx.phil import parse as parse_phil
    phil_object = parse_phil(file_name=self.file_name, process_includes=True)
    assert (len(phil_object.objects) > 0), "Empty parameter file."
    self.file_type = "phil"
    self.file_object = phil_object

  def try_as_seq (self) :
    from iotbx.bioinformatics import any_sequence_format
    objects, non_compliant = any_sequence_format(self.file_name)
    assert (objects is not None), "No sequence data found in file."
    assert (len(non_compliant) == 0), "Misformatted data in file."
    for seq_obj in objects :
      assert (not "-" in seq_obj.sequence)
    self.file_object = objects
#    self.try_as_txt()
#    assert len(self.file_object) != 0
#    for _line in self.file_object.splitlines() :
#      assert not _line.startswith(" ")
#      line = re.sub(" ", "", _line)
#      assert ((len(line) == 0) or
#              (line[0] == ">") or
#              (line == "*") or
#              ((line[-1] == '*') and line[:-1].isalpha()) or
#              line.isalpha())
    self.file_type = "seq"

  def try_as_hhr (self) :
    from iotbx.bioinformatics import any_hh_file
    hh_object = any_hh_file(self.file_name)
    assert (not hh_object.query in ["", None])
    self.file_object = hh_object
    self.file_type = "hhr"

  def try_as_aln (self) :
    from iotbx.bioinformatics import any_alignment_file
    aln_object = any_alignment_file(self.file_name)
    self.file_object = aln_object
    self.file_type = "aln"

  def try_as_xplor_map (self) :
    import iotbx.xplor.map
    map_object = iotbx.xplor.map.reader(file_name=str(self.file_name))
    self.file_type = "xplor_map"
    self.file_object = map_object

  def try_as_ccp4_map (self) :
    import iotbx.ccp4_map
    map_object = iotbx.ccp4_map.map_reader(file_name=str(self.file_name))
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

  def try_as_xml (self) :
    import xml.dom.minidom
    xml_in = xml.dom.minidom.parse(self.file_name)
    self.file_type = "xml"
    self.file_object = xml_in

  def try_as_img (self) :
    from iotbx.detectors import ImageFactory
    img = ImageFactory(self.file_name)
    img.read()
    self.file_type = "img"
    self.file_object = img

  def try_all_types (self) :
    for filetype in self.valid_types :
      if (filetype in self._tried_types) : continue
      read_method = getattr(self, "try_as_%s" % filetype)
      try :
        read_method()
      except KeyboardInterrupt :
        raise
      except Exception, e :
        self._errors[filetype] = str(e)
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
      return "%s%s" % (self.__descriptions__[self.file_type],
        file_size_str)

  def assert_file_type (self, expected_type) :
    if (expected_type is None) :
      return self
    elif (self.file_type == expected_type) :
      return self
    else :
      raise Sorry(("Expected file type '%s' for %s, got '%s'.  This is " +
        "almost certainly a bug; please contact the developers.") %
        (expected_type, str(self.file_name), str(self.file_type)))

  def check_file_type (self, expected_type=None, multiple_formats=()) :
    if (expected_type is not None) :
      if (self.file_type != expected_type) :
        raise Sorry(("This file format ('%s') is not supported as input for "+
          "this field; only files of type '%s' are allowed.") % (
          standard_file_descriptions.get(self.file_type, "Unknown"),
          standard_file_descriptions.get(expected_type, "Unknown")))
    else :
      assert (len(multiple_formats) > 0)
      if (not self.file_type in multiple_formats) :
        raise Sorry(("This file format ('%s') is not supported as input for "+
          "this field; only the following types are supported:\n  %s") % (
          standard_file_descriptions.get(self.file_type, "Unknown"),
          "\n  ".join([ standard_file_descriptions.get(f, "Unknown")
                        for f in multiple_formats ])))
    return self

  def show_summary (self, out=sys.stdout) :
    if (self.file_type is None) :
      print >> out, "File type could not be determined."
    else :
      print >> out, "File name: %s" % self.file_name
      print >> out, "Format: %s (%s)" % (self.file_type,
        standard_file_descriptions.get(self.file_type, "unknown"))
    if (self.file_type == "pdb") :
      print >> out, "Atoms in file: %d" % (len(self.file_object.atoms()))
      title = "\n".join(self.file_object.title_section())
      if (title != "") :
        print >> out, "Title section:"
        print >> out, title
    elif (self.file_type == "hkl") :
      for array in self.file_server.miller_arrays :
        print >> out, ""
        array.show_comprehensive_summary(f=out)
    elif (self.file_type == "ccp4_map") :
      self.file_object.show_summary(out)

# mimics any_file, but without parsing - will instead guess the file type from
# the extension.  for most output files produced by cctbx/phenix this is
# relatively safe.
def any_file_fast (file_name,
              get_processed_file=False,
              valid_types=supported_file_types,
              allow_directories=False,
              force_type=None,
              input_class=None) :
  assert (not get_processed_file) and (force_type is None)
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
      input_class = any_file_fast_input
    return input_class(file_name=file_name,
      valid_types=valid_types)

class any_file_fast_input (object) :
  __extensions__ = standard_file_extensions
  __descriptions__ = standard_file_descriptions
  def __init__ (self, file_name, valid_types) :
    self.valid_types = valid_types
    self.file_name = file_name
    self.file_object = None
    self.file_type = None
    self.file_server = None
    self.file_description = None
    file_base, file_ext = splitext(file_name)
    for file_type in valid_types :
      if file_ext[1:] in self.__extensions__[file_type] :
        self.file_type = file_type
    if self.file_type is not None :
      self.file_description = self.__descriptions__[self.file_type]
    # XXX this is a huge hole in this method - for reflection files, the
    # PHENIX GUI displays the specific format, which requires actually reading
    # in the file.  the stub class below will work for obvious formats.
    if (self.file_type == "hkl") :
      class fake_data_object (object) :
        def __init__ (self, ext) :
          self.ext = ext
        def file_type (self) :
          if (self.ext == ".mtz") : return "CCP4 MTZ"
          elif (self.ext == ".sca") : return "Scalepack"
          else : return "Data (other)"
      self.file_object = fake_data_object(file_ext)

class directory_input (object) :
  def __init__ (self, dir_name) :
    self.file_name = dir_name
    self.file_object = dircache.listdir(dir_name)
    self.file_server = None
    self.file_type = "dir"
    self.file_size = os.path.getsize(dir_name)

  def file_info (self, show_file_size=False) :
    return "Folder"

class group_files (object) :
  def __init__ (self,
                file_names,
                template_format="pdb",
                group_by_directory=True) :
    import iotbx.pdb
    self.file_names = file_names
    self.grouped_files = []
    self.ungrouped_files = []
    self.ambiguous_files = []
    templates = []
    other_files = []
    template_dirs = []
    for file_name in file_names :
      file_type = guess_file_type(file_name)
      if (file_type == template_format) :
        base, ext = splitext(file_name)
        templates.append(base)
        template_dirs.append(os.path.dirname(file_name))
        self.grouped_files.append([file_name])
      else :
        other_files.append(file_name)
    if (len(templates) == 0) :
      raise Sorry("Can't find any files of the expected format ('%s')." %
        template_format)
    if (len(set(templates)) != len(templates)) :
      raise Sorry("Multiple files with identical root names.")
    for file_name in other_files :
      group_name = find_closest_base_name(
        file_name=file_name,
        base_name=splitext(file_name)[0],
        templates=templates)
      if (group_name == "") :
        self.ambiguous_files.append(file_name)
      elif (group_name is not None) :
        i = templates.index(group_name)
        self.grouped_files[i].append(file_name)
      else :
        if group_by_directory :
          dir_name = os.path.dirname(file_name)
          group_name = find_closest_base_name(
            file_name=dir_name,
            base_name=dir_name,
            templates=template_dirs)
          if (group_name == "") :
            self.ambiguous_files.append(file_name)
          elif (group_name is not None) :
            i = template_dirs.index(group_name)
            self.grouped_files[i].append(file_name)
          else :
            self.ungrouped_files.append(file_name)
        else :
          self.ungrouped_files.append(file_name)

def find_closest_base_name (file_name, base_name, templates) :
  groups = []
  for base in templates :
    if file_name.startswith(base) or base.startswith(base_name) :
      groups.append(base)
  if (len(groups) == 1) :
#    print file_name, groups[0]
    return groups[0]
  elif (len(groups) > 1) :
#    print file_name, groups
    prefix_len = [ os.path.commonprefix([g, file_name]) for g in groups ]
    max_common_prefix = max(prefix_len)
    if (prefix_len.count(max_common_prefix) > 1) :
      return ""
    else :
      return groups[ prefix_len.index(max_common_prefix) ]
  return None

#---end
