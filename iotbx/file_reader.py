"""
Generic file input module, used in Phenix GUI and elsewhere.  This trades some
loss of efficiency for a simplified API for loading any file type more or less
automatically.  It will first try to guess the format based on the extension;
if this fails, it will then try other formats.  This is used on the command
line and the Phenix GUI to process bulk file input.  In most other cases a
specific file type is desired, and the force_type argument will ensure that
only this format is attempted.

Examples
--------
>>> from iotbx.file_reader import any_file
>>> input_file = any_file(sys.argv[1:])
>>> file_data = input_file.file_object

>>> pdb_in = any_file("model.pdb", force_type="pdb")
>>> pdb_in.assert_file_type("pdb")
>>> hierarchy = pdb_in.file_object.hierarchy

>>> mtz_in = any_file("data.mtz", force_type="hkl")
>>> miller_arrays = mtz_in.file_server.miller_arrays
"""
from __future__ import absolute_import, division, print_function

# MTZ file handling is kludgy, but unfortunately there are circumstances
# where only an MTZ file will do, so it requires some extra code to work
# around the automatic behavior

from libtbx.utils import Sorry, to_str
from six.moves import cPickle as pickle
import os
import sys
import six

standard_file_types = ["hkl", "ccp4_map", "xplor_map", "pdb", "cif", "phil",
  "hhr", "ncs", "aln", "seq", "xml", "pkl", "txt",]

standard_file_extensions = {
  'pdb'  : ["pdb", "ent"],
  'hkl'  : ["mtz", "hkl", "sca", "cns", "xplor", "cv", "ref", "fobs"],
  'cif'  : ["cif", "mmcif"],
  'seq'  : ["fa", "faa", "seq", "pir", "dat", "fasta"],
  'xplor_map' : ["xplor", "map"],
  'ccp4_map'  : ["ccp4", "map", "mrc"],
  'map'  : ["xplor", "map", "ccp4"],
  'phil' : ["params", "eff", "def", "phil", "param"],
  'xml'  : ["xml"],
  'pkl'  : ["pickle", "pkl"],
  'txt'  : ["txt", "log", "html", "geo"],
  'mtz'  : ["mtz"],
  'aln'  : ["aln", "ali", "clustal"],
  'a3m'  : ["a3m"],
  'hhr'  : ["hhr"],
  'ncs'  : ["ncs","ncs_spec"],
  'img'  : ["img", "osc", "mccd", "cbf", "nxs", "h5", "hdf5"],
  # XXX these are not supported by this module, but are used in Phenix, so we
  # need to be able to include them in GUI tools in wxtbx.
  'smi' : ['smi'],
  'sdf' : ['sdf'],
  'rosetta' : ['gz'],
}
compression_extensions = ["gz", "Z", "bz2", "zip"] # gz and bz2 work with maps... gz, Z and bz2 work with models

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
  'a3m'  : "MSA (a3m)",
  'hhr'  : "HHpred alignment",
  'ncs'  : "NCS information file",
  'img'  : "Detector image",
  # XXX see comment above
  'smi' : "SMILES",
  'sdf' : "Structure data file",
  'rosetta' : "Rosetta fragments",
}

supported_file_types = ["pdb","hkl","cif","pkl","ncs","seq","phil",
  "aln", "a3m", "txt", "xplor_map", "ccp4_map"]

binary_types = ["hkl","ccp4_map","img","pkl"]
ascii_types = ["hkl","xplor_map","pdb","cif","phil","hhr", "ncs", "aln", "a3m",
   "seq", "xml", "txt"]

# Try files with these extensions only with their associated file types
extensions_absolutely_defining_type = ['ccp4','mrc']

def get_wildcard_string(format):
  assert (format in standard_file_extensions), format
  wildcards = [ "*.%s" % ext for ext in standard_file_extensions[format] ]
  # Add wildcard without breaking other features by having cif in two places
  if format == 'pdb':
    cif_format = "*.%s" % "cif"
    if not cif_format in wildcards:
      wildcards.append(cif_format)

  wildcard_str = "%s file (%s)|%s" % (standard_file_descriptions[format],
    ", ".join(wildcards), ";".join(wildcards))
  return wildcard_str

def get_wildcard_strings(formats, include_any=True):
  wildcards = [ get_wildcard_string(format) for format in formats ]
  if (include_any):
    if (sys.platform != "darwin"):
      wildcards.insert(0, "All files (*.*)|*.*")
    else :
      wildcards.append("All files (*.*)|*.*")
  wildcards_str = "|".join(wildcards)
  return wildcards_str

class FormatError(Sorry):
  pass

def strip_shelx_format_extension(file_name):
  if (file_name.endswith("=hklf3") or file_name.endswith("=hklf4") or
      file_name.endswith("=amplitudes") or file_name.endswith("=intensities")):
    file_name = "".join(file_name.split("=")[:-1])
  return file_name

def splitext(file_name):
  """
  Args:
    file_name: A plain text string of the file path

  Returns: A tuple of the base filename, the file format extension, and possibly a compression extension
  """

  folder_name, file_only = os.path.split(file_name)
  period_split = file_only.split(".")
  compressed = period_split[-1] in compression_extensions
  if compressed:
    file_ext = '.'+period_split[-2]
    compress_ext = '.'+period_split[-1]
    file_base = file_only.replace('.' + file_ext + '.' + compress_ext, '')
  else:
    file_ext = '.'+period_split[-1]
    compress_ext = None
    file_base = file_only.replace('.' + file_ext, '')
  file_base = os.path.join(folder_name, file_base)

  return (file_base, file_ext, compress_ext)


def guess_file_type(file_name, extensions=standard_file_extensions):
  base, ext, compress_ext = splitext(file_name)
  if (ext == ""):
    return None
  if (ext == ".mtz") : # XXX gross
    return "hkl"
  for known_type, known_extensions in six.iteritems(extensions):
    if ext[1:] in known_extensions :
      return known_type
  return None

def sort_by_file_type(file_names, sort_order=None):
  if (sort_order is None):
    sort_order = standard_file_types
  def _score_extension(ext):
    for n, format in enumerate(sort_order):
      extensions = standard_file_extensions.get(format, [])
      if (ext[1:] in extensions):
        return len(sort_order) - n
    return 0
  def _sort(f1, f2):
    base1, ext1, compress_ext1  = splitext(f1)
    base2, ext2, compress_ext2 = splitext(f2)
    s1 = _score_extension(ext1)
    s2 = _score_extension(ext2)
    if (s1 > s2):
      return -1
    elif (s2 > s1):
      return 1
    else :
      return 0
  from functools import cmp_to_key
  file_names.sort(key = cmp_to_key(_sort))
  return file_names

def any_file(file_name,
              get_processed_file=False,
              valid_types=supported_file_types,
              extensions_absolutely_defining_type=
                  extensions_absolutely_defining_type,
              allow_directories=False,
              force_type=None,
              input_class=None,
              raise_sorry_if_errors=False,
              raise_sorry_if_not_expected_format=False):
  """
  Main input method, wrapper for any_file_input class.

  :param file_name: path to file (relative or absolute)
  :param get_processed_file: TODO
  :param valid_types: file types to consider
  :param allow_directories: process directory if given as file_name
  :param force_type: read as this format, don't try any others
  :param input_class: optional substitute for any_file_input, with additional
    parsers
  :param raise_sorry_if_errors: raise a Sorry exception if parsing fails (used
    with force_type)
  :param raise_sorry_if_not_expected_format: raise a Sorry exception if the
    file extension does not match the parsed file type
  :param extensions_absolutely_defining_type: if the file has one of these
    extensions, only try the associated file type
  :returns: any_file_input object, or an instance of the input_class param
  """
  file_name_raw = file_name
  file_name = strip_shelx_format_extension(file_name)
  if (file_name != file_name_raw) and (force_type is None):
    force_type = "hkl"
  # See if extension is in extensions_absolutely_defining_type
  _, ext = os.path.splitext(file_name_raw)
  ext = ext[1:]  # remove . in .ccp4
  if (force_type is None) and extensions_absolutely_defining_type and (
      ext in extensions_absolutely_defining_type):
    for ft in supported_file_types:
      if ext in standard_file_extensions[ft]:
        force_type = ft
        break
    assert force_type is not None

  if not os.path.exists(file_name):
    raise Sorry("Couldn't find the file %s" % file_name)
  elif os.path.isdir(file_name):
    if not allow_directories :
      raise Sorry("This application does not support folders as input.")
    else :
      return directory_input(file_name)
  elif not os.path.isfile(file_name):
    raise Sorry("%s is not a valid file.")
  else :
    if input_class is None :
      input_class = any_file_input
    return input_class(file_name=file_name_raw,
      get_processed_file=get_processed_file,
      valid_types=valid_types,
      force_type=force_type,
      raise_sorry_if_errors=raise_sorry_if_errors,
      raise_sorry_if_not_expected_format=raise_sorry_if_not_expected_format)



class any_file_input(object):
  """
  Container for file data of any supported type.  Usually obtained via the
  any_file() function rather than being instantiated directly.
  """

  __extensions__ = standard_file_extensions
  __descriptions__ = standard_file_descriptions

  def __init__(self,
      file_name,
      get_processed_file,
      valid_types,
      force_type,
      raise_sorry_if_errors=False,
      raise_sorry_if_not_expected_format=False) : # XXX should probably be True
    self.valid_types = valid_types
    self._file_name = file_name
    self._file_object = None
    self._file_type = None
    self._file_server = None
    self.file_description = None
    self._cached_file = None # XXX: used in phenix.file_reader
    self._tried_types = []
    self._errors = {}
    file_name_clean = strip_shelx_format_extension(file_name)
    self.file_size = os.path.getsize(file_name_clean)
    self.get_processed_file = get_processed_file
    file_base, file_ext, compress_ext = splitext(file_name)
    file_ext = file_ext.lower()
    if (force_type not in [None, "None"]):
      read_method = getattr(self, "_try_as_%s" % force_type, None)
      if (read_method is None):
        raise Sorry("Couldn't force file type to '%s' - unrecognized format." %
                    force_type)
      else :
        if (raise_sorry_if_errors):
          try :
            read_method()
          except Exception as e :
            raise Sorry("Couldn't read '%s' as file type '%s': %s" %
              (file_name, force_type, str(e)))
        else :
          read_method()
    else :
      # XXX this is probably not the best way to do this - if the file format
      # is obviously something we don't want, this should be determined first
      # isntead of trying the limited set of parsers which won't work.
      for file_type in valid_types :
        if ((file_ext[1:] in self.__extensions__[file_type]) and
            (not file_ext in [".txt"])):
          read_method = getattr(self, "_try_as_%s" % file_type)
          self._tried_types.append(file_type)
          try :
            read_method()
          except KeyboardInterrupt :
            raise
          except FormatError as e :
            raise e
          except Exception as e :
            # XXX we need to be a little careful about what extensions we
            # do this for - they are not necessarily all unambiguous!
            if ((raise_sorry_if_not_expected_format) and
                (file_ext in [".pdb",".mtz",".cif",".sca",".xml",".phil"])):
              raise Sorry("File format error:\n" + str(e))
            self._errors[file_type] = str(e)
            self._file_type = None
            self._file_object = None
          else :
            break
      if self._file_type is None :
        self.try_all_types()
    if self._file_type is not None :
      self.file_description = self.__descriptions__[self.file_type]

  @property
  def file_name(self):
    return self._file_name

  @property
  def file_type(self):
    """
    Return a string representing the generic data type, for example 'pdb' or
    'hkl'.  Note that this is not necessarily the same as the underlying
    format, for example 'pdb' can mean either PDB or mmCIF format, and 'hkl'
    could mean MTZ, CIF, XDS, Scalepack, or SHELX format.
    """
    return self._file_type

  def set_file_type(self, file_type):
    self._file_type = file_type

  @property
  def file_object(self):
    """Synonym for file_content()"""
    return self._file_object

  @property
  def file_content(self):
    """Return the underlying format-specific object containing file data."""
    return self._file_object

  @property
  def file_server(self):
    """
    For reflection files only, returns an
    :py:class:`iotbx.reflection_file_utils.reflection_file_server` object
    containing the extracted Miller arrays.  Note that this will implicitly
    merge any non-unique observations.
    """
    from iotbx.reflection_file_utils import reflection_file_server
    if (self._file_server is None):
      if (self._file_type == "hkl"):
        self._file_server = reflection_file_server(
          crystal_symmetry=None,
          force_symmetry=True,
          reflection_files=[self._file_object],
          err=sys.stderr)
    return self._file_server

  def _try_as_pdb(self):
    """
    PDB parser, actually tries both 'classic' PDB and mmCIF formats.
    """
    import iotbx.pdb.hierarchy
    try :
      pdb_inp = iotbx.pdb.hierarchy.input(self.file_name)
    except ValueError as e :
      raise Sorry(str(e))
    if (pdb_inp.hierarchy.models_size() == 0):
      raise ValueError("No ATOM or HETATM records found in %s."%self.file_name)
    self._file_type = "pdb"
    self._file_object = pdb_inp

  def _try_as_hkl(self):
    from iotbx.reflection_file_reader import any_reflection_file
    # XXX this is unfortunate, but unicode breaks Boost.Python extensions
    hkl_file = any_reflection_file(str(self.file_name))
    assert (hkl_file.file_type() is not None), "Not a valid reflections file."
    self._file_type = "hkl"
    self._file_object = hkl_file

  def _try_as_cif(self):
    # XXX hack to avoid choking on CCP4 maps and images
    file_ext = os.path.splitext(self.file_name)[1]
    assert (not file_ext in [".ccp4", ".img", ".osc", ".mccd"])
    import iotbx.cif
    from iotbx.reflection_file_reader import any_reflection_file
    cif_file = any_reflection_file(str(self.file_name))
    if cif_file.file_type() is not None:
      self._file_object = cif_file
      self._file_type = "hkl"
    else:
      # Try to read as simple cif model file.  If it fails, use the
      #  input_hierarchy_pair reader as previously (totally unknown
      #  function or reason for this reader. See:
      #    modules/cctbx_project/iotbx/pdb/hierarchy.py
      try :
        pdb_inp = iotbx.pdb.hierarchy.input(self.file_name)
        if len(pdb_inp.input.atoms()) > 0:
          self._file_object = pdb_inp
          self._file_type = "pdb"
      except Exception as e :
        pass
      if not self._file_object:
        from iotbx.pdb.mmcif import cif_input
        from iotbx.pdb.hierarchy import input_hierarchy_pair
        try:
          cif_in = cif_input(file_name=self.file_name)
          pdb_inp = input_hierarchy_pair(cif_in, cif_in.hierarchy)
          self._file_object  = pdb_inp
          self._file_type = "pdb"
        except Exception as e:
          if (str(e).startswith("Space group is incompatible") or
              str(e).startswith("The space group") ):
            raise
          else:
            pdb_inp = iotbx.cif.reader(file_path=self.file_name,
              strict=False)
            self._file_type = "cif"
            self._file_object = pdb_inp

  def _try_as_phil(self):
    from iotbx.phil import parse as parse_phil
    phil_object = parse_phil(file_name=self.file_name, process_includes=True)
    assert (len(phil_object.objects) > 0), "Empty parameter file."
    self._file_type = "phil"
    self._file_object = phil_object

  def _try_as_seq(self):
    # XXX hack to avoid choking on CCP4 maps
    assert (not self.file_name.endswith(".ccp4"))
    # XXX hack to avoid choking on NCS files:
    assert (not self.file_name.endswith(".ncs"))
    assert (not self.file_name.endswith(".ncs_spec"))

    from iotbx.bioinformatics import any_sequence_format
    objects, non_compliant = any_sequence_format(self.file_name)
    assert (objects is not None), "No sequence data found in file."
    assert (len(non_compliant) == 0), "Misformatted data in file."
    for seq_obj in objects :
      assert (not "-" in seq_obj.sequence)
    self._file_object = objects
#    self._try_as_txt()
#    assert len(self._file_object) != 0
#    for _line in self._file_object.splitlines():
#      assert not _line.startswith(" ")
#      line = re.sub(" ", "", _line)
#      assert ((len(line) == 0) or
#              (line[0] == ">") or
#              (line == "*") or
#              ((line[-1] == '*') and line[:-1].isalpha()) or
#              line.isalpha())
    # Filter out empty sequences
    tmp = []
    for o in self._file_object:
      if(len(o)!=0): tmp.append(o)
    self._file_object = tmp
    self._file_type = "seq"

  def _try_as_hhr(self):
    from iotbx.bioinformatics import any_hh_file
    hh_object = any_hh_file(self.file_name)
    assert (not hh_object.query in ["", None])
    self._file_object = hh_object
    self._file_type = "hhr"

  def _try_as_a3m(self):
    from iotbx.bioinformatics import any_a3m_file
    a3m_object = any_a3m_file(self.file_name)
    self._file_object = a3m_object
    self._file_type = "a3m"

  def _try_as_aln(self):
    from iotbx.bioinformatics import any_alignment_file
    aln_object = any_alignment_file(self.file_name)
    self._file_object = aln_object
    self._file_type = "aln"

  def _try_as_xplor_map(self):
    import iotbx.xplor.map
    map_object = iotbx.xplor.map.reader(file_name=str(self.file_name))
    self._file_type = "xplor_map"
    self._file_object = map_object

  def _try_as_ccp4_map(self):
    from iotbx.map_manager import map_manager
    from libtbx.utils import null_out
    map_object=map_manager(file_name=str(self.file_name),log=null_out())
    self._file_type = "ccp4_map"
    self._file_object = map_object

  def _try_as_pkl(self):
    with open(self.file_name, "rb") as fh:
      pkl_object = pickle.load(fh)
    self._file_type = "pkl"
    self._file_object = pkl_object

  def _try_as_txt(self):
    with open(self.file_name) as fh:
      file_as_string = fh.read()
    file_as_ascii = to_str(file_as_string)
    self._file_type = "txt"
    self._file_object = file_as_string

  def _try_as_xml(self):
    import xml.dom.minidom
    xml_in = xml.dom.minidom.parse(self.file_name)
    self._file_type = "xml"
    self._file_object = xml_in

  def _try_as_img(self):
    from iotbx.detectors import ImageFactory
    img = ImageFactory(self.file_name)
    img.read()
    self._file_type = "img"
    self._file_object = img

  def _try_as_ncs(self):
    from mmtbx.ncs.ncs import ncs
    from libtbx.utils import null_out
    ncs_object=ncs()
    try: # see if we can read biomtr records
      pdb_inp=iotbx.pdb.input(file_name=self.file_name)
      ncs_object.ncs_from_pdb_input_BIOMT(pdb_inp=pdb_inp,log=null_out(),
        quiet=True)
    except Exception as e: # try as regular ncs object
      ncs_object.read_ncs(file_name=self.file_name,log=sys.stdout,quiet=True)
    assert ncs_object.max_operators() > 0
    self._file_object = ncs_object
    self._file_type = "ncs"

  def try_all_types(self):
    for filetype in self.valid_types :
      if (filetype in self._tried_types) : continue
      read_method = getattr(self, "_try_as_%s" % filetype)
      try :
        read_method()
      except KeyboardInterrupt :
        raise
      except Exception as e :
        if (str(e).startswith("Space group is incompatible") or
            str(e).startswith("The space group") ):
          raise
        self._errors[filetype] = str(e)
        self._file_type = None
        self._file_object = None
        continue
      else :
        if self._file_type is not None :
          break

  def crystal_symmetry(self):
    """
    Extract the crystal symmetry (if any).  Only valid for model (PDB/mmCIF)
    and reflection files.
    """
    from cctbx import crystal
    if(self._file_type == "pdb"):
      return self._file_object.input.crystal_symmetry()
    elif(self._file_type == "hkl"):
      try:
        return self._file_object.file_content().crystal_symmetry()
      except AttributeError:
        return self._file_object.as_miller_arrays()[0].crystal_symmetry()
    elif(self._file_type == "ccp4_map"):
      return self._file_object.crystal_symmetry()
    else:
      raise NotImplementedError()

  def file_info(self, show_file_size=True):
    """
    Format a string containing the file type and size.
    """
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
    if self._file_type == None :
      return "Unknown file%s" % file_size_str
    else :
      return "%s%s" % (self.__descriptions__[self._file_type],
        file_size_str)

  def assert_file_type(self, expected_type):
    """
    Verify that the automatically determined file type is the expected format.
    """
    if (expected_type is None):
      return self
    elif (self._file_type == expected_type):
      return self
    else :
      raise Sorry(("Expected file type '%s' for %s, got '%s'.  This is " +
        "almost certainly a bug; please contact the developers.") %
        (expected_type, str(self.file_name), str(self._file_type)))

  def check_file_type(self, expected_type=None, multiple_formats=()):
    """
    Verify that the automatically determined file type is the expected format,
    with the option to consider multiple formats.
    """
    if (expected_type is not None):
      if (self._file_type != expected_type):
        raise Sorry(("This file format ('%s') is not supported as input for "+
          "this field; only files of type '%s' are allowed.") % (
          standard_file_descriptions.get(self._file_type, "Unknown"),
          standard_file_descriptions.get(expected_type, "Unknown")))
    else :
      assert (len(multiple_formats) > 0)
      if (not self._file_type in multiple_formats):
        raise Sorry(("This file format ('%s') is not supported as input for "+
          "this field; only the following types are supported:\n  %s") % (
          standard_file_descriptions.get(self._file_type, "Unknown"),
          "\n  ".join([ standard_file_descriptions.get(f, "Unknown")
                        for f in multiple_formats ])))
    return self

  def show_summary(self, out=sys.stdout):
    """
    Print out some basic information about the file.
    """
    if (self._file_type is None):
      print("File type could not be determined.", file=out)
    else :
      print("File name: %s" % self.file_name, file=out)
      print("Format: %s (%s)" % (self._file_type,
        standard_file_descriptions.get(self._file_type, "unknown")), file=out)
    if (self._file_type == "pdb"):
      print("Atoms in file: %d" % (len(self._file_object.input.atoms())), file=out)
      title = "\n".join(self._file_object.input.title_section())
      if (title != ""):
        print("Title section:", file=out)
        print(title, file=out)
    elif (self._file_type == "hkl"):
      for array in self.file_server.miller_arrays :
        print("", file=out)
        array.show_comprehensive_summary(f=out)
    elif (self._file_type == "ccp4_map"):
      self._file_object.show_summary(out)

def any_file_fast(file_name,
              get_processed_file=False,
              valid_types=supported_file_types,
              allow_directories=False,
              force_type=None,
              input_class=None):
  """
  mimics any_file, but without parsing - will instead guess the file type from
  the extension.  for most output files produced by cctbx/phenix this is
  relatively safe; for files of unknown provenance it is less effective.
  """
  assert (not get_processed_file) and (force_type is None)
  if not os.path.exists(file_name):
    raise Sorry("Couldn't find the file %s" % file_name)
  elif os.path.isdir(file_name):
    if not allow_directories :
      raise Sorry("This application does not support folders as input.")
    else :
      return directory_input(file_name)
  elif not os.path.isfile(file_name):
    raise Sorry("%s is not a valid file.")
  else :
    if input_class is None :
      input_class = any_file_fast_input
    return input_class(file_name=file_name,
      valid_types=valid_types)

class any_file_fast_input(object):
  __extensions__ = standard_file_extensions
  __descriptions__ = standard_file_descriptions
  def __init__(self, file_name, valid_types):
    self.valid_types = valid_types
    self.file_name = file_name
    self.file_object = None
    self.file_type = None
    self.file_server = None
    self.file_description = None
    file_base, file_ext, compress_ext = splitext(file_name)
    for file_type in valid_types :
      if file_ext[1:] in self.__extensions__[file_type] :
        self.file_type = file_type
    if self.file_type is not None :
      self.file_description = self.__descriptions__[self.file_type]
    # XXX this is a huge hole in this method - for reflection files, the
    # PHENIX GUI displays the specific format, which requires actually reading
    # in the file.  the stub class below will work for obvious formats.
    if (self.file_type == "hkl"):
      class fake_data_object(object):
        def __init__(self, ext):
          self.ext = ext
        def file_type(self):
          if (self.ext == ".mtz") : return "CCP4 MTZ"
          elif (self.ext == ".sca") : return "Scalepack"
          else : return "Data (other)"
      self.file_object = fake_data_object(file_ext)

class directory_input(object):
  def __init__(self, dir_name):
    self.file_name = dir_name
    self.file_object = dircache.listdir(dir_name)
    self.file_server = None
    self.file_type = "dir"
    self.file_size = os.path.getsize(dir_name)

  def file_info(self, show_file_size=False):
    return "Folder"

class group_files(object):
  def __init__(self,
                file_names,
                template_format="pdb",
                group_by_directory=True):
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
      if (file_type == template_format):
        base, ext, compress_ext = splitext(file_name)
        templates.append(base)
        template_dirs.append(os.path.dirname(file_name))
        self.grouped_files.append([file_name])
      else :
        other_files.append(file_name)
    if (len(templates) == 0):
      raise Sorry("Can't find any files of the expected format ('%s')." %
        template_format)
    if (len(set(templates)) != len(templates)):
      raise Sorry("Multiple files with identical root names.")
    for file_name in other_files :
      group_name = find_closest_base_name(
        file_name=file_name,
        base_name=splitext(file_name)[0],
        templates=templates)
      if (group_name == ""):
        self.ambiguous_files.append(file_name)
      elif (group_name is not None):
        i = templates.index(group_name)
        self.grouped_files[i].append(file_name)
      else :
        if group_by_directory :
          dir_name = os.path.dirname(file_name)
          group_name = find_closest_base_name(
            file_name=dir_name,
            base_name=dir_name,
            templates=template_dirs)
          if (group_name == ""):
            self.ambiguous_files.append(file_name)
          elif (group_name is not None):
            i = template_dirs.index(group_name)
            self.grouped_files[i].append(file_name)
          else :
            self.ungrouped_files.append(file_name)
        else :
          self.ungrouped_files.append(file_name)

def find_closest_base_name(file_name, base_name, templates):
  groups = []
  for base in templates :
    if file_name.startswith(base) or base.startswith(base_name):
      groups.append(base)
  if (len(groups) == 1):
#    print file_name, groups[0]
    return groups[0]
  elif (len(groups) > 1):
#    print file_name, groups
    prefix_len = [ os.path.commonprefix([g, file_name]) for g in groups ]
    max_common_prefix = max(prefix_len)
    if (prefix_len.count(max_common_prefix) > 1):
      return ""
    else :
      return groups[ prefix_len.index(max_common_prefix) ]
  return None

#---end
