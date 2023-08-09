
from __future__ import absolute_import, division, print_function
from libtbx import slots_getstate_setstate
from libtbx.str_utils import format_value
import sys
import six
import json

class entity(slots_getstate_setstate):
  """
  Base class for all validation results.  This includes a boolean outlier flag,
  the information used to zoom in the Phenix GUI (optional, but strongly
  recommended), and some kind of numerical score (also optional, but strongly
  recommended - although some analyses may require multiple distinct scores).
  """
  __slots__ = [
    "atom_selection",
    "xyz",
    "outlier",
    "score",
  ]
  score_format = "%s"
  molprobity_table_labels = []

  def __init__(self, **kwds):
    for name in self.__slots__:
      setattr(self, name, None)
    for name, value in six.iteritems(kwds):
      assert (name in self.__slots__), name
      setattr(self, name, value)

  @staticmethod
  def header():
    """
    Format for header in result listings.
    """
    raise NotImplementedError()

  def format_score(self, replace_none_with="None"):
    return format_value(self.score_format, self.score,
      replace_none_with=replace_none_with)

  def is_outlier(self):
    return self.outlier

  def as_string(self, prefix=""):
    raise NotImplementedError()

  def __str__(self):
    return self.as_string()

  def id_str(self, ignore_altloc=None):
    """
    Returns a formatted (probably fixed-width) string describing the molecular
    entity being validation, independent of the analysis type.
    """
    raise NotImplementedError()

  def __hash__(self):
    return self.id_str().__hash__()

  def as_list(self):
    """
    Optional; returns old format used by some tools in mmtbx.validation.
    """
    raise NotImplementedError()

  def as_table_row_molprobity(self):
    """
    Returns a list of formatted table cells for display by MolProbity.
    """
    raise NotImplementedError()

  def as_table_row_phenix(self):
    """
    Returns a list of formatted table cells for display by Phenix.
    """
    raise NotImplementedError()

  def as_kinemage(self):
    """
    Returns a kinemage string for displaying an outlier.
    """
    raise NotImplementedError()

  def format_old(self):
    raise NotImplementedError()

  def as_JSON(self):
    """
    Returns a (empty) JSON object representing a single validation object.
    Should be overwritten by each validation script to actually output the data.
    """
    return json.dumps({})

  def __eq__(self, other):
    """
    Compare two validation results to determine whether they correspond to the
    same molecular entity and analysis type.  This is intended to be used for
    analysis of a structure before-and-after refinement (etc.).
    """
    return self.score == other.score

  def __cmp__(self, other):
    return cmp(self.score, other.score)

  def __ne__(self, other):
    return self.score != other.score

  def __lt__(self, other):
    return self.score < other.score

  def __le__(self, other):
    return self.score <= other.score

  def __gt__ (self, other):
    return self.score > other.score

  def __ge__(self, other):
    return self.score >= other.score

  def is_single_residue_object(self):
    raise NotImplementedError()

  def as_selection_string(self):
    """
    Returns PDB atom selection string for the atom(s) involved.
    """
    return None

  def zoom_info(self):
    """
    Returns data needed to zoom/recenter the graphics programs from the Phenix
    GUI.
    """
    return [ self.as_selection_string(), self.xyz ]

__residue_attr__ = [
  "chain_id",
  "resseq",
  "icode",
  "resname",
  "altloc",
  "segid",
]

class residue(entity):
  """
  Base class for validation information about a single residue, which depending
  on context could mean either any one of the residue_group, atom_group, or
  residue objects from the PDB hierarchy.
  """
  __slots__ = entity.__slots__ + __residue_attr__ + ["occupancy"]
  def _copy_constructor(self, other):
    for attr in __residue_attr__ :
      setattr(self, attr, getattr(other, attr))

  def assert_all_attributes_defined(self):
    for name in self.__slots__ :
      assert (getattr(self, name) is not None) or (name == "segid")

  def id_str(self, ignore_altloc=False):
    base = "%2s%4s%1s" % (self.chain_id, self.resseq, self.icode)
    if (not ignore_altloc):
      base += "%1s" % self.altloc
    else :
      base += " "
    base += "%3s" % self.resname
    if (self.segid is not None):
      base += " segid='%4s'" % self.segid
    return base

  def resseq_as_int(self):
    from iotbx.pdb import hybrid_36
    return hybrid_36.hy36decode(len(self.resseq), self.resseq)

  @property
  def resid(self):
    return "%4s%1s" % (self.resseq, self.icode)

  def residue_id(self, ignore_altloc=False):
    return self.id_str(ignore_altloc=ignore_altloc)

  def simple_id(self):
    return ("%s%s%s" % (self.chain_id, self.resseq, self.icode)).strip()

  # XXX probably needs to be flexible about altloc...
  def is_same_residue(self, other, ignore_altloc=False):
    if hasattr(other, "residue_id"):
      return (self.residue_id(ignore_altloc=ignore_altloc) ==
              other.residue_id(ignore_altloc=ignore_altloc))
    return (self.id_str(ignore_altloc=ignore_altloc) ==
            other.id_str(ignore_altloc=ignore_altloc))

  def is_same_residue_group(self, other):
    return ((self.chain_id == other.chain_id) and
            (self.resseq == other.resseq) and
            (self.icode == other.icode) and
            (self.segid == other.segid))

  def atom_group_id_str(self):
    return "%1s%3s%2s%4s%1s" % (self.altloc, self.resname, self.chain_id,
      self.resseq, self.icode)

  def residue_group_id_str(self):
    return "%2s%4s%1s" % (self.chain_id, self.resseq, self.icode)

  def __eq__(self, other):
    assert type(self) == type(other), type(other)
    return self.id_str() == other.id_str()

  def is_single_residue_object(self):
    return True

  def atom_selection_string(self):
    return "(chain '%s' and resid '%s' and resname %s and altloc '%s')" % \
      (self.chain_id, self.resid, self.resname, self.altloc)

  def set_coordinates_from_hierarchy(self, pdb_hierarchy,
      atom_selection_cache=None):
    if (atom_selection_cache is None):
      atom_selection_cache = pdb_hierarchy.atom_selection_cache()
    sel = atom_selection_cache.selection(self.atom_selection_string())
    assert (len(sel) > 0)
    self.xyz = pdb_hierarchy.atoms().select(sel).extract_xyz().mean()

  # for creating the hierarchical json output
  def nest_dict(self, level_list, upper_dict):
    inner_dict = {}
    if len(level_list) > 0:
      next_level = level_list[0]
      upper_dict[getattr(self, next_level)] = self.nest_dict(level_list[1:], inner_dict)
    else:
      return json.loads(self.as_JSON())
    return upper_dict

class atoms(entity):
  """
  Base class for validation results involving a specific set of atoms, such
  as covalent geometry restraints, clashes, etc.
  """
  __atoms_attr__ = [
    "atoms_info",
  ]
  __slots__ = entity.__slots__ + __atoms_attr__

  def n_atoms(self):
    return len(self.atoms_info)

  def __eq__(self, other):
    assert type(self) == type(other), type(other)
    return sorted(self.atoms_info) == sorted(other.atoms_info)

  def is_single_residue_object(self):
    return False

  def get_altloc(self):
    consensus_altloc = ''
    for atom in self.atoms_info :
      if (atom.altloc.strip() != ''):
        if (consensus_altloc == ''):
          consensus_altloc = atom.altloc
        else :
          assert (atom.altloc == consensus_altloc)
    return consensus_altloc

  def sites_cart(self):
    return [ a.xyz for a in self.atoms_info ]

  def is_in_chain(self, chain_id):
    if (chain_id == None):
      return True
    for a in self.atoms_info :
      if (a.chain_id == chain_id):
        return True
    return False

  def nest_dict(self, level_list, upper_dict):
    inner_dict = {}
    if len(level_list) > 0:
      next_level = level_list[0]
      try:
        upper_dict[getattr(self, next_level)] = self.nest_dict(level_list[1:], inner_dict)
      except AttributeError:
        upper_dict[getattr(self.atoms_info[0], next_level)] = self.nest_dict(level_list[1:], inner_dict)
    else:
      return [json.loads(self.as_JSON())]
    return upper_dict

  def merge_two_dicts(self, x, y):
    """Given two dictionaries, merge them into a new dict as a shallow copy, for json output."""
    z = x.copy()
    z.update(y)
    return z

class atom_base(slots_getstate_setstate):
  """
  Container for metadata for a single atom, in the context of validation
  results involving multiple atoms.  Intended to be used as-is inside various
  atoms classes.
  """
  __atom_attr__ = [
    "name",
    "element", # XXX is this necessary?
    "xyz",
    "symop",
    "occupancy",
    "b_iso",
    "model_id"
  ]
  __atom_slots__ = __residue_attr__ + __atom_attr__
  # XXX __slots__ should be left empty here

  def __init__(self, **kwds):
    pdb_atom = kwds.get("pdb_atom", None)
    if (pdb_atom is not None):
      del kwds['pdb_atom']
    for name in self.__slots__:
      setattr(self, name, None)
    for name, value in six.iteritems(kwds):
      assert (name in self.__slots__), name
      setattr(self, name, value)

    if (pdb_atom is not None):
      labels = pdb_atom.fetch_labels()
      self.model_id = labels.model_id
      self.chain_id = labels.chain_id
      self.resseq = labels.resseq
      self.icode = labels.icode
      self.resname = labels.resname
      self.altloc = labels.altloc
      self.segid = labels.segid
      self.name = pdb_atom.name
      self.xyz = pdb_atom.xyz
      self.occupancy = pdb_atom.occ
      self.b_iso = pdb_atom.b
      self.element = pdb_atom.element

  @property
  def resid(self):
    return "%4s%1s" % (self.resseq, self.icode)

  def __cmp__(self, other):
    return cmp(self.id_str(), other.id_str())

  def __eq__(self, other):
    assert isinstance(other, atom_base), type(other)
    return self.id_str() == other.id_str()

  def __ne__(self, other):
    return self.id_str() != other.id_str()

  def __lt__(self, other):
    return self.id_str() < other.id_str()

  def __le__(self, other):
    return self.id_str() <= other.id_str()

  def __gt__ (self, other):
    return self.id_str() > other.id_str()

  def __ge__(self, other):
    return self.id_str() >= other.id_str()

  def id_str(self, ignore_altloc=False, ignore_segid=False):
    base = "%2s%4s%1s" % (self.chain_id, self.resseq, self.icode)
    if (not ignore_altloc):
      base += "%1s" % self.altloc
    else :
      base += " "
    base += "%3s %4s" % (self.resname, self.name)
    if ( (self.segid is not None) and (not ignore_segid) ):
      base += " segid='%4s'" % self.segid
    return base

  def residue_group_id_str(self):
    return "%2s%4s%1s" % (self.chain_id, self.resseq, self.icode)

  def atom_group_id_str(self):
    return "%1s%3s%2s%4s%1s" % (self.altloc, self.resname, self.chain_id,
      self.resseq, self.icode)

class atom_info(atom_base):
  """
  Container for metadata for a single atom, in the context of validation
  results involving multiple atoms.  Intended to be used as-is inside various
  atoms classes.
  """
  __slots__ = atom_base.__atom_slots__ + ["symop"]

def get_atoms_info(pdb_atoms, iselection,
      use_segids_in_place_of_chainids=False):
  proxy_atoms = []
  for n, i_seq in enumerate(iselection):
    atom = pdb_atoms[i_seq]
    labels = atom.fetch_labels()
    if use_segids_in_place_of_chainids:
      chain_id = atom.segid
    else:
      chain_id = labels.chain_id
    info = atom_info(
      name=atom.name,
      element=atom.element,
      model_id=labels.model_id,
      chain_id=chain_id,
      resseq=labels.resseq,
      icode=labels.icode,
      resname=labels.resname,
      altloc=labels.altloc,
      occupancy=atom.occ,
      #segid=atom.segid,
      xyz=atom.xyz)
    proxy_atoms.append(info)
  return proxy_atoms

class atom(atom_base, entity):
  """
  Base class for validation results for a single atom.  This is distinct from
  the atom_info class above, which is used to track individual atoms within
  a multi-atom validation result.
  """
  __slots__ = entity.__slots__ + atom_base.__atom_slots__

  def is_single_residue_object(self):
    return True

#-----------------------------------------------------------------------
class validation(slots_getstate_setstate):
  """
  Container for a set of results from a single analysis (rotamers, clashes,
  etc.).  This is responsible for the console display of these results and
  associated statistics.  Individual modules will subclass this and override
  the unimplemented methods.
  """
  __slots__ = ["n_outliers", "n_total", "results", "_cache", "n_outliers_by_model", "n_total_by_model"]
  program_description = None
  output_header = None
  gui_list_headers = [] # for Phenix GUI ListCtrl widgets
  gui_formats = []      # for Phenix GUI ListCtrl widgets
  wx_column_widths = []

  def __init__(self):
    self.n_outliers = 0
    self.n_total = 0
    self.n_outliers_by_model = {}
    self.n_total_by_model = {}
    self.results = []
    self._cache = None
    assert (len(self.gui_list_headers) == len(self.gui_formats) ==
            len(self.wx_column_widths))

  def get_outliers_count_and_fraction(self):
    if (self.n_total != 0):
      fraction = float(self.n_outliers) / self.n_total
      assert fraction <= 1.0
      return self.n_outliers, fraction
    return 0, 0.

  def get_outliers_fraction_for_model(self, model_id):
    if (self.n_total_by_model[model_id] != 0):
      fraction = float(self.n_outliers_by_model[model_id]) / self.n_total_by_model[model_id]
      assert fraction <= 1.0
      return fraction
    return 0.

  @property
  def percent_outliers(self):
    n_outliers, frac_outliers = self.get_outliers_count_and_fraction()
    return frac_outliers * 100.

  def outlier_selection(self):
    """
    Return a flex.size_t object containing the i_seqs of atoms flagged as
    outliers (either individually or as part of an atom_group).  This needs
    to be implemented in the underlying classes unless they include a
    pre-built _outlier_i_seqs attribute.
    """
    if hasattr(self, "_outlier_i_seqs"):
      return self._outlier_i_seqs
    raise NotImplementedError()

  def get_outliers_goal(self):
    raise NotImplementedError()

  def get_result_class(self):
    raise NotImplementedError()

  def coot_todo(self):
    return ""

  def show_old_output(self, out=sys.stdout, verbose=False):
    """
    For backwards compatibility with output formats of older utilities
    (phenix.ramalyze et al.).
    """
    if (verbose):
      assert (self.output_header is not None)
      print(self.output_header, file=out)
    for result in self.results :
      print(result.format_old(), file=out)
    if (verbose):
      self.show_summary(out)

  def show_summary(self, out=sys.stdout, prefix=""):
    raise NotImplementedError()

  def show(self, out=sys.stdout, prefix="  ", outliers_only=True,
      verbose=True):
    if (len(self.results) > 0):
      print(prefix + self.get_result_class().header(), file=out)
      for result in self.iter_results(outliers_only):
        print(prefix + str(result), file=out)
    self.show_summary(out=out, prefix=prefix)

  def iter_results(self, outliers_only=True):
    for result in self.results :
      if (not outliers_only) or (result.is_outlier()):
        yield result

  def as_kinemage(self):
    return None

  def as_coot_data(self):
    """
    Return results in a format suitable for unpickling in Coot.
    """
    raise NotImplementedError()

  def as_gui_table_data(self, outliers_only=True, include_zoom=False):
    """
    Format results for display in the Phenix GUI.
    """
    table = []
    for result in self.iter_results(outliers_only):
      extra = []
      if (include_zoom):
        extra = result.zoom_info()
      row = result.as_table_row_phenix()
      assert (len(row) == len(self.gui_list_headers))
      table.append(row + extra)
    return table

  def save_table_data(self, file_name=None):
    """
    Save all results as a comma separated, text file
    """
    if (file_name is not None):
      outliers_only = False
      table = self.as_gui_table_data(outliers_only=outliers_only)
      if (len(table) > 0):
        f = open(file_name, 'w')
        f.write('%s\n' % ', '.join(self.gui_list_headers))
        for row in table:
          f.write('%s\n' % ', '.join([str(x) for x in row]))
        f.close()
        return True
    return False

  def merge_dict(self, a, b, path=None):
    """
    Recursive function for merging two dicts, merges b into a
    Mainly used to build hierarchical JSON outputs
    """
    if path is None: path = []
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                self.merge_dict(a[key], b[key], path + [str(key)])
            elif a[key] == b[key]:
                pass # same leaf value
            else:
                #raise Exception('Conflict at %s' % '.'.join(path + [str(key)]))
                # should only get called for JSONs that have a list
                # of different validations for a residue (e.g. clashes)
                a[key] = a[key]+b[key]
        else:
            a[key] = b[key]
    return a

  def find_residue(self, other=None, residue_id_str=None):
    assert ([other, residue_id_str].count(None) == 1)
    if (other is not None):
      if hasattr(other, "residue_group_id_str"):
        residue_id_str = other.residue_group_id_str()
      elif hasattr(other, "id_str"):
        residue_id_str = other.id_str()
      else :
        residue_id_str = str(other)
    if (self._cache is None):
      self._cache = {}
      for i_res, result in enumerate(self.results):
        result_id_str = result.residue_group_id_str()
        self._cache[result_id_str] = i_res
    i_res = self._cache.get(residue_id_str, None)
    if (i_res is not None):
      return self.results[i_res]
    return None

  def find_atom_group(self, other=None, atom_group_id_str=None):
    """
    Attempt to locate a result corresponding to a given atom_group object.
    """
    assert ([other, atom_group_id_str].count(None) == 1)
    if (other is not None):
      if hasattr(other, "atom_group_group_id_str"):
        atom_group_id_str = other.atom_group_group_id_str()
      elif hasattr(other, "id_str"):
        atom_group_id_str = other.id_str()
      else :
        atom_group_id_str = str(other)
    if (self._cache is None):
      self._cache = {}
      for i_res, result in enumerate(self.results):
        result_id_str = result.atom_group_group_id_str()
        self._cache[result_id_str] = i_res
    i_res = self._cache.get(atom_group_id_str, None)
    if (i_res is not None):
      return self.results[i_res]
    return None

class rna_geometry(validation):
  #__slots__ = validation.__slots__ + ["n_outliers_small_by_model", "n_outliers_large_by_model"]

  def show(self, out=sys.stdout, prefix="  ", verbose=True):
    if (len(self.results) > 0):
      print(prefix + self.get_result_class().header(), file=out)
      for result in self.results :
        print(result.as_string(prefix=prefix), file=out)
    self.show_summary(out=out, prefix=prefix)

class test_utils(object):

  def count_dict_values(prod, count_key, c=0):
    """
    for counting hierarchical values for testing hierarchical jsons
    """
    for mykey in prod:
      if prod[mykey] == count_key:
        c += 1
      if isinstance(prod[mykey], dict):
        # calls repeatedly
        c = test_utils.count_dict_values(prod[mykey], count_key, c)
      elif isinstance(prod[mykey], list):
        for d in prod[mykey]:
          if isinstance(d, dict):
            c = test_utils.count_dict_values(d, count_key, c)
    return c

class dummy_validation(object):
  """
  Placeholder for cases where values may be undefined because of molecule type
  (e.g. all-RNA structures) but we want to substitute None automatically.
  """
  def __getattr__(self, name):
    return None

  def __bool__(self):
    return False

  __nonzero__ = __bool__

molprobity_cmdline_phil_str = """
  model = None
      .type = path
      .help = "Model file (PDB or mmCIF)"

  outliers_only = False
    .type = bool
    .help = "Only display outliers"

  verbose = True
    .type = bool
    .help = '''Verbose'''
"""
