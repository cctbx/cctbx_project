
from __future__ import division
from libtbx import slots_getstate_setstate
from libtbx.str_utils import format_value
import sys

class entity (slots_getstate_setstate) :
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
  phenix_table_labels = []
  molprobity_table_labels = []

  def __init__ (self, **kwds) :
    for name in kwds.keys() :
      assert (name in self.__slots__), name
    for name in self.__slots__ :
      setattr(self, name, kwds.get(name, None))

  @staticmethod
  def header () :
    raise NotImplementedError()

  def format_score (self, replace_none_with="None") :
    return format_value(self.score_format, self.score,
      replace_none_with=replace_none_with)

  def is_outlier (self) :
    return self.outlier

  def as_string (self, prefix="") :
    raise NotImplementedError()

  def __str__ (self) :
    return self.as_string()

  def id_str (self, ignore_altloc=None) :
    """
    Returns a formatted (probably fixed-width) string describing the molecular
    entity being validation, independent of the analysis type.
    """
    raise NotImplementedError()

  def __hash__ (self) :
    return self.id_str().__hash__()

  def as_list (self) :
    """
    Optional; returns old format used by some tools in mmtbx.validation.
    """
    raise NotImplementedError()

  def as_table_row_molprobity (self) :
    """
    Returns a list of formatted table cells for display by MolProbity.
    """
    raise NotImplementedError()

  def as_table_row_phenix (self) :
    """
    Returns a list of formatted table cells for display by Phenix.
    """
    raise NotImplementedError()

  def as_kinemage (self) :
    """
    Returns a kinemage string for displaying an outlier.
    """
    raise NotImplementedError()

  def format_old (self) :
    raise NotImplementedError()

  def __eq__ (self, other) :
    """
    Compare two validation results to determine whether they correspond to the
    same molecular entity and analysis type.  This is intended to be used for
    analysis of a structure before-and-after refinement (etc.).
    """
    raise NotImplementedError()

  def __cmp__ (self, other) :
    return cmp(self.score, other.score)

  def is_single_residue_object (self) :
    raise NotImplementedError()

__residue_attr__ = [
  "chain_id",
  "resseq",
  "icode",
  "resname",
  "altloc",
  "segid",
]

class residue (entity) :
  """
  Base class for validation information about a single residue, which depending
  on context could mean either any one of the residue_group, atom_group, or
  residue objects from the PDB hierarchy.
  """
  __slots__ = entity.__slots__ + __residue_attr__ + ["occupancy"]

  def assert_all_attributes_defined (self) :
    for name in self.__residue_attr__ :
      assert (getattr(self, name) is not None) or (name == "segid")

  def id_str (self, ignore_altloc=False) :
    base = "%2s%4s%1s" % (self.chain_id, self.resseq, self.icode)
    if (not ignore_altloc) :
      base += "%1s" % self.altloc
    else :
      base += " "
    base += "%3s" % self.resname
    if (self.segid is not None) :
      base += " segid='%4s'" % self.segid
    return base

  def resseq_as_int (self) :
    from iotbx.pdb import hybrid_36
    return hybrid_36.hy36decode(len(self.resseq), self.resseq)

  def residue_id (self, ignore_altloc=False) :
    return self.id_str(ignore_altloc=ignore_altloc)

  def simple_id (self) :
    return ("%s%s%s" % (self.chain_id, self.resseq, self.icode)).strip()

  # XXX probably needs to be flexible about altloc...
  def is_same_residue (self, other, ignore_altloc=False) :
    if hasattr(other, "residue_id") :
      return (self.residue_id(ignore_altloc=ignore_altloc) ==
              other.residue_id(ignore_altloc=ignore_altloc))
    return (self.id_str(ignore_altloc=ignore_altloc) ==
            other.id_str(ignore_altloc=ignore_altloc))

  def is_same_residue_group (self, other) :
    return ((self.chain_id == other.chain_id) and
            (self.resseq == other.resseq) and
            (self.icode == other.icode) and
            (self.segid == other.segid))

  def residue_group_id_str (self) :
    return "%2s%4s%1s" % (self.chain_id, self.resseq, self.icode)

  def __eq__ (self, other) :
    assert type(self) == type(other), type(other)
    return self.id_str() == other.id_str()

  def is_single_residue_object (self) :
    return True

class atoms (entity) :
  """
  Base class for validation results involving a specific set of atoms, such
  as covalent geometry restraints, clashes, etc.
  """
  __atoms_attr__ = [
    "atoms_info",
  ]
  __slots__ = entity.__slots__ + __atoms_attr__

  def n_atoms (self) :
    return len(self.atoms_info)

  def __eq__ (self, other) :
    assert type(self) == type(other), type(other)
    return sorted(self.atoms_info) == sorted(other.atoms_info)

  def is_single_residue_object (self) :
    return False

  def get_altloc (self) :
    consensus_altloc = ''
    for atom in self.atoms_info :
      if (atom.altloc.strip() != '') :
        if (consensus_altloc == '') :
          consensus_altloc = atom.altloc
        else :
          assert (atom.altloc == consensus_altloc)
    return consensus_altloc

  def sites_cart (self) :
    return [ a.xyz for a in self.atoms_info ]

class atom (entity) :
  """
  Base class for validation results for a single atom.  This is distinct from
  the atom_info class above, which is used to track individual atoms within
  a multi-atom validation result.
  """
  __atom_attr__ = __residue_attr__ + [
    "name",
    "element",
    "b_iso",
    "occupancy",
    "i_seq", # XXX is this necessary?
    # XXX do we need u_aniso too?  anything else?
  ]
  __slots__ = entity.__slots__ + __atom_attr__

  def __init__ (self, **kwds) :
    pdb_atom = kwds.get("pdb_atom", None)
    if (pdb_atom is not None) :
      del kwds['pdb_atom']
    for name in kwds.keys() :
      assert (name in self.__slots__), name
    for name in self.__slots__ :
      setattr(self, name, kwds.get(name, None))
    if (pdb_atom is not None) :
      labels = pdb_atom.fetch_labels()
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

  def __eq__ (self, other) :
    assert type(self) == type(other), type(other)
    return self.id_str() == other.id_str()

  def residue_group_id_str (self) :
    return "%2s%4s%1s" % (self.chain_id, self.resseq, self.icode)

  def is_single_residue_object (self) :
    return True

#-----------------------------------------------------------------------
# Generic utility classes
class atom_info (slots_getstate_setstate) :
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
  ]
  __slots__ = __residue_attr__ + __atom_attr__

  def __init__ (self, **kwds) :
    for name in kwds.keys() :
      assert (name in self.__slots__), name
    for name in self.__slots__ :
      setattr(self, name, kwds.get(name, None))

  def id_str (self, ignore_altloc=True) :
    base = "%2s%4s%1s" % (self.chain_id, self.resseq, self.icode)
    if (not ignore_altloc) :
      base += "%1s" % self.altloc
    else :
      base += " "
    base += "%3s %4s" % (self.resname, self.name)
    if (self.segid is not None) :
      base += " segid='%4s'" % self.segid
    return base

  def __str__ (self) :
    raise NotImplementedError() # FIXME

  def __eq__ (self, other) :
    assert type(self) == type(other), type(other)
    return self.id_str() == other.id_str()

  def __cmp__ (self, other) :
    return cmp(self.id_str(), other.id_str())

  def residue_group_id_str (self) :
    return "%2s%4s%1s" % (self.chain_id, self.resseq, self.icode)

#-----------------------------------------------------------------------
# SPECIFIC RESULTS
# FIXME move these to specific modules

# rna_validate is a real mess, I've been meaning to work on this for a while
# maybe split up pucker and suite analysis into separate classes?
class rna_pucker (residue) : # TODO
  pass

class rna_suite (residue) : # TODO
  pass

#class rna_bond (bond) : # TODO
#  pass

#class rna_angle (angle) : # TODO
#  pass

# XXX others?
class water (atom) : # TODO
  pass

class residue_bfactor (residue) : # TODO
  pass

#-----------------------------------------------------------------------
class validation (slots_getstate_setstate) :
  __slots__ = ["n_outliers", "n_total", "results"]
  program_description = None
  output_header = None
  def __init__ (self) :
    self.n_outliers = 0
    self.n_total = 0
    self.results = []

  def get_outliers_count_and_fraction(self):
    if (self.n_total != 0):
      fraction = float(self.n_outliers) / self.n_total
      assert fraction <= 1.0
      return self.n_outliers, fraction
    return 0, 0.

  def get_outliers_goal (self) :
    raise NotImplementedError()

  def get_result_class (self) :
    raise NotImplementedError()

  def coot_todo (self) :
    return ""

  def show_old_output (self, out=sys.stdout, verbose=False) :
    if (verbose) :
      assert (self.output_header is not None)
      print >> out, self.output_header
    for result in self.results :
      print >> out, result.format_old()
    if (verbose) :
      self.show_summary(out)

  def show_summary (self, out=sys.stdout, prefix="") :
    raise NotImplementedError()

  def show (self, out=sys.stdout, prefix="  ", outliers_only=True,
      verbose=True) :
    if (len(self.results) > 0) :
      print >> out, prefix + self.get_result_class().header()
      for result in self.results :
        if (not outliers_only) or (result.is_outlier()) :
          print >> out, prefix + str(result)
    self.show_summary(out=out, prefix=prefix)

  def as_kinemage (self) :
    return None

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
