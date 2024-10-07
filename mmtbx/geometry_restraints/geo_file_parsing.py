from __future__ import division
from collections import defaultdict
from itertools import chain
import shlex
import re
from libtbx.utils import Sorry
from cctbx import geometry_restraints
from cctbx.geometry_restraints.linking_class import linking_class
origin_ids = linking_class()

from libtbx import group_args

class Entry:
  """
  Base class for an 'entry' in a geo file. Analogous to both a proxy and restraint
    Subclass Entry for each type of entry that may be encountered during geo file parsing.

  Roles:
    1. Store restraint-specific parameters for parsing
    2. Implement interpretation of raw .geo file lines to common attributes
      Accessing the 'entry.record' attribute will return a dict of all available data

    3.
        Atom labels can appear in two lists:
        a. self.atom_labels: The raw atom label strings,  often atom.id_str()
        b. self.i_seqs: The integer i_seqs coming from:
                          a. The .geo file if all atom labels are also ints
                          b. A model file provided upon initialization and id_str labels

    4. Implement a method to convert entry object to cctbx proxy object
  """
  name = None                 # Name or label for the class

  def __init__(self,lines,line_idx=None,origin_id=0,origin_label='covalent'):
    """
    An entry is initialized first, then data is added with entry.lines.append()

    Attributes:
      origin_id (int): the integer origin id for an entry
    """
    # Parsing data structures
    self.lines = lines               # raw .geo lines
    # self.line_idxs =[line_idx]
    self.i_seqs = []                 # list of integer i_seqs (if possible)
    self.atom_labels  = []           # list of string atom label from .geo
    self._numerical = None           # a dict of numerical geo data
    self.labels_are_i_seqs = None    # boolean, atom labels are i_seqs
    self.labels_are_id_strs = None   # boolean, atom labels are id_strs
    self.origin_id = origin_id
    # self.origin_label = origin_label

    # Initialize result data structures
    self._proxy = None
    self._record = None

    self._prepare()

    # Check if labels are i_seqs (integers)
    self.labels_are_i_seqs, self.i_seqs = self._check_labels_are_i_seqs(self.atom_labels)

  def _prepare(self):
    # Parse Atom labels
    values = []
    for line in self.lines[:-2]:
      if not line.startswith(" "):
        line = line.replace(line.split()[0],"") # remove name like 'bond', 'angle'
      #line = line.replace("pdb=","")
      value = line.strip()
      values.append(value)

    self.atom_labels = values

    # Numerical labels
    labels = self.lines[-2]
    numerical_labels = labels.split()

    # Numerical values
    values = self.lines[-1]
    numerical_values =  values.split()
    numerical_values = [self._coerce_type(v) for v in numerical_values]
    self._numerical = dict(zip(numerical_labels,numerical_values))

  def _check_labels_are_i_seqs(self,atom_labels):
    """
    If all the labels are integers, assume i_seqs
    """
    i_seqs = [self._try_int(v) for v in atom_labels]
    check =   all([isinstance(i_seq,int)  for i_seq in i_seqs])
    if not check:
      i_seqs = []
    return check, i_seqs

  @property
  def has_i_seqs(self):
    return len(self.i_seqs)>0

  @property
  def record(self):
    """
    A dictionary representation of an entry
      Atom labels are split up into single atom key:value pairs
    """
    if not self._record:
      d = {
        "i_seqs":self.i_seqs,
        "atom_labels":self.atom_labels,
        }
      d.update(self._numerical)
      d["origin_id"] = self.origin_id
      # d["origin_label"] = self.origin_label
      self._record = d
    return self._record



  @property
  def ideal(self):
    """
    Provide the restraint ideal value
    """
    return float(self._numerical["ideal"])

  @property
  def weight(self):
    """
    Provide the restraint weight value
    """
    return float(self._numerical["weight"])


  @property
  def proxy(self):
    """
    Only create a proxy object if necessary, and if so only do it once
    """
    if not self._proxy and self.has_i_seqs:
      self._proxy = self.to_proxy()
    return self._proxy


  def _try_int(self, val):
    try:
      return int(val)
    except (ValueError, TypeError):
      return None

  def _try_float(self, val):
    try:
      return float(val)
    except (ValueError, TypeError):
      return None

  def _try_numeric(self,val):
    out = self._try_int(val)
    if out:
      return out
    out = self._try_float(val)
    if out:
      return out
    return None

  def _coerce_type(self,val):
    out = self._try_numeric(val)
    if out:
      return out
    return val # input (probably string)

### Start of Entry subclasses


class NonBondedEntry(Entry):
  name = "nonbonded"

class AngleEntry(Entry):
  name = "angle"

  def to_proxy(self):
    proxy = geometry_restraints.angle_proxy(
    i_seqs=self.i_seqs,
    angle_ideal=self.ideal,
    weight=self.weight,
    origin_id=self.origin_id)
    return proxy


class BondEntry(Entry):
  name = "bond"

  def to_proxy(self):
    proxy = geometry_restraints.bond_simple_proxy(
            i_seqs=self.i_seqs,
            distance_ideal=self.ideal,
            weight=self.weight,
            origin_id=self.origin_id,
            )
    return proxy

class DihedralEntry(Entry):
  name = "dihedral"
  @property
  def is_harmonic(self):
    return "harmonic" in self._numerical.keys()

  @property
  def is_sinusoidal(self):
    return "sinusoidal" in self._numerical.keys()

  @property
  def periodicity(self):
    if self.is_harmonic:
      return int(self._numerical["harmonic"])
    else:
      return int(self._numerical["sinusoidal"])

  def to_proxy(self):
    proxy = geometry_restraints.dihedral_proxy(
      i_seqs=self.i_seqs,
      angle_ideal=self.ideal,
      weight=self.weight,
      periodicity=self.periodicity,
      alt_angle_ideals=None,
      origin_id=self.origin_id)
    return proxy


class ChiralityEntry(Entry):
  name = "chirality"

  @property
  def both_signs(self):
    return bool(self._numerical["both_signs"])

  def to_proxy(self):
    proxy = geometry_restraints.chirality_proxy(
      i_seqs=self.i_seqs,
      volume_ideal=self.ideal,
      weight=self.weight,
      both_signs=self.both_signs,
      origin_id=self.origin_id
      )
    return proxy

class PlaneEntry(Entry):
  """
  Planes are very different because they can have varying number of atoms
    Most methods must be overridden
  """
  name = "plane"


  def _prepare(self):
    """
    Interpret lines from a Plane entry.
    """

    self.atom_labels = []
    nums = [[None]*5 for l in range(len(self.lines)-1)] # 5 values
    for i,line in enumerate(self.lines[1:]):
      line = line.replace(self.name,"")
      pdb_part = re.search(r'pdb="([^"]*)"', line)
      if pdb_part:
        pdb_value = pdb_part.group(0)  # Preserve the whole pdb="..." string
        remaining_line = line.replace(pdb_value, "")

        parts = shlex.split(remaining_line)
        parts.insert(0, pdb_value)
        # check seg id
        if len(parts)>1:
          if "segid=" in parts[1]:
            parts = [" ".join(parts[0:2])] + parts[2:]

      else:
        # No pdb="..." part
        parts = shlex.split(line)
      comp_value = parts[0]

      self.atom_labels.append(comp_value.strip())
      for j,p in enumerate(parts[1:]):
        nums[i][j] = p

    line = self.lines[0]
    numerical_labels = line.strip().split()

    # fill empty values down columns
    for i,row in enumerate(nums):
      if i>0:
        nums[i][-2] = nums[0][-2]
        nums[i][-1] = nums[0][-1]
    nums_T = list(map(list, zip(*nums))) # transpose
    numerical_values = nums_T

    numerical_values = [[self._coerce_type(v) for v in vals] for vals in numerical_values]
    self._numerical = dict(zip(numerical_labels,numerical_values))

  @property
  def n_atoms(self):
    # Override for planes
    return len(self.atom_labels)


  @property
  def weights(self):
    return [float(w) for w in self._numerical["weight"]]

  def to_proxy(self):
    proxy = geometry_restraints.planarity_proxy(
        i_seqs=self.i_seqs,
        weights=self.weights,
        origin_id=self.origin_id,
    )
    return proxy

class ParallelityEntry(Entry):
  name = "parallelity"

  def __init__(self,*args,**kwargs):
    self._atom_labels = []

    # add extra fields for 'j' atoms
    self.j_seqs = []
    self.atom_labels_i = []
    self.atom_labels_j = []
    super().__init__(*args,**kwargs)
    self.labels_are_i_seqs, self.i_seqs = self._check_labels_are_i_seqs(self.atom_labels)
    self.i_seqs, self.j_seqs = self.i_seqs # unpack tuple


  @property
  def atom_labels(self):
    return self.atom_labels_i + self.atom_labels_j

  @atom_labels.setter
  def atom_labels(self,value):
    self._atom_labels = value

  def _prepare(self):
    """
    Interpret lines from a Parallelity entry.
    """
    line0 = self.lines[0]
    line1 = self.lines[1]
    plane_2_idx = line0.index("plane 2")

    all_parts = []
    for line in self.lines[1:]:
      pdb_parts = re.findall(r'pdb="([^"]*)"', line)

      if len(pdb_parts) >= 1:
        pdb_value_i = f'pdb="{pdb_parts[0]}"'

        remaining_line = line.replace(pdb_value_i, "", 1)

        if len(pdb_parts) == 2:
          pdb_value_j = f'pdb="{pdb_parts[1]}"'

          remaining_line = remaining_line.replace(pdb_value_j, "", 1)
          parts = shlex.split(remaining_line)

          parts.insert(0, pdb_value_i)
          parts.insert(1, pdb_value_j)
        else:
          # Only one pdb="..." found, handle just pdb_value_i
          parts = shlex.split(remaining_line)
          parts.insert(0, pdb_value_i)

      else:
        # No pdb="..." part found, just split the line
        parts = shlex.split(line)
      all_parts.append(parts)

    for i,row in enumerate(all_parts):
      val = row[0].replace(self.name,"").strip()
      self.atom_labels_i.append(val)

    for j,row in enumerate(all_parts):
      if len(row)>1:
        val = row[1].replace(self.name,"").strip()
        self.atom_labels_j.append(val)

    num_idx = plane_2_idx+len("plane 2")
    num_labels = shlex.split(line0[num_idx:])
    num_values = all_parts[0][2:]

    num_values = [self._coerce_type(v) for v in num_values]
    self._numerical = dict(zip(num_labels,num_values))

  def _check_labels_are_i_seqs(self,*args):
    """
    If all the labels are integers, assume i_seqs
    """
    check_func = super()._check_labels_are_i_seqs
    check_i, i_seqs = check_func(self.atom_labels_i)
    check_j, j_seqs = check_func(self.atom_labels_j)
    check = check_i and check_j
    return check, (i_seqs, j_seqs)

    if not check:
      i_seqs = []
    return check, i_seqs

  @property
  def n_atoms(self):
    return len(self.atom_labels_i) + len(sel.atom_labels_j)

  def to_proxy(self):
    proxy = geometry_restraints.parallelity_proxy(
        i_seqs=self.i_seqs,
        j_seqs=self.j_seqs,
        weight=10, # Can get from .geo?
        origin_id=self.origin_id,
    )
    return proxy

  @property
  def record(self):
    """
    A dictionary representation of an entry
      Atom labels are split up into single atom key:value pairs
    """
    d = {
      "i_seqs":self.i_seqs,
      "j_seqs":self.j_seqs,
      "atom_labels_i":self.atom_labels_i,
      "atom_labels_j":self.atom_labels_j,
      }
    d.update(self._numerical)
    d["origin_id"] = self.origin_id
    return d

### End Entry classes

entry_class_config_default = (
  # (Entry subclass, Entry title, entry trigger)
  #
  # Entry title: A value returned by origin_ids.get_origin_label_and_internal, determines Entry subclass
  # Entry trigger: a text field in a .geo file which indicates the start of a new entry section
  (NonBondedEntry,"Nonbonded",'nonbonded'),
  (AngleEntry,"Bond angle","angle"),
  (BondEntry,"Bond","bond"),
  (DihedralEntry, "Dihedral angle",'dihedral'),
  (ChiralityEntry,"Chirality",'chirality'),
  (ParallelityEntry,"Parallelity", "plane 1"),
  (PlaneEntry,"Planarity", "delta"),
  (PlaneEntry,"Plane","delta"),
  )


class GeoParser:
  """
  A container class to hold parsed geometry entries, and to implement
    the parsing functions. The full functionality is to go:
    .geo text file ==> cctbx proxy objects

  Usage:
    parser = GeoParser(geo_lines)     # Initialize
    entries = container.entries       # Access data as Entry instances
    records = container.records       # Access data as lists of dicts


    proxies = container.proxies       # Access data as lists of proxies, if possible
  """


  def __init__(self,geo_lines,model=None,entry_class_config=None):
    """
     Initialize with a list of Entry subclasses
    """

    # Set initial arguments
    if not entry_class_config:
      entry_class_config = entry_class_config_default
    self.entry_class_config = entry_class_config
    self.entry_trigger_dict = {value[2]:value[0] for value in self.entry_class_config}
    self.entry_class_trigger_dict = {value[1]:value[0] for value in self.entry_class_config}
    self.model = model

    # Initialize parsing variables
    self.lines = geo_lines + ["\n"]
    self.current_entry= None
    self.current_entry_class = None
    self.current_origin_id = 0
    self.current_origin_label = 'covalent'

    # Initialize result variables
    self.entries = defaultdict(list) # Entry instances
    self._proxies = None             # cctbx geometry proxies
    self._records = None             # dictionaries

    # Parse the file
    self._parse()

    # If model present, add i_seqs
    if self.model:
      self._fill_labels_from_model(self.model)

    # Debug
    print(self.entries_list[-1].record)

  @property
  def proxies(self):
    return self._proxies

  @property
  def has_proxies(self):
    return self.proxies is not None

  def _fill_labels_from_model(self, model):
    """
    Add i_seq attributes for each atom on each entry
    """
    if not self.labels_are_i_seqs:
      # make i_seq:id_str mapping
      map_iseq_to_idstr = {
        atom.i_seq:atom.id_str()
        for atom in model.get_atoms()}

      # Reverse
      map_idstr_to_iseq = {id_str:i_seq for i_seq,id_str in map_iseq_to_idstr.items()}

      for entry in self.entries_list:
        # Case where entry was loaded with id_strs, can load i_seqs
        entry.labels_are_id_strs = True
        for val in entry.atom_labels:
          if val not in map_idstr_to_iseq:
            entry.labels_are_id_strs = False

        if entry.labels_are_id_strs:
          entry.i_seqs = [map_idstr_to_iseq[idstr] for idstr in entry.atom_labels]

  @property
  def entries_list(self):
    """
    Coalesce all entries (of all restraint entry types) into ao single list
    """
    return list(chain.from_iterable(self.entries.values()))

  @property
  def records_list(self):
    """
    Coalesce all records into ao single list
    """
    return list(chain.from_iterable(self.records.values()))

  @property
  def proxies_list(self):
    """
    Coalesce all proxies into ao single list
    """
    if self.proxies is not None:
      return list(chain.from_iterable(self.proxies.values()))

  def _parse_geo_file_header_full(self, line):
    """
    TODO: Move to linking class

    Collect three items:
      1. Entry type: ('Bond', 'Angle', etc... The class of restraint)
      2. Origin ID: (Integer id for Secondary Structure, Metal coordination, etc)
      3. entry_type_start_word: ('bond', 'angle', etc... indication of a new entry)

    The primary functionality here is 'origin_ids.get_origin_label_and_internal()',
    but it returns 'covalent' for all unrecognized inputs.
    So in order to 'trust' the result as a real switch back to covalent, we need to know
    which lines are headers. This is done using indentifying strings in entry_class_trigger_dict
    """
    class_trigger_to_entry = {value[1]:value[2] for value in self.entry_class_config}
    answer = None
    entry_class_trigger  = self._startswith_plural(line.split("|")[0],self.entry_class_trigger_dict.keys())
    if entry_class_trigger:
      origin_id, _ = origin_ids.get_origin_label_and_internal(line)
      restraint_type = self.entry_class_trigger_dict[entry_class_trigger]
      answer = group_args(
          entry_type = restraint_type,
          origin_id = origin_id,
          entry_type_start_word = class_trigger_to_entry[entry_class_trigger])
    return answer

  def _end_entry(self, entry_start_line_number, i, entries_info):
    # check that this is not the first line of the block
    if entry_start_line_number != -1:
      # Create an entry
      lines_for_entry = self.lines[entry_start_line_number:i]
      new_entry = entries_info.entry_type(
          lines=lines_for_entry,
          origin_id=entries_info.origin_id)
      # add new_entry to somewhere
      self.entries[entries_info.entry_type.name].append(new_entry)
    return i

  def _parse(self):
    entries_info=None # entries_info = group_args(origin_id, entry_type, entry_type_start_word)
    entry_start_line_number = -1
    for i, l in enumerate(self.lines + ["\n"]):
      if entries_info is None:
        # new block starts because we don't know entries_info
        entries_info = self._parse_geo_file_header_full(l) # "get_origin_label_and_internal"
        entry_start_line_number = -1
      elif l.startswith("Sorted by"):
        pass
      elif l.strip() == "":
        # block ends
        entry_start_line_number = self._end_entry(entry_start_line_number, i, entries_info)
        entries_info=None
      elif l.strip().startswith(entries_info.entry_type_start_word):
        # entry ends, the next is coming
        entry_start_line_number = self._end_entry(entry_start_line_number, i, entries_info)
      else:
        # entry continues, do nothing
        pass

  def _startswith_plural(self,text, labels, strip=False):
    """
    Utility function to compare a list of labels with a text string.
      Returns the first matched label.
    """
    if strip:
      text = text.strip()
    for label in labels:
      if text.startswith(label):
        return label
    return None

  @property
  def records(self):
    """
    Express the restraint entries as a dict of lists of dicts (records).
      This ONLY uses the entry object, not the proxies/grm.
      Meaning the only data present will be what was parsed.
    """
    if not self._records:
      record_dict = {}
      for entry_name, entries in self.entries.items():
        if len(entries)>0:
          record_dict[entry_name] = [entry.record for entry in entries]
      self._records =  record_dict
    return self._records

  @property
  def labels_are_i_seqs(self):
    return all([entry.labels_are_i_seqs for entry in self.entries_list])

  def build_proxies(self):
    """
    Convert the entry objects to cctbx proxy objects.
      Collect into a dict of lists of proxies.
    """
    if not self.model and not self.labels_are_i_seqs:
      raise Sorry("Cannot build proxies without instantiating with a model.")
    self._proxies = defaultdict(list)
    for entries in self.entries.values():
      if len(entries)>0:
        entry_class = entries[0].__class__
      if hasattr(entry_class,"to_proxy") and not hasattr(entry_class,"ignore"):
        self._proxies[entry_class.name] = [entry.proxy for entry in entries]
    return self.proxies
