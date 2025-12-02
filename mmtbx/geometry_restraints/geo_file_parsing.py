from __future__ import division
from collections import defaultdict
from itertools import chain
import shlex
import re
from libtbx import group_args
from libtbx.utils import Sorry
from cctbx import geometry_restraints
from cctbx.geometry_restraints.linking_class import linking_class
from cctbx.array_family import flex
from cctbx import crystal
from cctbx.crystal import direct_space_asu

origin_ids = linking_class()


class Entry:
  """
  1. Base class for an 'entry' in a geo file. Analogous to both a proxy and restraint
    Entry is subclasses for each type of entry that may be encountered during geo file parsing.
    The functions for parsing and value type conversion are implemented on the subclasses

  2. Accessing the 'entry.record' attribute will return a plain dict with all relevant data

  3. Atom labels can appear in two lists:
    a. self.atom_labels: The raw atom label strings,  often atom.id_str()
    b. self.i_seqs: The integer i_seqs coming from:
                      a. The .geo file if all atom labels are also ints
                      b. A model file provided upon initialization and id_str labels

  4. Implement a method to convert entry object to cctbx proxy object (self.to_proxy())
  """

  def __init__(self,lines,origin_id=0,origin_label="covalent"):
    """
    An entry is initialized first, then data is added with entry.lines.append()

    Attributes:
      origin_id (int): the integer origin id for an entry
    """
    # Parsing data structures
    self.lines = lines               # raw .geo lines
    self.i_seqs = []                 # list of integer i_seqs (if possible)
    self.sites_cart = None           # Flex vec3 array, coord for each i_seq if available
    self.atom_labels  = []           # list of string atom label from .geo
    self._numerical = None           # a dict of numerical geo data
    self.origin_id = origin_id
    self.origin_label = origin_label

    # Initialize result data structures
    self._proxy = None
    self._record = None

    self._prepare()

    # Check if labels are i_seqs (integers)
    labels_are_i_seqs, self.i_seqs = self._check_labels_are_i_seqs(self.atom_labels)
    if labels_are_i_seqs:
      self.atom_labels = []

  def labels_are_available(self):
    return len(self.atom_labels)>0

  def i_seqs_are_available(self):
    return len(self.i_seqs)>0

  def _prepare(self):
    # Parse Atom labels
    values = []
    for line in self.lines[:-2]:
      if not line.startswith(" "):
        line = line.replace(line.split()[0],"") # remove name like 'bond', 'angle'
      value = line.strip()
      values.append(value)

    self.atom_labels = values

    # Numerical labels
    labels = self.lines[-2]
    numerical_labels = labels.split()

    # Numerical values
    values = self.lines[-1]
    numerical_values =  values.split()
    numerical_values = [coerce_type(v) for v in numerical_values]
    self._numerical = dict(zip(numerical_labels,numerical_values))

  def _check_labels_are_i_seqs(self,atom_labels):
    """
    If all the labels are integers, assume i_seqs
    """
    i_seqs = [try_int(v) for v in atom_labels]
    check =   all([isinstance(i_seq,int)  for i_seq in i_seqs])
    if not check:
      i_seqs = []
    return check, i_seqs

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
    if not self._proxy and self.i_seqs_are_available():
      self._proxy = self.to_proxy()
    return self._proxy

### Start of Entry subclasses

class NonBondedEntry(Entry):
  pass

class AngleEntry(Entry):

  def to_proxy(self):
    proxy = geometry_restraints.angle_proxy(
    i_seqs=self.i_seqs,
    angle_ideal=self.ideal,
    weight=self.weight,
    origin_id=self.origin_id)
    return proxy


class BondEntry(Entry):

  def has_sym_op(self):
    if "sym.op." in self._numerical.keys() and self._numerical["sym.op."] not in ["",None]:
      return True
    else:
      return False

  def to_proxy(self):
    if self.has_sym_op():
      asu_mappings = direct_space_asu.non_crystallographic_asu_mappings(
      sites_cart=self.sites_cart)
      pair_generator = crystal.neighbors_fast_pair_generator(
      asu_mappings=asu_mappings,
      distance_cutoff=5)
      pair = geometry_restraints.bond_asu_proxy(
        pair=next(pair_generator),
        distance_ideal=self.ideal,
        weight=self.weight,
        origin_id=self.origin_id)
      proxy = geometry_restraints.bond_asu_proxy(pair=pair, params=pair)
      return proxy

    else:
      proxy = geometry_restraints.bond_simple_proxy(
              i_seqs=self.i_seqs,
              distance_ideal=self.ideal,
              weight=self.weight,
              origin_id=self.origin_id,
              )
    return proxy

class DihedralEntry(Entry):

  def is_harmonic(self):
    return "harmonic" in self._numerical.keys()

  def is_sinusoidal(self):
    return "sinusoidal" in self._numerical.keys()

  @property
  def periodicity(self):
    if self.is_harmonic():
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
      line = line.replace("plane","")
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

    numerical_values = [[coerce_type(v) for v in vals] for vals in numerical_values]
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

  def __init__(self,*args,**kwargs):
    """
    Parallelity is special because there are two sets of atoms for each plane i, j
    Here, self.i_seqs is only the first plane, not all atoms as in other entries.
    """
    self._atom_labels = []

    # add extra fields for 'j' atoms
    self.j_seqs = []
    self.atom_labels_i = []
    self.atom_labels_j = []
    super().__init__(*args,**kwargs)
    labels_are_i_seqs, self.i_seqs = self._check_labels_are_i_seqs(self.atom_labels)
    if labels_are_i_seqs:
      self.atom_labels_i = []
      self.atom_labels_j = []
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
        pdb_value_i = 'pdb="{}"'.format(pdb_parts[0])

        remaining_line = line.replace(pdb_value_i, "", 1)

        if len(pdb_parts) == 2:
          pdb_value_j = 'pdb="{}"'.format(pdb_parts[1])

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
      val = row[0].replace("parallelity","").strip()
      self.atom_labels_i.append(val)

    for j,row in enumerate(all_parts):
      if len(row)>1:
        val = row[1].replace("parallelity","").strip()
        self.atom_labels_j.append(val)

    num_idx = plane_2_idx+len("plane 2")
    num_labels = shlex.split(line0[num_idx:])
    num_values = all_parts[0][2:]

    num_values = [coerce_type(v) for v in num_values]
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

entry_config_default = {
  # Keys: The name of the restraint class in the geo header
  #
  # entry_class: An Entry subclass, defines parsing function and handles value type conversion
  # entry_trigger: a unique text field in a .geo file which indicates the start of a new entry
  "Nonbonded": {
    "entry_class": NonBondedEntry,
    "entry_trigger": "nonbonded"
  },
  "Bond angle": {
    "entry_class": AngleEntry,
    "entry_trigger": "angle"
  },
  "Bond": {
    "entry_class": BondEntry,
    "entry_trigger": "bond"
  },
  "Dihedral angle": {
    "entry_class": DihedralEntry,
    "entry_trigger": "dihedral"
  },
  "Chirality": {
    "entry_class": ChiralityEntry,
    "entry_trigger": "chirality"
  },
  "Parallelity": {
    "entry_class": ParallelityEntry,
    "entry_trigger": "plane 1"
  },
  "Planarity": {
    "entry_class": PlaneEntry,
    "entry_trigger": "delta"
  }
}


class GeoParser:
  """
  A container class to hold parsed geometry entries, and to implement
    the parsing functions. The full functionality is to go:
    .geo text file ==> cctbx proxy objects

  Usage:
    parser = GeoParser(geo_lines)     # Initialize
    entries = container.entries       # Access data as Entry instances
    records = container.records       # Access data as plain dictionaries


    proxies = container.proxies       # Access data as proxy objects, if possible
  """


  def __init__(self,geo_lines,model=None,entry_config=None):
    """
     Initialize with a list of Entry subclasses

     Params:
       geo_lines (list):             List of line strings from .geo file
       model (mmtbx.model.manager):  Optional model file, a source of atom_label: i_seq matching.
       entry_config (dict):          Configuration dict
    """

    # Set initial arguments
    self.entry_config = (entry_config if entry_config else entry_config_default)
    self.model = model

    # Initialize parsing variables
    self.lines = geo_lines + ["\n"]

    # Initialize result variables
    self.entries = defaultdict(list) # Entry instances
    self._proxies = None             # cctbx geometry proxies
    self._records = None             # dictionaries

    # Parse the file
    self._parse()

    # If model present, add i_seqs
    if self.model:
      self._fill_labels_from_model(self.model)
      self._fill_sites_cart_from_model(self.model)

  @property
  def proxies(self):
    return self._proxies

  def has_proxies(self):
    return self.proxies is not None

  def i_seqs_are_available(self):
    return all([entry.i_seqs_are_available() for entry in self.entries_list])

  def _fill_sites_cart_from_model(self,model):
    # Verify iseqs
    if not self.i_seqs_are_available():
      self._fill_labels_from_model(model)
    # get sites cart
    sites_cart = model.get_sites_cart()
    for entry in self.entries_list:
      entry.sites_cart = sites_cart.select(flex.size_t(entry.i_seqs))


  def _fill_labels_from_model(self, model):
    """
    Add i_seq attributes for each atom on each entry
    """
    if not self.i_seqs_are_available():
      # make i_seq:id_str mapping
      map_idstr_to_iseq = {
        atom.id_str():atom.i_seq for atom in model.get_atoms()}

      for entry in self.entries_list:
        # Case where entry was loaded with id_strs, can load i_seqs
        labels_are_id_strs = True
        for val in entry.atom_labels:
          if val not in map_idstr_to_iseq:
            labels_are_id_strs = False

        if labels_are_id_strs:
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

  def _end_entry(self, entry_start_line_number, i, entries_info):
    # check that this is not the first line of the block
    if entry_start_line_number != -1:
      # Create an entry
      lines_for_entry = self.lines[entry_start_line_number:i]
      new_entry = entries_info.entry_class(
          lines=lines_for_entry,
          origin_id=entries_info.origin_id,
          origin_label=entries_info.origin_label,
          )
      # add new_entry to somewhere
      self.entries[entries_info.origin_header].append(new_entry)
    return i

  def _parse(self):
    entries_info=None
    entry_start_line_number = -1
    for i, l in enumerate(self.lines):
      if entries_info is None:
        # new block starts because we don't know entries_info

        # Query linking class with line
        result = origin_ids.get_origin_label_and_internal(l)
        if result:
          # if recognized as header, unpack result, store in entries_info
          origin_id, header, label, num = result
          entry_class = self.entry_config[header]["entry_class"]
          entry_trigger = self.entry_config[header]["entry_trigger"]
          entries_info = group_args(
              entry_class = entry_class,
              origin_id = origin_id,
              origin_label = label,
              origin_header = header,
              entry_trigger = entry_trigger,
              )
          entry_start_line_number = -1

      elif l.startswith("Sorted by"):
        pass
      elif l.strip() == "":
        # block ends
        entry_start_line_number = self._end_entry(entry_start_line_number, i, entries_info)
        entries_info=None
      elif l.strip().startswith(entries_info.entry_trigger):
        # entry ends, the next is coming
        entry_start_line_number = self._end_entry(entry_start_line_number, i, entries_info)
      else:
        # entry continues, do nothing
        pass

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
      self._records = record_dict
    return self._records


  def build_proxies(self):
    """
    Convert the entry objects to cctbx proxy objects.
      Collect into a dict of lists of proxies.
    """
    if not self.model and not self.i_seqs_are_available():
      raise Sorry("Cannot build proxies without instantiating with a model.")
    self._proxies = defaultdict(list)
    for entries in self.entries.values():
      if len(entries)>0:
        entry_class = entries[0].__class__
      if hasattr(entry_class,"to_proxy") and not hasattr(entry_class,"ignore"):
        self._proxies[entry_class.__name__] = [entry.proxy for entry in entries]
    return self.proxies


def try_int(val):
  """
  Try to convert to int
  """
  try:
    return int(val)
  except (ValueError, TypeError):
    return None

def try_float(val):
  """
  Try to convert to flaot
  """
  try:
    return float(val)
  except (ValueError, TypeError):
    return None

def try_numeric(val):
  """
  Try to convert input to 1) int, 2) float
  Otherwise return None
  """
  out = try_int(val)
  if out is not None:
    return out
  out = try_float(val)
  if out is not None:
    return out
  return None

def coerce_type(val):
  """
  If input can be numeric, convert it
  Else return unmodified
  """
  out = try_numeric(val)
  if out is not None:
    return out
  return val # input (probably string)

if __name__ == '__main__':
  import sys
  print(sys.argv)
  f=open(sys.argv[1])
  geo_lines=f.readlines()
  del f
  parser = GeoParser(geo_lines)     # Initialize
  print(dir(parser))
  entries = parser.entries_list       # Access data as Entry instances
  records = parser.records_list
  proxies = parser.proxies
  for record in records:
    print(record)
    # print(dir(record))
    # print(type(record))
    break
    if record['origin_id']!=0:
      break
  for entry in entries:
    print(entry)
    print(dir(entry))
    print(type(entry))
    print(entry.ideal)
    assert 0

