from collections import defaultdict
from pathlib import Path
from itertools import chain
import shlex
import re
from libtbx import group_args
from libtbx.utils import Sorry
from cctbx import geometry_restraints
from cctbx.geometry_restraints import bond, angle, dihedral, chirality, parallelity, planarity
from cctbx.geometry_restraints.linking_class import linking_class
origin_ids = linking_class()

def build_origin_maps():
  """
  Collect origin id data and return two dicts. 

  
  1. For a given restraint type, get list of possible origin headers

      origin_map1 = 
        {restraint_label: concatenated_header_list
        }
      
  2. For a (restraint,header) pair, get the origin label
  
      origin_map2 = 
        {(restraint_label, header): origin_label}

  """
  origin_map1 = {}
  origin_map2 = {}
  internals = ["bonds", "angles", "torsions", "planes", "chirals", "parallelities"]
  
  for origin_label, origins in origin_ids.data.items():
    for i in origins.internals:
      if len(origins) >= 4:
        header = origins[3]
        if header:
          restraint_label = internals[i]
          valid_headers = [v for v in header if v]
          
          if restraint_label in origin_map1:
            origin_map1[restraint_label].extend(valid_headers)
          else:
            origin_map1[restraint_label] = valid_headers.copy()
          
          for v in valid_headers:
            origin_map2[(restraint_label, v)] = origin_label
  
  for restraint_label, header_list in origin_map1.items():
    origin_map1[restraint_label] = list(set(header_list))
    
  return origin_map1, origin_map2

  
class Entry:
  """
  Base class for an 'entry' in a geo file. Analogous to both a proxy and restraint
    Subclass Entry for each type of entry that may be encountered during geo file parsing.

  Roles:
    1. Store restraint-specific parameters for parsing
    2. Implement interpretation method of raw .geo file lines to common attributes
      Accessing the 'entry.record' attribute will return a dict of all available data 

    3.
        Atom labels can appear in two ways:
        a. atom_1: The raw atom label read in, often atom.id_str()
        b. i_seq_1: The integer i_seq

    4. Implement a method to convert entry object to cctbx proxy object
  """

  # Implement these attriutes in subbclasses
  name = None                 # Name or label for the class
  entry_class_trigger = None  # Text trigger to identify entry class
  entry_trigger = None        # Text trigger to identify new entry instance
  n_atoms = None              # explicit number of atoms to expect
  n_values = None             # explicit number of values to expect
  default_origin_id = 0       # default origin_id integer for an entry
  
  # Only build origin maps once, in class definition
  origin_map1, origin_map2 = build_origin_maps()

  def __init__(self,lines,origin_id=None):
    """
    An entry is initialized first, then data is added with entry.lines.append()

    Attributes:
      origin_id (int): the integer origin id for an entry
    """
    # Parsing data structures
    self.lines = lines               # raw .geo lines
    self._i_seqs = None              # a dict of atom i_seqs
    self._compositional= None        # a dict of unspecified atom labels
    self._numerical = None           # a dict of numerical geo data
    self._labels_are_i_seqs = None   # boolean, atom labels are i_seqs
    self._labels_are_id_strs = None  # boolean, atom labels are id_strs

    # Set an origin id (if not provided)
    if not origin_id:                 
      origin_id = self.default_origin_id
    self.origin_id = origin_id
    self.origin_label = origin_ids.get_origin_key(origin_id)
    

    # Initialize result data structures
    self._proxy = None  
    self._record = None


    # parse lines
    self._prepare()

  @classmethod
  @property
  def origin_headers(cls):
    assert cls.internals, "Must define the 'internals' class variable"
    return cls.origin_map1[cls.internals]
    

  def _prepare(self):

    # Compositional labels
    stem = "atom"
    comp_labels =  [f"{stem}_{i+1}" for i in range(len(self.lines[:-2]))]
    if self.n_atoms:
      # test labels
      n_atoms =  len(comp_labels)
      assert n_atoms == self.n_atoms, (comp_labels,n_atoms,self.n_atoms)
    comp_labels

    # Compositional values
    values = []
    for ln,line in self.lines[:-2]:
      if not line.startswith(" "):
        line = line.replace(line.split()[0],"") # remove name like 'bond', 'angle'
      line = line.replace("pdb=","")
      values.append(line)
      
    comp_values = [self._remove_outer_quotes(value.strip()) for value in values]
    if self.n_atoms:
      n_atoms =  len(comp_values)
      assert n_atoms == self.n_atoms, (comp_values,n_atoms,self.n_atoms)

    comp_values
    self._compositional = dict(zip(comp_labels,comp_values))

    # Numerical labels
    ln, labels = self.lines[-2]
    numerical_labels = labels.split()
    if self.n_values:
      # test labels
      n_values =  len(numerical_labels)
      assert n_values == self.n_values or n_values == self.n_values-1, (numerical_labels,n_values,self.n_values)
    
    # Numerical values
    ln, values = self.lines[-1]
    numerical_values =  values.split()
    if self.n_values:
      n_values =  len(numerical_values)
      assert n_values == self.n_values or n_values == self.n_values-1, (numerical_values,n_values,self.n_values)
    numerical_values = [self._coerce_type(v) for v in numerical_values]
    self._numerical = dict(zip(numerical_labels,numerical_values))
    

  @property
  def record(self):
    """
    A dictionary representation of an entry
      Atom labels are split up into single atom key:value pairs
    """
    d = {}
    if self._i_seqs:
      d.update(self._i_seqs)
    d.update(self._compositional)
    d.update(self._numerical)
    d["origin_id"] = self.origin_id
    return d

  @property
  def labels_are_i_seqs(self):
    """
    If all the labels are integers, assume i_seqs
    """
    if not self._labels_are_i_seqs:
      if all([isinstance(self._try_int(v),int)  for v in self._compositional.values()]):
        self._labels_are_i_seqs = True
        self._form_i_seqs(self._compositional)
      else:
        self._labels_are_i_seqs = False
    return self._labels_are_i_seqs

  @labels_are_i_seqs.setter
  def labels_are_i_seqs(self,value):
    self._labels_are_i_seqs = value
    if self._labels_are_i_seqs and self._labels_are_id_strs:
      raise Sorry("Atom labels cannot be both i_seqs and id_strs")

  @property
  def labels_are_id_strs(self):
    return self._labels_are_id_strs

  @labels_are_id_strs.setter
  def labels_are_id_strs(self,value):
    self._labels_are_id_strs = value
    if self._labels_are_i_seqs and self._labels_are_id_strs:
      raise Sorry("Atom labels cannot be both i_seqs and id_strs")

  @property
  def i_seqs(self):
    """
    Provide the i_seqs as a list
    """
    if self._i_seqs:
      return list(self._i_seqs.values())
  
  def _form_i_seqs(self,compositional_dict):
    self._i_seqs =  {f"i_seq_{i+1}":int(v) for i,v in enumerate(compositional_dict.values())}
  
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

  @staticmethod
  def _remove_outer_quotes(s):
    """
    Utility function to clean up quotes when present in .geo text
    """
    if (s.startswith('"') and s.endswith('"')) or (s.startswith("'") and s.endswith("'")):
      return s[1:-1]
    return s
        
  @property
  def proxy(self):
    """
    Only create a proxy object if necessary, and if so only do it once
    """
    if not self._proxy:
      self._proxy = self.to_proxy()
    return self._proxy


  def _try_int(self,val):
    try:
      return int(val)
    except:
      return None

  def _try_float(self,val):
    try:
      return float(val)
    except:
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
  entry_class_trigger = "Nonbonded interactions"
  entry_trigger = "nonbonded"
  n_atoms = 2
  n_values = 3
  internals = "bonds"

  
class AngleEntry(Entry):
  name = "angle"
  entry_class_trigger = "Bond angle restraints"
  entry_trigger = "angle"
  n_atoms = 3
  n_values = 6
  internals = "angles"

    
  def to_proxy(self):
    proxy = geometry_restraints.angle_proxy(
    i_seqs=self.i_seqs,
    angle_ideal=self.ideal,
    weight=self.weight,
    origin_id=self.origin_id)
    return proxy


class BondEntry(Entry):
  name = "bond"
  entry_class_trigger = "Bond restraints"
  entry_trigger = "bond"
  n_atoms = 2
  n_values = 6
  internals = "bonds"
  
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
  entry_class_trigger = "Dihedral angle restraints"
  entry_trigger = "dihedral"
  n_atoms = 4
  n_values = 7
  internals = "torsions"

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
  name = "chiral"
  entry_class_trigger= "Chirality restraints"
  entry_trigger = "chirality"
  n_atoms = 4
  n_values = 7
  internals = "chirals"
    
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
  entry_class_trigger = "Planarity restraints"
  entry_trigger = "plane"
  n_values = 5
  internals = "planes"



  def _prepare(self):
    """
    Interpret lines from a Plane entry.
    """
    compositional_labels = []
    compositional_values = []

    nums = [[None]*self.n_values for l in range(len(self.lines)-1)]
    for i,(ln,line) in enumerate(self.lines[1:]):
      line = line.replace(self.name,"")
      parts = shlex.split(line)
      comp_value = parts[0]
      stem = "atom"
      
      comp_label = f"{stem}_{i}"
      comp_value = comp_value.replace("pdb=","")
      comp_value = self._remove_outer_quotes(comp_value)
      compositional_labels.append(comp_label)
      compositional_values.append(comp_value)
      for j,p in enumerate(parts[1:]):
        nums[i][j] = p

    ln,line = self.lines[0]
    numerical_labels = line.strip().split()
    
    # fill empty values down columns
    for i,row in enumerate(nums):
      if i>0:
        nums[i][-2] = nums[0][-2]
        nums[i][-1] = nums[0][-1]
    nums_T = list(map(list, zip(*nums))) # transpose
    numerical_values = nums_T
    
    numerical_values = [[self._coerce_type(v) for v in vals] for vals in numerical_values]
    self._compositional = dict(zip(compositional_labels,compositional_values))
    self._numerical = dict(zip(numerical_labels,numerical_values))

  @property
  def n_atoms(self):
    # Override for planes
    return len(self._compositional)


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
  entry_trigger = "plane 1"
  n_atoms = None
  n_values = None
  internals = "parallelities"

  def _prepare(self):
    """
    Interpret lines from a Parallelity entry. 
    """
    ln0,line0 = self.lines[0]
    ln1,line1 = self.lines[1]
    plane_2_idx = line0.index("plane 2")

    comp_labels = []
    comp_values = []
    all_parts = []
    for ln,line in self.lines[1:]:
      #parts = [line[:plane_2_idx]]+shlex.split(line[plane_2_idx:])
      parts  = shlex.split(line)
      all_parts.append(parts)

    # Determine label type
    test_label = all_parts[0][0]
    stem = "atom"

    for i,row in enumerate(all_parts):
      comp_labels.append(f"{stem}_i_{i}")
      val = self._remove_outer_quotes(row[0].replace(self.name,"").replace("pdb=",""))
      comp_values.append(val)
      
    for j,row in enumerate(all_parts):
      if len(row)>1:
        comp_labels.append(f"{stem}_j_{j}")
        val = self._remove_outer_quotes(row[1].replace(self.name,"").replace("pdb=",""))
        comp_values.append(val)

    num_idx = plane_2_idx+len("plane 2")
    num_labels = shlex.split(line0[num_idx:])
    num_values = all_parts[0][2:]

    self._compositional = dict(zip(comp_labels,comp_values))
    num_values = [self._coerce_type(v) for v in num_values]
    self._numerical = dict(zip(num_labels,num_values))

  @property
  def i_seqs(self):
    """
    Provide the i_seqs as a list
    """
    if self._i_seqs:
      return [v for k,v in self._i_seqs.items() if "_i_" in k]
  
  def _form_i_seqs(self,compositional_dict):
      self._i_seqs =  {}
      for k,v in compositional_dict.items():
        if "_i_" in k:
          suffix = "_i_"+k.split("_i_")[-1]
        elif "_j_" in k:
          suffix = "_j_"+k.split("_j_")[-1]

        self._i_seqs[f"i_seq{suffix}"] = int(v)

  @property
  def j_seqs(self):
    """
    Provide the i_seqs as a list
    """
    if self._i_seqs:
      return [int(v) for k,v in self._i_seqs.items() if "_j_" in k]

  @property
  def n_atoms(self):
    return len(self._compositional)
    
  def to_proxy(self):
    proxy = geometry_restraints.parallelity_proxy(
        i_seqs=self.i_seqs,
        j_seqs=self.j_seqs,
        weight=10, # Can get from .geo?
        origin_id=self.origin_id,
    )
    return proxy

class StackingParallelityEntry(ParallelityEntry):
  entry_class_trigger = "Stacking parallelity restraints"

class BasepairParallelityEntry(ParallelityEntry):
  entry_class_trigger = "Basepair parallelity restraints"

### End Entry classes


class GeoParseContainer:
  """
  A container class to hold parsed geometry entries, and to implement
    the parsing functions. The full functionality is to go:
    .geo text file ==> cctbx proxy objects

  Usage:
    container = GeoParseContainer()   # Initialize
    container.parse_file("model.geo") # Parse .geo file
    entries = container.entries       # Access data as Entry instances
    records = container.records       # Access data as lists of dicts

                                      # Associate a model
    container.model = model           # Not neccessary if i_seqs in .geo

    proxies = container.proxies       # Access data as lists of proxies
  """

  entry_classes_default = [
    NonBondedEntry,
    AngleEntry,
    BondEntry,
    DihedralEntry,
    ChiralityEntry,
    BasepairParallelityEntry,
    StackingParallelityEntry,
    PlaneEntry,
  ]

  @classmethod
  def from_geo_file(cls,geo_file,model=None,entry_classes=None):
    """
    Parse .geo file
    """
    with open(geo_file,"r") as f:
      return cls.from_geo_str(f.read(),model=model,entry_classes=entry_classes)

  @classmethod
  def from_geo_str(cls,geo_str,model=None,entry_classes=None):
    """
    Parse .geo file str
    """
    return cls(geo_str=geo_str,model=model,entry_classes=entry_classes)
    
  def __init__(self,geo_str=None,model=None,entry_classes=None):
    """
     Initialize with a list of Entry subclasses
    """
    # Set initial arguments
    if not entry_classes:
      entry_classes = self.entry_classes_default
    self.entry_classes = entry_classes
    self._model = model
    self.geo_str = geo_str
    assert isinstance(geo_str,str), "Must initialize with a .geo file string."

    # Initialize parsing variables
    self.lines = self.geo_str.split("\n")
    self.line_labels = []
    self.current_entry_lines = None
    self.current_entry_class = None
    self.current_origin_id = None

    # Names of each entry class ('bond','angle', etc)
    self.entry_names = [entry_class.name for entry_class in self.entry_classes]

    # Text triggers for start of new entry class in .geo file
    self.entry_class_trigger_dict = {
      entry_class.entry_class_trigger:entry_class for entry_class in self.entry_classes
    }
    # text triggers for start of new entry class instance in .geo file
    self.entry_triggers = [entry_class.entry_trigger for entry_class in self.entry_classes]
    


    # Initialize result variables
    self.entries = defaultdict(list) # Entry instances
    self._proxies = None             # cctbx geometry proxies
    self._records = None             # dictionaries

    # Parse the file
    self._parse()

    # If model present, add i_seqs
    if self.model:
      self._fill_labels_from_model(self.model)

    if self.model or self.labels_are_i_seqs:
      self._build_proxies()

    
  @property
  def model(self):
    return self._model

  @property
  def proxies(self):
    return self._proxies

  @property
  def has_proxies(self):
    return self.proxies is not None

  def _fill_labels_from_model(self,model):
    """
    Add i_seq attributes for each atom on each entry
    """
    if not self.labels_are_i_seqs:
      # make i_seq:id_str mapping
      map_iseq_to_idstr = {
        atom.i_seq:Entry._remove_outer_quotes(atom.id_str().replace("pdb=","")) 
        for atom in model.get_atoms()}

      # Reverse
      map_idstr_to_iseq = {id_str:i_seq for i_seq,id_str in map_iseq_to_idstr.items()}

      for entry in self.entries_list:
        # Case where entry was loaded with id_strs, can load i_seqs
        entry.labels_are_id_strs = True
        for val in entry._compositional.values():
          if val not in map_idstr_to_iseq:
            entry.labels_are_id_strs = False
          
        if entry.labels_are_id_strs:
          entry._form_i_seqs({k:map_idstr_to_iseq[idstr] for k,idstr in entry._compositional.items()})

  @property
  def entries_list(self):
    """
    Coalesce all entries (of all restraint entry types) into ao single list
    """
    return list(chain.from_iterable([entries for name,entries in self.entries.items()]))

  @property
  def records_list(self):
    """
    Coalesce all records into ao single list
    """
    return list(chain.from_iterable([records for name,records in self.records.items()]))

  @property
  def proxies_list(self):
    """
    Coalesce all proxies into ao single list
    """
    return list(chain.from_iterable([proxies for name,proxies in self.proxies.items()]))
  

  def _parse(self):
    
    last_line_label = "init"    # assign a label to every line
    for i,line in enumerate(self.lines):
      
      # Line labeling accounting
      assert last_line_label, f"Failed to set a line label at line: {i-1}"
      if i>0:
        self.line_labels.append(last_line_label)
      last_line_label = None

      # check blank
      line_strip = line.strip()
      if len(line_strip)==0:
        last_line_label = "blank"
        self._end_entry()
        continue

      # Check for start of new entry
      entry_trigger = self._startswith_plural(line, self.entry_triggers,strip=True)

      # Check for data
      if self.current_entry_lines and not entry_trigger:
        self.current_entry_lines.append((i,line))
        last_line_label = "data"
        continue
      

      if entry_trigger:
        # Start new entry
        self._end_entry()
        last_line_label = "entry_trigger"
        self.current_entry_lines = []
        if self.current_entry_class == PlaneEntry:
          j = i-1 # The prior line
          self.current_entry_lines.append((j,self.lines[j]))
        
        self.current_entry_lines.append((i,line))
        continue

      # Check for start of new entry section
      entry_class_trigger = self._startswith_plural(line, self.entry_class_trigger_dict.keys())
      if entry_class_trigger:
        # Start a new entry class
        last_line_label = "entry_class_trigger"
        self._end_entry()
        entry_class = self.entry_class_trigger_dict[entry_class_trigger]
        self.current_entry_class = entry_class
        self.current_origin_id = None
        continue

      # Check origin id trigger
      if self.current_entry_class:
        origin_headers = self.current_entry_class.origin_headers
        origin_trigger = self._startswith_plural(line, origin_headers)
        if origin_trigger:
          last_line_label = "origin_trigger"
          internals = self.current_entry_class.internals
          origin_label = self.current_entry_class.origin_map2[(internals,origin_trigger)]
          origin_id = origin_ids.get_origin_id(origin_label)
          self.current_origin_id = origin_id
          continue

      # Catch all
      last_line_label = "unknown"

    # End of function
    self.line_labels.append(last_line_label)


  def _end_entry(self):
    """
    Called at the end of an entry when parsing line-by-line
      Adds the entry to the container
    """
    if self.current_entry_lines:
      if self.current_entry_class == PlaneEntry:
        self.current_entry_lines = self.current_entry_lines[:-1]
      new_entry = self.current_entry_class(
        lines=self.current_entry_lines,
        origin_id=self.current_origin_id
      )
      self.entries[self.current_entry_class.name].append(new_entry)
      self.current_entry_lines = None

  def _startswith_plural(self,text, labels, strip=False):
    """
    Utility function to compare a list of labels with a text stringn.
      Returns the first matched label.
    """
    if strip:
      text = text.strip()
    for label in labels:
      if text.startswith(label):
        return label
    return None

  def _debug_print(self,start=0,stop=-1):
    if stop == -1:
      stop = len(self.lines)
    for i in range(start,stop):
      line = self.lines[i]
      label = self.line_labels[i]
      print(i,label,line)


  @property
  def records(self):
    """
    Express the restraint entries as a dict of lists of dicts (records).
      This ONLY uses the entry object, not the proxies/grm.
      Meaning the only data present will be what was parsed.
    """
    record_dict = {}
    for entry_name,entries in self.entries.items():
      if len(entries)>0:
        record_dict[entry_name] = [entry.record for entry in entries]
    return record_dict

  @property
  def labels_are_i_seqs(self):
    return all([entry.labels_are_i_seqs for entry in self.entries_list])

  def _build_proxies(self):
    """
    Convert the entry objects to cctbx proxy objects. 
      Collect into a dict of lists of proxies.
    """
    if not self.model and not self.labels_are_i_seqs: 
      raise Sorry("Cannot build proxies without instantiating with a model.")
    self._proxies = defaultdict(list)
    for entry_name,entries in self.entries.items():
      if len(entries)>0:
        entry_class = entries[0].__class__
      if hasattr(entry_class,"to_proxy") and not hasattr(entry_class,"ignore"):
        self._proxies[entry_class.name] = [entry.proxy for entry in entries]

