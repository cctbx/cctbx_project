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


  
class Entry:
  """
  Base class for an 'entry' in a geo file. Analogous to both a proxy and restraint
    Subclass Entry for each type of entry that may be encountered during geo file parsing.

  Roles:
    1. Store restraint-specific parameters for parsing
    2. Implement interpretation method of raw .geo file lines to common attributes
      Example attributes:
        compositional_labels: the atom label labels: (i_seq_1, id_str_2, etc)
        compositional_values: the atom label values: (10, ' O   ASN A   3 ')
        numerical_labels: 'ideal', 'delta', 'weight', 'residual', etc
        numerical_values: (The values associated with numerical labels)

      Accessing the 'entry.record' attribute will return a dict of all available data 
        Atom labels can appear in three ways:
        a. atom_1: The raw atom label read in
        b. i_seq_1: The integer i_seq
        c. id_str_1: The atom.id_str() value. Often written to .geo

    3. Implement a method to convert entry object to cctbx proxy object
  """

  # Implement these attriutes in subbclasses
  name = None                 # Name or label for the class
  entry_class_trigger = None  # Text trigger to identify entry class
  entry_trigger = None        # Text trigger to identify new entry instance
  n_atoms = None              # explicit number of atoms to expect
  n_values = None             # explicit number of values to expect
  origin_label_func = None    # string name of method on origin_ids object
  default_origin_id = 0       # default origin_id integer for an entry
  internals = None            # Internal filter for origin ids
  restraint_class = None      # the corresponding cctbx.geometry_restraints class
  

  def __init__(self,origin_id=None):
    """
    An entry is initialized first, then data is added with entry.lines.append()

    Attributes:
      origin_id (int): the integer origin id for an entry
    """
    self.lines = []
    self.proxy_seq = None
    self._proxy = None
    self._i_seqs = None
    self._id_strs = None
    self._compositional_dict = None
    self._numerical_dict = None
    self._labels_are_i_seqs = None
    self._labels_are_id_strs = None
    if not origin_id:
      origin_id = self.default_origin_id
    self.origin_id = origin_id

  @classmethod
  @property
  def origin_header_to_id(cls):
    """
    Create mapping between origin_header and origin id for a single entry type
    """
    if not hasattr(cls,"_origin_header_to_id"):
      cls._origin_header_to_id = {}
      if cls.origin_label_func:
        origin_label_func = getattr(origin_ids,cls.origin_label_func)
        origin_labels = origin_label_func()
        for origin_label in origin_labels:
          origin_id=origin_ids.get_origin_id(origin_label)
          origin_header=origin_ids.get_geo_file_header(origin_label, internals=cls.internals)
          cls._origin_header_to_id[origin_header] = origin_id
    return cls._origin_header_to_id

  @property
  def record(self):
    """
    A dictionary representation of an entry
      Atom labels are split up into single atom key:value pairs
    """
    d = {}
    if self._i_seqs:
      d.update(self._i_seqs)
    if self._id_strs:
       d.update(self._id_strs)
    d.update(self.compositional_dict)
    d.update(self.numerical_dict)
    return d

  @property
  def labels_are_i_seqs(self):
    if not self._labels_are_i_seqs:
      if all([self._is_int(v) for v in self.compositional_values]):
        self._labels_are_i_seqs = True
        self._form_i_seqs(self.compositional_dict)
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
    if value:
      self._form_id_strs(self.compositional_dict)
    if self._labels_are_i_seqs and self._labels_are_id_strs:
      raise Sorry("Atom labels cannot be both i_seqs and id_strs")

  @property
  def i_seqs(self):
    """
    Provide the i_seqs as a list
    """
    if self.labels_are_i_seqs and self._i_seqs:
      return list(self._i_seqs.values())
  
  def _form_i_seqs(self,compositional_dict):
    self._i_seqs =  {f"i_seq_{i+1}":int(v) for i,v in enumerate(compositional_dict.values())}
    
  def _form_id_strs(self,compositional_dict):
    self._id_strs =  {f"id_str_{i+1}":v for i,v in enumerate(compositional_dict.values())}

  @property
  def id_strs(self):
    """
    Provide the id_strs as a list
    """
    if self.labels_are_id_strs and self._id_strs:
      return list(self._id_strs.values())
  
  @property
  def ideal(self):
    """
    Provide the restraint ideal value
    """
    return float(self.numerical_dict["ideal"])

  @property
  def weight(self):
    """
    Provide the restraint weight value
    """
    return float(self.numerical_dict["weight"])

  @property
  def compositional_labels(self):
    """
    Generate atom labels titles from .geo text. 
      Assert consistent with n_atom expectation
    """
    stem = "atom"
    comp_labels =  [f"{stem}_{i+1}" for i in range(len(self.lines[:-2]))]
    if self.n_atoms:
      # test labels
      n_atoms =  len(comp_labels)
      assert n_atoms == self.n_atoms, (comp_labels,n_atoms,self.n_atoms)
    return comp_labels


  @property
  def compositional_values(self):
    """
    Generate the atom label values from .geo text
      Assert consistent with n_atom expectation
    """
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
    return comp_values

  @property
  def compositional_dict(self):
    """
    Provide a dictionary mapping compositional_labels: compositional_values
      So to get residual value, can do compositional_dict["id_str_1"]
    """
    if not self._compositional_dict:
      self._compositional_dict = dict(zip(self.compositional_labels,self.compositional_values))
    return self._compositional_dict

  @property
  def numerical_labels(self):
    """
    Generate the numerical labels from .geo text 
      ("ideal", "weight", "residual", etc)
      Assert consistent with n_atom expectation
    """
    ln, labels = self.lines[-2]
    numerical_labels = labels.split()
    if self.n_values:
      # test labels
      n_values =  len(numerical_labels)
      assert n_values == self.n_values or n_values == self.n_values-1, (numerical_labels,n_values,self.n_values)
    return numerical_labels

  @property
  def numerical_values(self):
    """
    Generate the numerical values from .geo text 
      Assert consistent with n_atom expectation
    """
    ln, values = self.lines[-1]
    values =  values.split()
    if self.n_values:
      n_values =  len(values)
      assert n_values == self.n_values or n_values == self.n_values-1, (values,n_values,self.n_values)
    return values
  
  @property
  def numerical_dict(self):
    """
    Provide a dictionary mapping numerical_labels: numerical_values
      So to get residual value, can do numerical_dict["residual"]
    """
    if not self._numerical_dict:
      self._numerical_dict = dict(zip(self.numerical_labels,self.numerical_values))
    return self._numerical_dict

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


  @staticmethod
  def _is_int(val):
    try:
      _ =int(val)
      return True
    except:
      return False

### Start of Entry subclasses


class NonBondedEntry(Entry):
  name = "nonbonded"
  entry_class_trigger = "Nonbonded interactions"
  entry_trigger = "nonbonded"
  n_atoms = 2
  n_values = 3
  internals = "bonds"
  origin_label_func = "get_bond_origin_id_labels" 

  
class AngleEntry(Entry):
  name = "angle"
  entry_class_trigger = "Bond angle restraints"
  entry_trigger = "angle"
  n_atoms = 3
  n_values = 6
  origin_label_func = "get_angle_origin_id_labels"
  internals = "angles"
  restraint_class = angle

    
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
  origin_label_func = "get_bond_origin_id_labels" 
  restraint_class = bond
  
  def to_proxy(self):
    d = self.numerical_dict
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
  internals = "angles"
  # internals = [
  #     'dihedrals',
  #      'torsion',
  # ]
  origin_label_func = "get_dihedral_origin_id_labels" 
  restraint_class = dihedral
    
  @property
  def is_harmonic(self):
    return "harmonic" in self.numerical_labels
  
  @property
  def is_sinusoidal(self):
    return "sinusoidal" in self.numerical_labels

  @property
  def periodicity(self):
    if self.is_harmonic:
      return int(self.numerical_dict["harmonic"])
    else:
      return int(self.numerical_dict["sinusoidal"])
    
  def to_proxy(self):
    proxy = geometry_restraints.dihedral_proxy(
      i_seqs=self.i_seqs,
      angle_ideal=self.ideal,
      weight=self.weight,
      periodicity=self.periodicity,
      alt_angle_ideals=None,
      origin_id=self.origin_id)
    return proxy


class CBetaEntry(DihedralEntry):
  name = "c-beta"
  entry_class_trigger = "C-Beta improper torsion"



class ChiralityEntry(Entry):
  name = "chirality"
  entry_class_trigger= "Chirality restraints"
  entry_trigger = "chirality"
  n_atoms = 4
  n_values = 7
  internals = "chirals"
  origin_label_func = "get_chiral_origin_id_labels"
  restraint_class = chirality
    
  @property
  def both_signs(self):
    return bool(self.numerical_dict["both_signs"])
    
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
  restraint_class = planarity
  internals = "planes"

  def __init__(self,*args,**kwargs):
    super().__init__(*args,**kwargs)
    self._compositional_labels = []
    self._compositional_values = []
    self._numerical_labels = []
    self._numerical_values = []
    self._prepared = False

  def _prepare(self):
    """
    Interpret lines from a Plane entry. Save results as hidden attributes
    Sets:
      self._compositional_labels
      self._compositional_values
      self._numerical_labels
      self._numerical_values
    """
    if not self._prepared:


      nums = [[None]*self.n_values for l in range(len(self.lines)-1)]
      for i,(ln,line) in enumerate(self.lines[1:]):
        line = line.replace(self.name,"")
        parts = shlex.split(line)
        comp_value = parts[0]
        stem = "atom"
        
        comp_label = f"{stem}_{i}"
        comp_value = comp_value.replace("pdb=","")
        comp_value = self._remove_outer_quotes(comp_value)
        self._compositional_labels.append(comp_label)
        self._compositional_values.append(comp_value)
        for j,p in enumerate(parts[1:]):
          nums[i][j] = p

      ln,line = self.lines[0]
      self._numerical_labels = line.strip().split()
      
      # fill empty values down columns
      for i,row in enumerate(nums):
        if i>0:
          nums[i][-2] = nums[0][-2]
          nums[i][-1] = nums[0][-1]
      nums_T = list(map(list, zip(*nums))) # transpose
      self._numerical_values = nums_T
      self._prepared = True


  @property
  def record(self):
    """
    A record is a dictionary representation with with simple 
      key:value attributes. 
    """
    records = []
    for i in range(self.n_atoms):
      d = {}
      if self._i_seqs:
        d.update(self._i_seqs)
      if self._id_strs:
        d.update(self._id_strs)
      d.update({
        self.compositional_labels[i]:self.compositional_values[i]
          })
      for key,value in self.numerical_dict.items():
        d[key] = value[i]
      d["proxy_seq"] = self.proxy_seq
      records.append(d)
    return records

  @property
  def n_atoms(self):
    # Override for planes
    return len(self.compositional_values)

  @property
  def compositional_labels(self):
    # Override for planes
    self._prepare()
    return self._compositional_labels

  @property
  def compositional_values(self):
    # Override for planes
    self._prepare()
    return self._compositional_values
      
  @property
  def numerical_labels(self):
    # Override for planes
    self._prepare()
    return self._numerical_labels

  @property
  def numerical_values(self):
    # Override for planes
    self._prepare()
    return self._numerical_values

  @property
  def weights(self):
    return [float(w) for w in self.numerical_dict["weight"]]
    
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
  origin_label_func = "get_parallelity_origin_id_labels"
  restraint_class = parallelity
  parsing_trigger = "plane 1"

  def __init__(self,*args,**kwargs):
    super().__init__(*args,**kwargs)
    self._compositional_labels = []
    self._compositional_values = []
    self._numerical_labels = []
    self._numerical_values = []
    self._prepared = False
    self._j_seqs = None


  def _prepare(self):
    """
    Interpret lines from a Plane entry. Save results as hidden attributes
    Sets:
      self._compositional_labels
      self._compositional_values
      self._numerical_labels
      self._numerical_values
    """
    if not self._prepared:
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

      self._compositional_labels = comp_labels
      self._compositional_values = comp_values
      self._numerical_labels = num_labels
      self._numerical_values = num_values
      self._prepared = True

  @property
  def i_seqs(self):
    """
    Provide the i_seqs as a list
    """
    if self.labels_are_i_seqs and not self._i_seqs:
      self._form_i_seqs(self.compositional_dict)
    if self._i_seqs:
      return [int(v) for k,v in self._i_seqs.items() if "_i_" in k]
  
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
    _ = self.i_seqs
    if self._i_seqs:
      return [int(v) for k,v in self._i_seqs.items() if "_j_" in k]

  def _form_id_strs(self,compositional_dict):
    self._id_strs =  {}
    for k,v in compositional_dict.items():
      if "_i_" in k:
        suffix = "_i_"+k.split("_i_")[-1]
      elif "_j_" in k:
        suffix = "_j_"+k.split("_j_")[-1]
      self._id_strs[f"id_str{suffix}"] = v
  @property
  def id_strs(self):
    if self._id_strs:
      return list(self._id_strs.values())

  @property
  def id_strs_i(self):
    if self._id_strs:
      return [v for k,v in self._id_strs.items() if "_i_" in k]

  @property
  def id_strs_j(self):
    if self._id_strs:
      return [v for k,v in self._id_strs.items() if "_j_" in k]

  @property
  def n_atoms(self):
    return len(self.compositional_values)

  @property
  def compositional_labels(self):
    self._prepare()
    return self._compositional_labels

  @property
  def compositional_values(self):
    self._prepare()
    return self._compositional_values
      
  @property
  def numerical_labels(self):
    self._prepare()
    return self._numerical_labels

  @property
  def numerical_values(self):
    self._prepare()
    return self._numerical_values

  @property
  def record(self):
    """
    A dictionary representation of an entry
      Atom labels are split up into single atom key:value pairs
    """
    d = {}
    # for planarity, need to specify i or j suffix
    if self.i_seqs:
      d.update(self._i_seqs)
    if self.id_strs:
      d.update(self._id_strs)
    d.update(self.compositional_dict)
    d.update(self.numerical_dict)
    return d
    
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
    CBetaEntry,
  ]
  def __init__(self,restraint_entry_classes=None):
    """
     Initialize with a list of Entry subclasses
    """
    if not restraint_entry_classes:
      restraint_entry_classes = self.entry_classes_default
    self.current_entry = None
    self.current_entry_class = None
    self.entry_classes = restraint_entry_classes
    self.entry_names = [entry_class.name for entry_class in self.entry_classes]
    self.entry_class_trigger_dict = {
      entry_class.entry_class_trigger:entry_class for entry_class in self.entry_classes
    }
    self.entry_triggers = [entry_class.entry_trigger for entry_class in self.entry_classes]
    self.entries = defaultdict(list)
    self._proxies = None
    self._model = None
    self.lines = None
    self.line_labels = []
    self.current_origin_id = None

  
  @property
  def origin_header_to_id_all(self):
    """
    Coalesce origin_header_to_id dicts from all entry types into one
    """
    cls = self.__class__
    if not hasattr(cls,"_origin_header_to_id_all"):
      cls._origin_header_to_id_all = {}
      for entry_class in self.entry_classes:
        cls._origin_header_to_id_all.update(entry_class.origin_header_to_id)
      if None in cls._origin_header_to_id_all:
        del cls._origin_header_to_id_all[None]
    return cls._origin_header_to_id_all
  
  @property
  def origin_headers(self):
    """
    The keys for origin_header_to_id_all, but with
      capitalized entries first for searching. A small optimization
    """
    if not hasattr(self,"_origin_headers"):
      self._origin_headers = []
      for key in self.origin_header_to_id_all.keys():
        if key[0].isupper():
          self._origin_headers.insert(0,key)
        else:
          self._origin_headers.append(key)
    return self._origin_headers

  def fill_labels_from_model(self,model):
    """
    Add i_seq or id_str attributes for each atom on each entry
    """
    # make i_seq:id_str mapping
    map_iseq_to_idstr = {
      atom.i_seq:Entry._remove_outer_quotes(atom.id_str().replace("pdb=","")) 
      for atom in model.get_atoms()}
    # Reverse
    map_idstr_to_iseq = {id_str:i_seq for i_seq,id_str in map_iseq_to_idstr.items()}

    for entry in self.entries_list:

      # Case where entry was loaded with i_seqs, can add id_strs
      if entry.labels_are_i_seqs:
        entry._form_id_strs({k:map_iseq_to_idstr[int(iseq)] for k,iseq in entry._i_seqs.items()})
        continue

      # Case where entry was loaded with id_strs, can load i_seqs
      entry.labels_are_id_strs = True
      for val in entry.compositional_values:
        if val not in map_idstr_to_iseq:
          entry.labels_are_id_strs = False
        
      if entry.labels_are_id_strs:
        entry._form_i_seqs({k:map_idstr_to_iseq[idstr] for k,idstr in entry._id_strs.items()})

  @property
  def model(self):
    """
    The mmtbx model manager that is associated with this geometry.
      Used primarily to determine i_seqs
    """
    return self._model

  @model.setter
  def model(self,value):
    """
    When setting a model, will attempt to fill out i_seqs/id_strs 
      automatically.
    """
    value.add_crystal_symmetry_if_necessary()
    self.fill_labels_from_model(value)
    self._model = value

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

  def parse_str(self,geo_str):
    """
    Parse .geo file lines
    """
    self.parse_lines(geo_str.split("\n"))

  def parse_file(self,filename):
    """
    Parse .geo file
    """
    with open(filename,"r") as fh:
      lines = fh.readlines()
      self.parse_lines(lines)
      
  def parse_lines(self,lines):
    """
    Parse .geo file lines
    """
    self.lines = lines
    self._parse()
  

  def _parse(self):
    started_first_entry = False # toggled after each new entry type
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
      if self.current_entry and not entry_trigger:
        self.current_entry.lines.append((i,line))
        last_line_label = "data"
        continue
      

      if entry_trigger:
        # Start new entry
        self._end_entry()
        last_line_label = "entry_trigger"
        self.current_entry = self.current_entry_class(origin_id=self.current_origin_id)
        if isinstance(self.current_entry,PlaneEntry):
          j = i-1 # The prior line
          self.current_entry.lines.append((j,self.lines[j]))
        
        self.current_entry.lines.append((i,line))
        started_first_entry = True
        continue

      # Check for start of new entry section
      entry_class_trigger = self._startswith_plural(line, self.entry_class_trigger_dict.keys())
      if entry_class_trigger:
        # Start a new entry class
        last_line_label = "entry_class_trigger"
        self._end_entry()
        entry_class = self.entry_class_trigger_dict[entry_class_trigger]
        self.current_entry_class = entry_class
        started_first_entry = False
        continue

      # Check origin id trigger
      origin_trigger = self._startswith_plural(line, self.origin_headers)
      if origin_trigger:
        last_line_label = "origin_trigger"
        origin_id = self.origin_header_to_id_all[origin_trigger]
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
    if self.current_entry:
      if isinstance(self.current_entry,PlaneEntry):
        self.current_entry.lines = self.current_entry.lines[:-1]
      self.current_entry.proxy_seq = len(self.entries[self.current_entry_class.name])
      self.entries[self.current_entry_class.name].append(self.current_entry)
      self.current_entry = None

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

  def debug_print(self,start=0,stop=-1):
    if stop == -1:
      stop = len(self.lines)
    for i in range(start,stop):
      line = self.lines[i]
      label = self.line_labels[i]
      print(i,label,line)

  @staticmethod
  def _replace_idstr_with_int(text,max_int=100):
    """
    Replace id_strs in a geo_file str with integers
      For debugging.
    """
    index = 0
    def replacement(match):
      nonlocal index # allow accessing index
      original_length = len(match.group(0))
      replacement_text = f'{index}'.ljust(original_length)
      index += 1
      if index>max_int:
        index = 0
      return replacement_text

    new_text = re.sub(r'pdb=".*?"', replacement, text)
    return new_text

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

  @property
  def proxies(self):
    """
    Convert the entry objects to cctbx proxy objects. 
      Collect into a dict of lists of proxies.
    """
    if not self.model and not self.labels_are_i_seqs: 
      raise Sorry("Cannot build proxies without associating a model. Use container.model = model")
    if not self._proxies:
      self._proxies = defaultdict(list)
      for entry_name,entries in self.entries.items():
        if len(entries)>0:
          entry_class = entries[0].__class__
        if hasattr(entry_class,"to_proxy") and not hasattr(entry_class,"ignore"):
          self._proxies[entry_class.name] = [entry.proxy for entry in entries]
    return self._proxies
