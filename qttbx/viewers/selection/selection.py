"""
The selection data structure to represent selections in the GUI in a multi-syntax way
"""
import json
import re
from dataclasses import dataclass, asdict
from typing import List, Optional

from qttbx.viewers.selection.parser import (
  Parser,
  attrs_map_to_phenix,
  core_map_to_mmcif,
  is_compatible_string,
  attrs_map_to_mmcif
)

from cctbx.array_family import flex
from iotbx.pdb.atom_selection import selection_string_from_selection



class Selection:
  """
  This class:
    1. Provides a container to store both an mmtbx.model.manager along with selections.
      The selections can take multiple forms:
        a. Phenix string selections
        b. boolean atom selections (convertable with integer selection)
        c. Atom record selections (key:value combinations that match a single atom)
        d. Molstar selection syntax
        e. Chimera selection syntax (incomplete)

    2. Stores selections as a sequential history. Meaning this code snippet: 

          selection = Selection(model) # implicit selection of 'all'
          selection = selection.select_from_string("chain A")
          selection = selection.select_from_string("resname GLY")

        Results in a selection object equivalent to:

          model.selection("chain A and resname GLY")

        But with the important difference that selection.bool or selection.int remain applicable 
          to the original model, rather than only the prior selected sub-model

  """
  model = None
  blanks = set([""," ",".","?"])


  def __init__(self,model=None,data=None,selection_string=None,selection_bool=None):
    """
    Initialize a selection from a string or boolean selection, and optionally a molecular model

    Args:
      model (mmtbx.model.manager): The model the selection references
      data (dict): Selection history. Only populated when instantiating a selection from a previous selection.
      selection_string (str): The Phenix string selection
      selection_bool (flex.bool): A boolean atom selection array

    Returns:
      Selection: instance of this class
    """

    assert model, "Must instantiate with mmtbx.model.manager"
    self.model = model
    self.is_validated = False
    if not selection_bool and not selection_string:
      selection_string = "all"
    elif not selection_bool:
      selection_bool = self.model.selection(selection_string)
    
    assert selection_bool.size() in [0,self._n_atoms_initial],(
      f"Cannot provide selection of different size ({selection_bool.size()}) than number of atoms ({self._n_atoms_initial})"
    )

    
    if not data:
      data = {"selections":[]}
    self.data = data
    self.data["selections"].append(
        {
            "phenix_string":selection_string,
            "selection_bool":selection_bool,
        }
    )
    compatible, fail_reason, sel_str = self.validate()
    self.data["selections"][-1]["phenix_string"] = sel_str
  
  def __eq__(self,other):
    if not isinstance(other,self.__class__):
      return False
    a = list(self.bool)
    b = list(other.bool)
    if a==b:
      print("Selections are equal")
      return True
    return False
  def __len__(self):
    return self.bool.count(True)
  @property
  def _h(self):
    # Shortcut to the iotbx.pdb.hierarchy that fails silently if model not set
    return self.model.get_hierarchy()

  @property
  def _atoms(self):
    # Shortcut to the iotbx shared atoms that fails silently if model not set
    return self._h.atoms()
  @property
  def _n_atoms_initial(self):
    # Shortcut to the total number of atoms if model is set
    return self.model.get_number_of_atoms()

  @classmethod
  def from_selection_string(cls,selection_string=None,model=None):
    """
    New instance of Selection from a selection string and optionally a model.
    """
    return cls(model=model,
              selection_string=selection_string)

  def select_from_string(self,selection_string):
    """
    Select a new selection from an existing selection using a selection string
    """
    return self.__class__(self.model,
                          data=self.data,
                          selection_string=selection_string)

  @classmethod
  def from_selection_bool(cls,selection_bool=None,model=None):
    """
    New instance of Selection from a boolean atom selection and optionally a model.
    """
    return cls(model=model,
              selection_bool=selection_bool)
  
  def select_from_bool(self,selection_bool):
    """
    Select a new selection from an existing selection using a boolean atom selection
    """
    if selection_bool.size() != self._n_atoms_initial:
      assert selection_bool.size() == self.int.size(), (
        "If selection bool does not equal original number of atoms, must equal number of atoms selected")
      sel_bool = flex.bool(len(self.int),False)
      sel_bool.set_selected(selection_bool,True)
      sel_int = self.int.select(sel_bool)
      return self.select_from_int(sel_int)
    return self.__class__(self.model,
                          data=self.data,
                          selection_bool=selection_bool)
  
  @classmethod
  def from_selection_int(cls,selection_int=None,model=None):
    """
    New instance of Selection from an integer atom selection and optionally a model.
    """
    sel_bool = flex.bool(self._n_atoms_initial,False)
    sel_bool.set_selected(selection_int,True)
    return cls(model=model,
              selection_bool=selection_bool)

  def select_from_int(self,selection_int):
    """
    Select a new selection from an existing selection using an integer atom selection
    """
    sel_bool = flex.bool(self._n_atoms_initial,False)
    sel_bool.set_selected(selection_int,True)
    return self.__class__(self.model,
                          data=self.data,
                          selection_bool=sel_bool)    
  @classmethod
  def from_empty(cls,model=None):
    return cls.from_selection_bool(selection_bool=model.selection(""),model=model)


  @staticmethod
  def atom_records_to_phenix_string(atom_records):
    # Create a set of preferred keys to avoid duplicates
    preferred_keys = {
        'auth_asym_id',    # for chain
        'label_comp_id',   # for resname
        'auth_seq_id',     # for resseq
        'label_atom_id'    # for name
    }
    
    # Filter by keys - only keep preferred keys that are in core_map_to_mmcif
    atom_records_reduced = [{k:v for k,v in atom_dict.items() 
                        if k in attrs_map_to_mmcif.values()} 
                       for atom_dict in atom_records]
    
    # Filter by values
    atom_records_reduced = [{k:v for k,v in atom_dict.items() 
                           if v not in Selection.blanks} 
                          for atom_dict in atom_records_reduced]
    # Remove id
    atom_records_reduced = [{k:v for k,v in atom_dict.items() if  k != "id"} for atom_dict in atom_records_reduced]

    
    output = ""
    for r_i, record in enumerate(atom_records_reduced):
        output += "("
        for i, (k,v) in enumerate(record.items()):
            output += str(attrs_map_to_phenix[k])
            output += " "
            output += str(v)
            if i < len(record)-1:
                output += " and "
        output+=")"
        if r_i < len(atom_records_reduced)-1:
            output += " or "
    return output
    
  @classmethod
  def from_atom_records(cls,atom_records,model):
    """
    New instance of Selection from a list of atom records. An atom record is a group of 
      key:value associations using mmcif keys that is unique to a single atom. 

    Example atom record:
    record = {
      "asym_id":"C",
      "comp_id":"GLY",
      "seq_id":4,
      "atom_id":"CA",
      "alt_id": "A",
    }

    TODO: Insertion codes and segid are not currently supported.
    """
    if len(atom_records)==0:
      return cls.from_empty(model=model)

    phenix_string = cls.atom_records_to_phenix_string(atom_records)
    return cls.from_selection_string(phenix_string,model=model)


  @property
  def _selections_bool(self):
    # Return all the stored boolean selections
    bools = [sel["selection_bool"] for sel in self.data["selections"]]
    if bools.count(None)>0:
      return None
    return bools
      
  @property
  def bool(self):
    # Join boolean selection history with 'and' 
    if self.model:
      sels = self._selections_bool
      if sels:
        combined_sel = sels[0]
        for sel in sels[1:]:
          combined_sel &= sel
        return combined_sel
      else:
        s = self._debug_string
        return self.model.selection(s)


  @property
  def int(self):
    # convert boolean selection to integer
    return self.bool.iselection()
      

  @property
  def string(self):
    # If single string, return that. Else convert boolean selection to string
    if len(self.data["selections"])==1 and self.data["selections"][0]["phenix_string"]:
      return self.data["selections"][0]["phenix_string"]
    else:
      if self.bool.count(True)>0:
        return selection_string_from_selection(self._h,self.bool)
      else:
        return "none"
        
  @property
  def _debug_string(self):
    # Join string selections together 
    sel_strings = []
    for sel in self.data["selections"]:
        if not sel["phenix_string"]:
          b = sel["selection_bool"]
          if b.count(True)>0:
            sel["phenix_string"] = selection_string_from_selection(self._h,b)
          else:
            sel["phenix_string"] = "none"
        sel_strings.append(sel["phenix_string"])
    output = " and ".join([f"({sel_string})" for sel_string in sel_strings])
    return output
      
  def validate_debug(self):
    # Validate for consistency during tests

    # Make sure self.bool and self.string are consistent
    s_bool = self.bool
    s_string = self.model.selection(self.string)
    s_debug_string = self.model.selection(self._debug_string)
    assert (s_bool==s_string).count(True) == self._n_atoms_initial
    assert (s_bool==s_debug_string).count(True) == self._n_atoms_initial

    return True
    
  def validate(self):
    # Validate that a string is consistent with what PhenixParser can process. 
    #  Requires that self.model be set
    compatible, fail_reason, sel_str = is_compatible_string(self.string,self.model)
    self.is_validated = compatible
    self.fail_reason = fail_reason
    return compatible, fail_reason, sel_str

  @property
  def is_empty(self):
    return self.bool.count(True)==0

  # Molstar
  @property
  def molstar_syntax(self):
    # Convert and return this selection as molstar syntax
    key = "molstar_syntax"
    if key not in self.data:
      if not self.is_validated:
        compatible, fail_reason, sel_str = self.validate()
        if not compatible:
            return None
      parser = Parser.from_string(self.string)
      ast = parser.parse()
      self.data[key] = ast.molstar_syntax()
    
    return self.data[key]
  
  @property
  def phenix_syntax(self):
    # The default string selection is Phenix string selection
    return self.string
    
  @property
  def phenix_string(self):
    # The default string selection is Phenix string selection
    return self.string

  @property
  def pandas_syntax(self):
    # Convert and return this selection in a form that can select
    #   atoms from a Pandas dataframe with mmcif column labels
    key = "pandas_syntax"
    if key not in self.data:
      if not self.is_validated:
        compatible, fail_reason, sel_str = self.validate()
        if not compatible:
            return None
      parser = PhenixParser(self.string)
      self.data[key] = parser.pandas_string
    return self.data[key]

  @property
  def chimerax_syntax(self):
    raise NotImplementedError
