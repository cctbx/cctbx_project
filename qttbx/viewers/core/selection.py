

import json
import re
from dataclasses import dataclass, asdict
from typing import List, Optional

from .parameters import attrs_map_to_phenix, blanks, core_map_to_mmcif, logic_map_to_phenix, protein_comp_ids
from .parameters import params
from .parser import PhenixParser


@dataclass(frozen=True)
class Selection:
  """
  The unifying selection object
  """
  # Molstar code in string form, 
  #   to specify a StructureSelectionQuery object
  molstar_syntax: str

  # A selection string that runs with model.selection()
  phenix_string: str

  # A selection string that runs with df.query()
  pandas_string: str

  # Debug
  # Dictionary with enough key:value pairs to
  #   uniquely identify 1 atom
  atom_records: Optional[List[dict]] = None
  atom_records_reduced: Optional[List[dict]] = None
  parser: Optional[PhenixParser] = None # Only use for debug

      

  @classmethod
  def from_phenix_string(cls,phenix_string,debug=True):
    parser = PhenixParser(phenix_string)
    pandas_string = parser.pandas_string
    molstar_syntax = parser.molstar_syntax
    if not debug:
        parser = None
    return cls(phenix_string=phenix_string,
              pandas_string=pandas_string,
              molstar_syntax=molstar_syntax,
              parser=parser,
    )
  @classmethod
  def from_pandas_string(cls,pandas_string,debug=True):
    parser = PhenixParser(pandas_string)
    pandas_string = parser.pandas_string
    molstar_syntax = parser.molstar_syntax
    phenix_string = None # TODO: Make this in parser
    if not debug:
        parser = None
    return cls(phenix_string=phenix_string,
              pandas_string=pandas_string,
              molstar_syntax=molstar_syntax,
              parser=parser,
    )
  # @staticmethod
  # def convert_atom_records_to_pandas_string(atom_records):
  #   # filter by keys
  #   atom_list = [{k:v for k,v in atom_dict.items() if k in core_map_to_mmcif.values()} for atom_dict in atom_records]
  #   # filter by values
  #   atom_list  = [{k:v for k,v in atom_dict.items() if  v not in blanks} for atom_dict in atom_list]
  #   # Convert to pandas string
  #   sel_str = ""
  #   for atom_dict in atom_list:
  #     atom_prefix = ""
  #     if len(sel_str)>1:
  #       atom_prefix = "or"
  #     atom_str = "("
  #     for key,value in atom_dict.items():
  #       if value != "":
        
  #         if len(atom_str)>1:
  #           prefix = "and"
  #         else:
  #           prefix = ""
  #         if isinstance(value,str):
  #           atom_str+=f"{prefix} {key} == '{value}' " # quote
  #         else:
  #           atom_str+=f"{prefix} {key} == {value} "
  #     atom_str+=") "
  #     sel_str+=f"{atom_prefix} {atom_str}"
  #   return atom_list, sel_str

  @classmethod
  def from_atom_records(cls,atom_records,debug=True):
    # a list of dictionaries with key:value matches
    # Here build tree directly

    # filter by keys
    atom_records_reduced = [{k:v for k,v in atom_dict.items() if k in core_map_to_mmcif.values()} for atom_dict in atom_records]
    # filter by values
    atom_records_reduced = [{k:v for k,v in atom_dict.items() if  v not in blanks} for atom_dict in atom_records_reduced]

    # remove id, no way to access with phenix
    atom_records_reduced = [{k:v for k,v in atom_dict.items() if  k != "id"} for atom_dict in atom_records_reduced]


    root = Or()
    for d in atom_records_reduced:
      and_node = And()
      for k,v in d.items():         
        comparison = Comparison(k,"==",str(v))
        and_node.children.append(comparison)
      root.children.append(and_node)
    tree = SelectionTree(root=root)

    pandas_string = tree.pandas_string
    molstar_syntax = tree.molstar_syntax
    phenix_string = tree.phenix_string
    return cls(phenix_string=phenix_string,
              pandas_string=pandas_string,
              molstar_syntax=molstar_syntax,
              atom_records = atom_records, # store the raw input
              atom_records_reduced = atom_records_reduced,
              parser=None,
    )    

  @classmethod 
  def from_i_seqs(cls,sites,i_seqs):
    sel_sites = sites.select_from_i_seqs(i_seqs)
    atom_records = sel_sites.to_records_compositional()
    return cls.from_atom_records(atom_records)

  def to_dict(self):
    d = asdict(self)
    # DEBUG: can be None?
    d["atom_list"] = [""]
    d["parser"] = None
    return d

  def to_json(self,indent=2):
    return json.dumps(self.to_dict(),indent=indent)

