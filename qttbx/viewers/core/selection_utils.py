

import json
import re
from dataclasses import dataclass, asdict
from typing import List, Optional

import networkx as nx
import numpy as np



"""
Selection API

Rationale:
The basic idea is to define a common way to specify selections among all the relevant other programs.
Each program will have corner cases that can't be translated. This api serves as a filter, if a selection can be
formatted using this api, it should be translatable to any other relevant program.

Basics:
1. The universal language is mmcif keys and logic understood by df.query()
    In cases where it is ambiguous and unspecified, (auth_ vs label_) there is a default choice.
"""

class SelectionQuery:
  pass # TODO: remove

@dataclass(frozen=True)
class Selection:
  # Molstar code in string form, 
  #   to specify a StructureSelectionQuery object
  molstar_syntax: Optional[str] = None

  # A selection string that runs with model.selection()
  phenix_string: Optional[str] = None

  # A selection string that runs with df.query()
  pandas_string: Optional[str] = None

  # Dictionary with enough key:value pairs to
  #   uniquely identify 1 atom
  atom_list: Optional[List[dict]] = None 

  


  @classmethod
  def from_phenix_string(cls,phenix_string):
    parser = PhenixParser(phenix_string)
    pandas_string = parser.pandas_string
    molstar_syntax = parser.molstar_syntax
    return cls(phenix_string=phenix_string,
              pandas_string=pandas_string,
              molstar_syntax=molstar_syntax,
    )
  def to_dict(self):
    return asdict(self)
  def to_json(self,indent=2):
    return json.dumps(self.to_dict(),indent=indent)

class SelConverterPhenix:
  """
  A basic conversion between selection syntax that doesn't require
  an abstract syntax tree. (Just simple word replacement)
  """
  def __init__(self):
      # Bidirectional phenix and pandas
      self.bidirectional_map = {
          "element":"type_symbol",
          'chain': 'asym_id',
          'resseq': 'seq_id',
          'name': 'atom_id',
          'bfactor': 'B_iso_or_equiv',
          'bfac': 'B_iso_or_equiv',
          'occupancy': 'occupancy',
          'resname': 'comp_id',
          'altloc': 'alt_id',
          'or': '|',
          'and': '&',
          'not': '~',
      }
      # Unidirectional mapping from Pandas to Phenix
      self.unidirectional_pandas_to_phenix = {
          '==': '',  # Remove '==' when converting from Pandas to Phenix
      }
      # Combining maps for specific direction
      self.pandas_to_phenix_map = {**{v: k for k, v in self.bidirectional_map.items()}, **self.unidirectional_pandas_to_phenix}
      auth_label_supplement = {}
      for key,value in self.pandas_to_phenix_map.items():
        auth_label_supplement["auth_"+value] = key
        auth_label_supplement["label_"+value] = key
      self.pandas_to_phenix_map.update(auth_label_supplement)

      # Phenix to pandas
      #self.phenix_to_pandas_map = self.bidirectional_map

      # Common
      #self.phenix_keyword_map_to_common = {k:v for k,v in copy.deepcopy(self.bidirectional_map).items() if k not in ['and','or','not']


  def replace_keys_in_str(self, input_str, replacements):
    # Separate dictionary entries into word and non-word items
    word_replacements = {}
    non_word_replacements = {}
    for key, value in replacements.items():
        # Check if key consists entirely of word characters
        if re.match(r'^\w+$', key):
            word_replacements[key] = value
        else:
            non_word_replacements[re.escape(key)] = value

    # Create patterns for word and non-word replacements
    # Word replacements use word boundaries
    if word_replacements:
        word_pattern = r'\b(' + '|'.join(re.escape(k) for k in word_replacements) + r')\b'
        input_str = re.sub(word_pattern, lambda match: word_replacements[match.group(0)], input_str)

    # Non-word replacements do not use word boundaries
    if non_word_replacements:
        non_word_pattern = '|'.join(non_word_replacements.keys())
        input_str = re.sub(non_word_pattern, lambda match: non_word_replacements[re.escape(match.group(0))], input_str)

    return input_str

  def phenix_postprocess(self,s):
    s = self.replace_dots(s)
    s = self.unquote_string(s)
    return s

  @staticmethod
  def unquote_string(s):
    return s.replace("'","").replace('"','')

  @staticmethod
  def replace_dots(input_string):
    # Replace single-quoted and double-quoted periods
    modified_string = input_string.replace("'.'", "'*'")
    modified_string = modified_string.replace('"."', '"*"')
    return modified_string

  def convert_phenix_to_pandas(self, sel_str):
      return self.replace_keys_in_str(sel_str, self.phenix_to_pandas_map)

  def convert_pandas_to_phenix(self, sel_str):
      return self.phenix_postprocess(self.replace_keys_in_str(sel_str, self.pandas_to_phenix_map))

  def convert_phenix_to_common(self,sel_str):
    core_names_str =  self.replace_keys_in_str(sel_str,self.phenix_keyword_map_to_common)
    return self.replace_keys_in_str(core_names_str,core_keys_to_mmcif_keys_default) # TODO: don't hardcode default

  def convert_common_to_phenix(self,sel_str):
    return self.phenix_postprocess(self.replace_keys_in_str(sel_str,self.reversed_kw_map(self.phenix_keyword_map_to_common)))


def df_nodes_group_to_query(df,composition_cols):
  # A new query building function, here we ignore seq ranges for simplicity

  def to_python_type(obj):
    if isinstance(obj, np.generic):
      return obj.item()
    return obj


  selection_dicts = []
  is_single = False
  if len(composition_cols)==1:
    composition_cols = composition_cols[0]
    is_single = True
  for name,group in df.groupby(composition_cols):
    if is_single:
      name = [name]
      labels = [composition_cols]
    else:
      labels = composition_cols
    d = dict(zip(labels,name)) # simple key value dict
    d = {k: {'ops': [{'op': '==', 'value': to_python_type(value)}]} for k,value in d.items() if value != "*"}
    selection_dicts.append(d)

  selections = []
  for data in selection_dicts:
    selections.append(Selection(data))

  query = SelectionQuery(selections=selections)
  return query


def form_simple_str_common(df_sel_simple):

  final_component_strs = []
  for row in df_sel_simple.itertuples(index=False):
    s = f'(asym_id {row.asym_id} and seq_id {row.seq_id_start}:{row.seq_id_stop} and atom_id {row.atom_id})'
    final_component_strs.append(s)

  sel_str_simple_common = ' or '.join(final_component_strs)
  return sel_str_simple_common
