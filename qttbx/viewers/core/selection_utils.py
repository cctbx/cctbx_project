

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

  # def __post_init__(self):
  #   # The parser does not read molstar code, need to 
  #   #   specify at least one other form of input
  #   assert [self.phenix_string,self.pandas_string,self.atom_list].count(None) >0, (
  #     "Specify at least one form of input"
  #   )
  #   if self.phenix_string is not None:
      
  #     self.pandas_string = parser.pandas_string
  #     self.molstar_syntax = parser.molstar_syntax
  #   elif self.pandas_string is not None:
  #     raise NotImplementedError
  #   elif self.atom_list is not None:
  #     raise NotImplementedError
    
    



# class SelectionQuery:
#   """
#   A class to specify a selection. 
#   """
#   default_params = DotDict({"keymap":core_map_to_mmcif,'refId':'','model_id':''})

#   @staticmethod
#   def get_default_params():
#     return copy.deepcopy(SelectionQuery.default_params)

#   def __init__(self,selections,params=None):

#     # set up params
#     if params is None:
#       params = self.get_default_params()
#     for key,value in params.items():
#       setattr(self,key,value)
#     self.params = DotDict(params)


#     self.selections = selections

#     # DEBUG
#     for sel in self.selections:
#       for key,value in list(sel.data.items()):
#         if key not in mmcif_keys_to_core_keys_default:
#           del sel.data[key]



#   def __eq__(self,other):
#     assert isinstance(other,self.__class__), 'Compare with same class'
#     return self.to_json()==other.to_json()

#   def print(self):
#     print(self.to_json(indent=2))

#   def to_dict(self):
#     d = {'selections':[s.resolve_keys(self).to_dict() for s in self.selections],"params":self.params}
#     return d

#   def to_json(self,indent=None):
#     d = self.to_dict()
#     return json.dumps(d,indent=indent)

#   @classmethod
#   def from_model_ref(cls,ref):
#     params = cls.get_default_params()
#     for key,value in ref.external_ids.items():
#       params[key] = value
#     params.refId = ref.id
#     selections = [Selection.from_default()]

#     return cls(selections=selections,params=params)


#   @classmethod
#   def from_i_seqs(cls,atom_sites,i_seqs):
#     try:
#       i_seqs = [int(i) for i in i_seqs]
#       sel_sites = atom_sites.__class__(atom_sites.iloc[i_seqs],attrs_map=atom_sites.attrs_map)
#       # Use machinery in atom_sites class. TODO: Move that here.
#       query = atom_sites._convert_sites_to_query(sel_sites)
#     except:
#       print("Failed to load from i_seqs:",i_seqs)
#     return query

#   @classmethod
#   def from_json(cls,json_str,ref_id=None):
#     d = json.loads(json_str)
#     assert 'selections' in d, 'Must have "selections" entry in json input'
#     if 'params' in d:
#       params = d['params']
#     else:
#       params = None
#     if ref_id is not None:
#       params["refId"] = ref_id
#     selections = [Selection(s) for s in d['selections']]
#     return cls(selections,params=params)

#   @classmethod
#   def from_all(cls,ref_id=None):
#     selection = Selection.from_default()
#     selections = [selection]
#     params = cls.get_default_params()
#     if ref_id is not None:
#       params["refId"] = ref_id
#     return cls(selections,params=params)

#   @property
#   def pandas_query(self):
#     full_query_list = []
#     for selection in self.selections:
#       query_list = []
#       for keyword, conditions in selection.data.items():
#         keyword_conditions = []
#         for op in conditions['ops']:
#           operator = op['op']
#           value = op['value']
#           if value == '*':
#             continue  # skip any
#           if isinstance(value, str):
#             value = f"'{value}'"  # Quote strings
#           keyword_conditions.append(f'{keyword} {operator} {value}')
#         if keyword_conditions:  # Only append if keyword_conditions is not empty
#           keyword_query = ' & '.join(keyword_conditions)
#           query_list.append(f'({keyword_query})')
#       if query_list:  # Only append if query_list is not empty
#         query = ' & '.join(query_list)
#         full_query_list.append(f'({query})')

#     if len(full_query_list)==0:
#       output =  "index >= 0" # all
#     else:
#       output =  ' | '.join(full_query_list)
#     return output


#   @property
#   def phenix_string(self):
#     # Use similar logic as for pandas query
#     converter = SelConverterPhenix()
#     return converter.convert_pandas_to_phenix(self.pandas_query)

#   @property
#   def common_string(self):
#     converter = SelConverterPhenix()
#     phenix_string =  converter.convert_pandas_to_phenix(self.pandas_query)
#     common_string = converter.convert_phenix_to_common(phenix_string)
#     return common_string
#   @property
#   def pandas_string(self):
#     return self.pandas_query # alias

# class Selection:
#   default_select_all = {
#     # keyword, op, value
#     #'entity_id':{'ops':[{'op': '==', 'value': '*'}]},
#     'asym_id':{'ops':[{'op': '==', 'value': '*'}]},
#     'comp_id':{'ops':[{'op': '==', 'value': '*'}]},
#     'seq_id':{'ops':[{'op': '==', 'value': '*'},{'op': '==', 'value': '*'}]},
#     'atom_id':{'ops':[{'op': '==', 'value': '*'}]},
#     # 'id':{'ops':[{'op': '==', 'value': '*'}]},
#     }

#   def __init__(self,data):
#     # data is a dict representing a selection, same structure as default_select_all
#     self.data = data
#     # check case for consistency
#     for keyword in ['comp_id','atom_id',"alt_id","asym_id","type_symbol"]:
#       if keyword in self.data:
#         for op in self.data[keyword]["ops"]:
#           op['value'] =  op['value'].upper()

    
#   # def _rename_keys_for_auth_label(self,d,auth_label_pref='auth'):
#   #   new_columns = {}
#   #   for col in d.keys():
#   #     if 'auth_' in col and auth_label_pref == 'label':
#   #       continue
#   #     elif 'label_' in col and auth_label_pref == 'auth':
#   #       continue
#   #     elif 'auth_' in col and auth_label_pref == 'auth':
#   #       newcol = col.replace('auth_','')
#   #     elif 'label_' in col and auth_label_pref == 'label':
#   #       newcol = col.replace('label_','')
#   #     else:
#   #       newcol = col
#   #     new_columns[col] = newcol

#   #   out = {}
#   #   for key,value in d.items():
#   #     if 'auth' in key and  auth_label_pref == 'label':
#   #       pass # skip
#   #     elif 'label' in key and auth_label_pref == 'auth':
#   #       pass
#   #     else:
#   #       out[new_columns[key]] = value
#   #   return out


#   @classmethod
#   def from_default(cls):
#     data = cls.default_select_all
#     return cls(data)

#   @classmethod
#   def from_df_record(cls,record,params=None):
#     d = {}
#     for key1,value1 in record.items():
#       if isinstance(value1,dict):

#         for key2,value2 in value1.items():
#           if key2 == 'ops':
#             d[key1] = value1
#             continue
#       # add op entries if not present (probably not)
#       d[key1] = {'ops':[{'op':'==','value':value1}]}

#     return cls(d,params=params)


#   def to_dict(self):
#     return self.data
  

#   @classmethod
#   def from_json(cls,json_str,params=None):
#     d = {"class_name":cls.__class__.__name__}
#     d.update(json.loads(json_str))

#     return cls(d,params=params)

#   def to_json(self,indent=None):
#     return json.dumps(self.to_dict(),indent=indent)

#   def resolve_keys(self,query):
#     # If a query is present, where auth_ label_ ambiguity is resolved (if needed),
#     # settle explicitly what the keys are
#     for k,v in list(self.data.items()):
#       if k in query.params.keymap:
#         newkey = query.params.keymap[k]
#         self.data[newkey] = self.data.pop(k)
#     return self
# @dataclass(frozen=True)
# class Selection:
#   asym_id: Dict[str, str] = field(default_factory=lambda: {'op': '==', 'value': '*'})
#   comp_id: Dict[str, str] = field(default_factory=lambda: {'op': '==', 'value': '*'})
#   seq_id: Dict[str, str] = field(default_factory=lambda: {'op': '==', 'value': '*'})
#   atom_id: str = '*'
#   id: Union[str, int] = '*'

#   def post_init(self):
#     # if self.seq_id not in ['*', '']:
#     #   if self.seq_id_start in ['*', '']:
#     #     self.seq_id_start = self.seq_id
#     #   if self.seq_id_stop in ['*', '']:
#     #     self.seq_id_stop = self.seq_id

#     for field_name in ['seq_id_start', 'seq_id_stop', 'id']:
#       field_value = getattr(self, field_name)
#       if field_value not in ['*', '']:
#         if not isinstance(field_value, int):
#           try:
#             int_value = int(str(field_value))
#           except ValueError:
#             raise ValueError(f'The field '{field_name}' must be convertible to int, or '*' or ''.')

#           setattr(self, field_name, int_value)

#   @classmethod
#   def from_dict(cls, init_dict: Dict[str, Union[str, int]], auth_label_pref: str = 'auth', fallback: bool = True):
#     instance = cls()

#     for field_name in cls.__annotations__.keys():
#       auth_key = f'auth_{field_name}'
#       label_key = f'label_{field_name}'

#       # If the field exists without any prefix, set it directly
#       if field_name in init_dict:
#         setattr(instance, field_name, init_dict[field_name])

#       # Otherwise, follow the logic for auth_ or label_ prefix
#       elif auth_key in init_dict and auth_label_pref == 'auth':
#         setattr(instance, field_name, init_dict[auth_key])
#       elif label_key in init_dict and auth_label_pref == 'label':
#         setattr(instance, field_name, init_dict[label_key])
#       elif fallback:
#         if auth_key in init_dict:
#           setattr(instance, field_name, init_dict[auth_key])
#         elif label_key in init_dict:
#           setattr(instance, field_name, init_dict[label_key])

#     instance.post_init()
#     return instance


#   def to_dict(self):
#     return asdict(self)

#   def to_json(self):
#     return json.dumps(self.to_dict())

# @dataclass(frozen=True)
# class SelectionQuery:
#   site_params: List[Selection] = field(default_factory=list)
#   model_key: str = '' # The key for the applicable model in State
#   viewer_key: str = '' # The key for the applicable model inside the external viewer. Often set by them.
#   auth_label_pref: str = 'auth'

#   @classmethod
#   def from_json(cls, json_str: str, auth_label_pref: str = 'auth', fallback: bool = True):
#     if len(json_str.strip())==0:
#       return None
#     d = json.loads(json_str)
#     return cls.from_dict(d,auth_label_pref=auth_label_pref,fallback=fallback)

#   @classmethod
#   def from_dict(cls, init_dict: Dict, auth_label_pref: str = 'auth', fallback: bool = True):
#     instance = cls()
#     instance.site_params = [Selection.from_dict(d, auth_label_pref, fallback) for d in init_dict.get('site_params', [])]
#     instance.model_key = init_dict.get('model_key', '')
#     instance.viewer_key = init_dict.get('viewer_key', '')
#     instance.auth_label_pref = init_dict.get('auth_label_pref', 'auth')
#     return instance

#   def to_javascript(self):
#     # replace with json, naming should be consistent across datatypes
#     s = repr(self)
#     s = s.replace('SelectionQuery(','{')
#     s = s.replace('Selection(','{')
#     s = s.replace(')','}')
#     s = s.replace('=',': ')
#     s = s.replace('selections','site_params') # TODO: refactor to consistent name
#     #s ='{query: '+s+'}'
#     s = s.replace('*','')
#     return s

#   def to_dict(self):
#     return asdict(self)

#   def to_json(self):
#     return json.dumps(self.to_dict())

#   def __repr__(self):
#     parsed_json = json.loads(self.to_json())
#     pretty_json_str = json.dumps(parsed_json, indent=2)
#     return pretty_json_str
# @dataclass(frozen=True)
# class Selection:
#   # A simple selection, a single value for each field
#   #entity_id: str = '*'
#   asym_id: str = '*'
#   comp_id: str = '*'
#   seq_id_start: Union[str, int] = '*' # must be convertable to int
#   seq_id_stop: Union[str, int] = '*'  # must be convertable to int
#   atom_id: str = '*'
#   id: Union[str, int] = '*' # must be convertable to int

#   def __post_init__(self):
#     for field_name in ['seq_id_start', 'seq_id_stop','id']:
#       field_value = getattr(self, field_name)
#       if not field_value in ['*','']:
#         if not isinstance(field_value,int):
#           try:
#             int_value = int(str(field_value))
#           except ValueError:
#             raise ValueError(f'The field '{field_name}' must be convertible to int, or all,none ('*', '') ')
#           setattr(self, field_name, int_value)

# @dataclass(frozen=True)
# class SelectionQuery:
#   # The highest level object to pass around to contain selection data
#   selections: List[Selection] = field(default_factory=list)
#   model_num: int = 1
#   auth_label_pref: str = 'auth'

#   def to_javascript(self):
#     # TODO: replace with json, naming should be consistent across datatypes
#     s = repr(self)
#     s = s.replace('SelectionQuery(','{')
#     s = s.replace('Selection(','{')
#     s = s.replace(')','}')
#     s = s.replace('=',': ')
#     s = s.replace('selections','site_params')
#     s ='{query: '+s+'}'
#     s = s.replace('*','')
#     return s


#   def to_dict(self):
#     return asdict(self)


#   @classmethod
#   def from_json(cls, json_str):
#     data = json.loads(json_str)
#     return cls(**data)


#   def to_json(self):
#     return json.dumps(self.to_dict())

#   def __repr__(self):
#     parsed_json = json.loads(self.to_json())
#     pretty_json_str = json.dumps(parsed_json, indent=2)
#     return pretty_json_str



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


# class SelConverterPhenix:
#   phenix_keyword_map_default = {
#         'chain':'asym_id',
#         'resseq':'seq_id',
#         'name':'atom_id',
#         'bfactor':'B_iso_or_equiv',
#         'bfac':'B_iso_or_equiv',
#         'occupancy':'occupancy',
#         'resname':'comp_id',
#         'altloc':'alt_id'
#       }

#   @staticmethod
#   def replace_keys_in_str(input_str, replacements):
#     # Sort keys by length in descending order to handle substrings correctly
#     sorted_keys = sorted(replacements.keys(), key=len, reverse=True)
#     pattern = r'\b(' + '|'.join(re.escape(key) for key in sorted_keys) + r')\b'
#     return re.sub(pattern, lambda match: replacements[match.group()], input_str)


#   def __init__(self,phenix_keyword_map=None):
#     if phenix_keyword_map is None:
#       phenix_keyword_map = self.phenix_keyword_map_default

#     # store forward and reverse mapping
#     self.phenix_keyword_map = phenix_keyword_map
#     self.phenix_keyword_map_rev = {v:k for k,v in phenix_keyword_map.items()}

#   def convert_phenix_to_common(self,sel_str):
#     return self.replace_keys_in_str(sel_str,self.phenix_keyword_map)

#   def convert_common_to_phenix(self,sel_str):
#     return self.replace_keys_in_str(sel_str,self.phenix_keyword_map_rev)



# def ast_to_pandas_query(ast):
#   # Handle the BinaryLogical nodes
#   if ast['type'] == 'BinaryLogical':
#     left = ast_to_pandas_query(ast['left'])
#     right = ast_to_pandas_query(ast['right'])
#     op = ast['op']['value']

#     if op == 'and':
#       op = '&'
#     elif op == 'or':
#       op = '|'

#     return f'({left} {op} {right})'

#   # Handle the Keyword nodes
#   elif ast['type'] == 'Keyword':
#     keyword = ast['keyword']['value']
#     value = ast['value']

#     if value['type'] == 'ID':
#       return f'{keyword} == '{value['value']}''
#     elif value['type'] == 'RANGE':
#       start, end = value['value'].split(':')
#       return f'{start}<={keyword}<={end}'
#     elif value['type'] == 'INT':
#       return f'{keyword} == {value['value']}'
#     elif value['type'] == 'WILDCARD':
#       return f'{keyword} == {keyword}'  # Match anything that's not NaN
#   # If none of the above, return an empty string
#   else:
#     return ''

def find_simplest_selected_nodes(graph, selected_atoms):
  '''
  graph: atom_sites.G, the nx hierarchy graph
  selected_atoms: atom 'id' field. Must be unique
  '''
  simplest_nodes = []

  # Convert the pandas Series to a set for faster lookup
  selected_atoms_set = set(selected_atoms)
  #print(f'Selected atom set: {sorted(list(selected_atoms_set))}')
  # Maintain a set of atoms that have already been covered by an included node
  covered_atoms = set()

  for node in nx.dfs_preorder_nodes(graph, source='root'):
    #print(f'Checking node: {node}')  # Debugging print

    downstream_atoms = set(graph.nodes[node].get('ids', []))
    #print(f'Downstream atoms: {list(downstream_atoms)[:6]} ... len: {len(list(downstream_atoms))}')  # Debugging print
    #print(f'Covered atoms: {list(covered_atoms)[:6]}... len: {len(list(covered_atoms))}')  # Debugging print

    # If all downstream atoms are in selected_atoms and none are covered yet
    if not downstream_atoms.issubset(selected_atoms_set):
      #print('Atoms in downstream_atoms not in selected_atoms_set:', list(downstream_atoms - selected_atoms_set)[:10])
      pass
    if downstream_atoms.issubset(selected_atoms_set) and not downstream_atoms & covered_atoms:
      #print(f'Adding node: {node}')  # Debugging print
      simplest_nodes.append(node)
      covered_atoms.update(downstream_atoms)

  return simplest_nodes

# def df_nodes_group_to_query(df):
#   '''
#   Nodes is a dataframe obtained from
#   simplifying a selection with a graph. Each
#   row is a selection, which is possible grouped by seq_id

#   An alternative to group_seq_range
#   '''
#   selections = []

#   # Group by all columns except 'seq_id'
#   excluded = ["seq_id"]
#   grouped = df.groupby([col for col in df.columns if col not in excluded])  # Add more columns if needed

#   for name, group in grouped:
#     # Sort and reset index for convenience
#     group = group.sort_values('seq_id').reset_index(drop=True)

#     # Initialize variables for sequence detection
#     start_seq = end_seq = group.loc[0, 'seq_id']
#     data = {}

#     for i in range(1, len(group)):
#       curr_seq = group.loc[i, 'seq_id']

#       if curr_seq == end_seq + 1:
#         # Continue the sequence
#         end_seq = curr_seq
#       else:
#         # Break the sequence
#         if start_seq == end_seq:
#           seq_ops = [{'op': '==', 'value': int(start_seq)}]
#         else:
#           seq_ops = [{'op': '>=', 'value': int(start_seq)}, {'op': '<=', 'value': int(end_seq)}]

#         data = {
#           'asym_id': {'ops': [{'op': '==', 'value': name[0]}]},
#           'comp_id': {'ops': [{'op': '==', 'value': name[1]}]},
#           'atom_id': {'ops': [{'op': '==', 'value': name[2]}]},
#           'seq_id': {'ops': seq_ops},
#           # ... other columns?
#         }

#         selections.append(Selection(data, None))  # Add the Selection object
#         start_seq = end_seq = curr_seq  # Reset sequence variables

#     # Don't forget the last sequence
#     if not (start_seq == "*" and end_seq =="*"):
#       if start_seq == end_seq:
#         seq_ops = [{'op': '==', 'value': int(start_seq)}]
#       else:
#         seq_ops = [{'op': '>=', 'value': int(start_seq)}, {'op': '<=', 'value': int(end_seq)}]


#       data = {
#         'asym_id': {'ops': [{'op': '==', 'value': name[0]}]},
#         'comp_id': {'ops': [{'op': '==', 'value': name[1]}]},
#         'atom_id': {'ops': [{'op': '==', 'value': name[2]}]},
#         'seq_id': {'ops': seq_ops},
#         # ... other columns?
#       }
#     else:
#       data = {
#         'asym_id': {'ops': [{'op': '==', 'value': name[0]}]},
#         'comp_id': {'ops': [{'op': '==', 'value': name[1]}]},
#         'atom_id': {'ops': [{'op': '==', 'value': name[2]}]},
#         # ... other columns?
#       }


#     selections.append(Selection(data, None))

#   query = SelectionQuery(selections=selections)
#   return query

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

# def group_seq_range(simplest_nodes, columns=['asym_id', 'seq_id', 'atom_id']):
#   ''''
#   Take a selection dataframe and transform it into a new dataframe with each row containing a range
#   of seq_id rather than a single value. It uses the output of find_simplest_selected_nodes() as the start
#   '''

#   # Create DataFrame
#   df = pd.DataFrame(simplest_nodes, columns=columns)

#   # Initialize an empty list to store transformed rows
#   transformed_rows = []

#   # Generate groupby columns by removing 'seq_id' from the columns list
#   groupby_columns = [col for col in columns if col != 'seq_id']

#   # Group by non-sequential columns
#   for _, group in df.groupby(groupby_columns):
#     group = group.sort_values('seq_id')  # Sort by 'seq_id'

#     seq_id_start = group['seq_id'].iloc[0]
#     seq_id_stop = group['seq_id'].iloc[0]

#     # Initialize a row dictionary with other columns
#     row_dict = {col: group[col].iloc[0] for col in groupby_columns}

#     for i in range(1, len(group)):
#       current_seq_id = group['seq_id'].iloc[i]

#       # Check if the sequence is continuous
#       if current_seq_id == seq_id_stop + 1:
#         seq_id_stop = current_seq_id
#       else:
#         # Add current sequence to transformed list
#         row_dict['seq_id_start'] = seq_id_start
#         row_dict['seq_id_stop'] = seq_id_stop
#         transformed_rows.append(row_dict.copy())

#         # Reset the sequence
#         seq_id_start = current_seq_id
#         seq_id_stop = current_seq_id

#     # Add the final sequence to the transformed list
#     row_dict['seq_id_start'] = seq_id_start
#     row_dict['seq_id_stop'] = seq_id_stop
#     transformed_rows.append(row_dict.copy())

#   # Convert transformed rows into a DataFrame
#   transformed_df = pd.DataFrame(transformed_rows)

#   return transformed_df

def form_simple_str_common(df_sel_simple):

  final_component_strs = []
  for row in df_sel_simple.itertuples(index=False):
    s = f'(asym_id {row.asym_id} and seq_id {row.seq_id_start}:{row.seq_id_stop} and atom_id {row.atom_id})'
    final_component_strs.append(s)

  sel_str_simple_common = ' or '.join(final_component_strs)
  return sel_str_simple_common

# def form_simple_str_phenix(df_sel_simple,keyword_map):
#   sel_str_simple_common = form_simple_str_common(df_sel_simple)
#   keyword_map = {v:k for k,v in phenix_keyword_map.items()}
#   sel_str_simple_phenix = replace_keys_in_str(sel_str_simple_common,keyword_map)
#   return sel_str_simple_phenix