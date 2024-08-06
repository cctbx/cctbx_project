

from io import StringIO
import warnings
from collections import UserDict, defaultdict
import re
from pathlib import Path

import pandas as pd
from pandas.errors import ParserWarning

from iotbx import cif
from iotbx.cif.builders import crystal_symmetry_builder
from iotbx.cif.model import block

from .python_utils import find_key_path, get_value_by_path


# Conversion functions
def convert_dict_to_dataframes(d):
  """
  Convert a nested dictionary with list leaves to nested dataframes
  Types are inferred from the contents of the lists
  """
  if isinstance(d, dict):
      # Check if all values are lists of the same length
      if all(isinstance(v, list) for v in d.values()) and len(set(len(v) for v in d.values())) == 1:
          return pd.DataFrame(d)
      elif all(not isinstance(child, dict) for child in d.values()):
          # Handle single-row dataframes
          return pd.DataFrame({k: [v] for k, v in d.items()})
      else:
          # Recursively apply the conversion to each child dictionary
          return {k: convert_dict_to_dataframes(v) for k, v in d.items()}
  else:
      return d  # In case the value is not a dictionary

def convert_dataframes_to_dict(obj):
  """
  Convert a nested dict with pd.DataFrame leaves to
  nested dict with list leaves

  """
  if isinstance(obj, pd.DataFrame):
      # Convert DataFrame to a dictionary of lists
      return obj.to_dict(orient='list')
  elif isinstance(obj, dict):
      # Recursively process dictionary
      return {k: convert_dataframes_to_dict(v) for k, v in obj.items()}
  else:
      # Return the object as-is if it's neither a DataFrame nor a dictionary
      return obj

def convert_iotbx_cif_to_dict(cif_model):
  d = remove_iotbx_cif(cif_model)
  d = nest_dict(d)
  d = clean_nested_dict(d)
  return d

def convert_dataframes_to_iotbx_cif(dfs):
  """
  Convert nested dataframes to a cctbx cif model
  Dataframes are converted to string here.
  """
  cif_out = cif.model.cif()
  for data_key,data_d in dfs.items():
    data_block = cif.model.block()
    for block_key, df in data_d.items():
      block = cif.model.block()
      assert isinstance(df,pd.DataFrame)

      # Convert all elements to strings
      df_as_strings = df.applymap(str)

      # Convert DataFrame to a list of lists (columns as lists of strings)
      list_of_columns = [df_as_strings[col].tolist() for col in df_as_strings.columns]

      if len(df)>1:
        list_as_cctbx = []
        for l in list_of_columns:
          #print("l",l)
          list_as_cctbx.append(l)
        loop = cif.model.loop()
        for i,column in enumerate(df.columns):

          k = block_key+"."+column
          d = list_as_cctbx[i]
          loop.add_column(k,d)
        block.add_loop(loop)
      else:
        data = [l[0] for l in list_of_columns]
        for column,value in zip(df.columns,data):
          c = block_key+"."+column
          block.add_data_item(c,value)

      data_block.update(block)
    cif_out[data_key] = data_block
  return cif_out

# General utilities for dicts, strings, and dataframes

def find_key_in_dict(d, target_key, counter=0, value=None):
  # search for a single key (like _atom_site) in a dict
  if not isinstance(d, dict):
    return counter, value

  for key, val in d.items():
    if key == target_key:
      counter += 1
      value = val
    if isinstance(val, dict):
      counter, value = find_key_in_dict(val, target_key, counter, value)
    elif isinstance(val, list):
      for item in val:
        counter, value = find_key_in_dict(item, target_key, counter, value)

  return counter, value



def remove_iotbx_cif(d, condition=lambda x: isinstance(x, int) and x > 10, new_value=0):
    new_dict = {}
    for key, value in d.items():
        if isinstance(value, dict):
            # Recursively call the function for nested dictionaries
            new_dict[key] = remove_iotbx_cif(value, condition, new_value)
        else:
            if "iotbx.cif" in str(type(value)):
              new_dict[key] = dict(value)
            else:
              new_dict[key] = value
    return new_dict

def nest_dict(flat_cif_dict):
  # data keys on top, block and value keys combined with "."
  for data_key,value in list(flat_cif_dict.items()):
    flat_cif_dict[data_key] = unpack_dict(value)
  return flat_cif_dict

def unpack_dict(d):
    def expand_key(key, value):
        keys = key.split('.')
        current = result = {}
        for k in keys[:-1]:
            current[k] = {}
            current = current[k]
        current[keys[-1]] = value
        return result

    if not isinstance(d, dict):
        return d

    result = {}
    for key, value in d.items():
        if '.' in key:
            expanded = expand_key(key, unpack_dict(value))
            for k, v in expanded.items():
                result.setdefault(k, {}).update(v)
        else:
            result[key] = unpack_dict(value)
    return result

def flatten_dict(d, parent_key='', sep='.'):
    items = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

def no_single_element_list(d):
    for key, value in d.items():
        if isinstance(value, dict):
            no_single_element_list(value)  # Recursive call for nested dictionaries
        elif isinstance(value, list) and len(value) == 1:
            d[key] = value[0]  # Replace the list with its first element


def clean_string(s):
    return s.strip().replace('\n', '').replace("\t","")

def clean_nested_dict(d):
    for key, value in list(d.items()):
        if isinstance(value, str):
            d[key] = clean_string(value)
        elif isinstance(value, list) or "scitbx" in str(type(value)):
            d[key] = [clean_string(item) if isinstance(item, str) else item for item in value]
        elif isinstance(value, dict):
            clean_nested_dict(value)
    return d


from io import StringIO
import warnings
from collections import UserDict, defaultdict
import re
from pathlib import Path

from iotbx import cif
from iotbx.cif.builders import crystal_symmetry_builder
from iotbx.cif.model import block
import pandas as pd
from pandas.errors import ParserWarning

from .python_utils import find_key_path, get_value_by_path
from .cif_io import convert_dict_to_dataframes



class CifInput(UserDict):

  # Define explicit dtypes, what isn't here become 'string'
  dtypes = defaultdict(lambda: 'string')
  dtypes.update({
  'auth_asym_id': 'string',
  'label_asym_id': 'string',
  'auth_seq_id': 'Int64',
  'label_seq_id': 'Int64',
  'label_comp_id': 'string',
  'type_symbol': 'string',
  'label_alt_id': 'string',
  'pdbx_PDB_model_num': 'string',
  'pdbx_PDB_ins_code': 'string',
  'Cartn_x': 'Float64',
  'Cartn_y': 'Float64',
  'Cartn_z': 'Float64',
  'B_iso_or_equiv': 'Float64',
  'occupancy': 'Float64',
  'pdbx_formal_charge': 'Int64',
  'group_PDB': 'string',
  'id': 'string',
  'label_atom_id': 'string',
  'label_entity_id': 'string'
  })

  @classmethod
  def from_mmcif_file(cls,filename):
    return cls(filename)

  @classmethod
  def from_mmtbx_model_via_mmcif(cls,model):
    # use mmcif string as intermediate
    string_obj = StringIO(model.model_as_mmcif())
    return cls(string_obj)


  def __init__(self,*args,method="iotbx",**kwargs):
    """
    Initialize with regular dictionary args,
      or a filename of a cif file
    """

    # Check if initialized with a filename
    d = None
    self.filename = None
    if len(args) == 1:
      if isinstance(args[0], (Path,str)):
        self.filename = args[0]
        d = self._read_mmcif_file_to_dataframes(self.filename,method=method)
        args = args[1:]
      elif isinstance(args[0],StringIO):
        d = self._read_mmcif_obj_to_dataframes(args[0],method=method)
        args = args[1:]

    super().__init__(*args,**kwargs)
    if d is not None:
      self.update(d)



  def read_mmcif_file_to_dataframes(self,file,method=None):
    d = self._read_mmcif_file_to_dataframes(file,method=method)
    self.update(d)

  def _read_mmcif_file_to_dataframes(self,file,method=None):
    with open(file, 'r') as fh:
      return self._read_mmcif_obj_to_dataframes(fh,method=method)

  def read_mmcif_obj_to_dataframes(self,file_obj,method=None):
    d = self._read_mmcif_obj_to_dataframes(file_obj,method=method)
    self.update(d)

  def _read_mmcif_obj_to_dataframes(self,file_obj,method=None):
    assert method in ["iotbx","cifpd"]
    if method == "cifpd":
      assert False, "Not fully functional"
      lines = file_obj.readlines()

      self.parser= CifParser(lines)
      self.parser.run()
      dataframes = self.parser.process_to_pandas()
      return dataframes

    elif method == "iotbx":
      reader = cif.reader(file_object=file_obj)
      cif_model = reader.model()
      dicts = convert_iotbx_cif_to_dict(cif_model)
      for key in list(dicts.keys()):
        if not key.startswith("data_"):
          dicts["data_"+key] = dicts.pop(key)
      dfs = convert_dict_to_dataframes(dicts)
      return dfs
  @property
  def atom_sites(self):
    # Search rather than assume location
    atom_site_path = find_key_path(self,"_atom_site")
    if atom_site_path is None:
      return None
    df =  get_value_by_path(self,atom_site_path)
    return df

  def extract_crystal_symmetry(self):
    """
    Extract the crystal symmetry as an object from a nested cif dict
    """

    # search the nested dict and retrieve value
    key_path = find_key_path(self,"_symmetry")
    d_symmetry = get_value_by_path(self,key_path)

    # value is dataframe, convert to dict
    d_cs =d_symmetry.to_dict(orient="records")[0]

    # join cif keys with dot
    d_cs = {"_symmetry."+key:value for key,value in d_cs.items()}

    # build the cctbx cif block
    b = block()
    for k,v in d_cs.items():
      b[k] = v

    # get crystal symmetry
    cs_builder = crystal_symmetry_builder(b)
    cs = cs_builder.crystal_symmetry
    return cs

  def find_dataframe(self,name):
    key_path = find_key_path(self,name)
    if key_path is None:
      # try with prefix
      name = "_"+name
      key_path = find_key_path(self,name)

    if key_path is None:
      return None
    d = get_value_by_path(self,key_path)
    return d
import pandas as pd
import re

class CifParser:
  def __init__(self, lines):
    self.lines = lines
    self.data_idxs = [i for i, line in enumerate(lines) if line.startswith("data_")]
    self.data_idxs.append(len(lines))

  def run(self):
    self.results = {}
    for i, start in enumerate(self.data_idxs[:-1]):
      block_name = self.lines[start].strip()
      end = self.data_idxs[i + 1]
      block_lines = self.lines[start:end]
      parser = CifBlockParser(block_lines)
      parser.run()
      dfs = parser.process_to_pandas()
      self.results[block_name] = dfs

  def process_to_pandas(self):
    return self.results

class CifBlockParser:
  def __init__(self, lines):
    self.lines = lines
    self.i = 0
    self.columns = []
    self.data = []
    self.results = {}
    self.loop_parsing = False
    self.block_parsing = False
    self.current_prefix = None
    self.ready_to_eat_loop = False
    self.ready_to_eat_block = False
    self.ready_to_eat_key_value = False
    self.loop_marker = "loop_"
    self.last_result_key = None
    self.current_result = {}
    self.in_loop = False

  @property
  def line(self):
    return self.lines[self.i]

  def run(self):
    while self.i < len(self.lines):
      self.read_next_line()
      self.end_line()
    self.finish_file()

  def clear(self):
    self.columns = []
    self.data = []
    self.loop_parsing = False
    self.block_parsing = False
    self.ready_to_eat_loop = False
    self.ready_to_eat_block = False
    self.ready_to_eat_key_value = False

  def eat(self, is_loop=False, is_block=False):
    key = self.current_prefix if is_loop or is_block else self.columns[0]
    self.current_result = {
      "columns": self.columns,
      "data": self.data,
      "is_loop": is_loop,
      "is_block": is_block,
      "is_single_value": not (is_loop or is_block),
      "prefix": self.current_prefix if is_loop or is_block else None,
      "stop": self.i
    }
    self.results[key] = self.current_result
    self.last_result_key = key
    self.current_result = {"start": self.i + 1}
    self.clear()

  def end_line(self):
    if self.ready_to_eat_key_value:
      self.eat(is_loop=False, is_block=False)
    elif self.ready_to_eat_loop:
      self.eat(is_loop=True, is_block=False)
    elif self.ready_to_eat_block:
      self.eat(is_loop=False, is_block=True)
    self.i += 1

  def finish_file(self):
    if self.loop_parsing:
      self.ready_to_eat_loop = True
    elif self.block_parsing:
      self.ready_to_eat_block = True
    elif self.columns:
      self.ready_to_eat_key_value = True
    self.end_line()

  def process_block_line(self, line):
    if line.startswith("_"):
      prefix, key, value = self.process_single_line_data(line)
      if prefix is not None:
        key = f"{prefix}.{key}"
      self.current_prefix = prefix
      self.columns.append(key)
      self.data.append(value)
    else:
      self.current_result["is_strange"] = True
      if self.data and self.data[-1] == "":
        self.data[-1] = line

  def start_block(self, prefix, key, value):
    self.current_prefix = prefix
    self.columns.append(key)
    self.data.append(value)
    self.block_parsing = True

  def start_loop(self):
    self.current_prefix = None
    self.loop_parsing = True
    self.in_loop = True

  def process_single_line_data(self, line):
    line_split = split_with_quotes(line)
    if len(line_split) == 1:
      line_split.append("")
    assert len(line_split) == 2, f"Failed to split key,value pair: {line_split}"
    prefix = None
    key, value = line_split
    prefix_split = key.split(".")
    if len(prefix_split) > 1:
      assert len(prefix_split) == 2
      prefix, key = prefix_split
    return prefix, key, value

  def read_next_line(self):
    line = self.lines[self.i]
    stripped_line = line.strip()
    is_key_row = stripped_line.startswith('_')
    is_break_row = not stripped_line or stripped_line.startswith("#")
    is_loop_entry_row = stripped_line.startswith(self.loop_marker)
    is_multi_line = stripped_line.startswith(";")

    if is_multi_line:
      multiline_value, last_index = handle_multiline(self.lines, self.i)
      if not self.columns:
        self.results[self.last_result_key]["data"][-1] = multiline_value
      else:
        self.data[-1] = multiline_value
      self.i = last_index
      return

    is_no_sep_break = False
    if self.loop_parsing or self.block_parsing:
      if not is_break_row and (is_key_row or is_loop_entry_row):
        if is_loop_entry_row:
          is_no_sep_break = True
        elif self.current_prefix and not stripped_line.startswith(self.current_prefix):
          is_no_sep_break = True

    if self.block_parsing:
      if is_break_row or is_no_sep_break:
        self.ready_to_eat_block = True
        return
      else:
        self.process_block_line(stripped_line)
        return

    if self.loop_parsing:
      if is_break_row or is_no_sep_break:
        self.ready_to_eat_loop = True
        return
      elif is_key_row:
        self.columns.append(stripped_line)
        if not self.current_prefix:
          prefix, _, _ = self.process_single_line_data(stripped_line)
          self.current_prefix = prefix
        return
      else:
        if self.data and isinstance(self.data[-1], str) and self.data[-1].startswith(";"):
          multiline_value, last_index = handle_multiline(self.lines, self.i)
          self.data[-1] = multiline_value
          self.i = last_index
        else:
          self.data.append(self.process_loop_data(stripped_line))
        return

    if is_loop_entry_row:
      self.start_loop()
      return

    if is_key_row:
      prefix, key, value = self.process_single_line_data(stripped_line)
      if prefix is not None:
        self.start_block(prefix, key, value)
        return
      else:
        self.columns.append(key)
        self.data.append(value)
        self.ready_to_eat_key_value = True
        return

    self.current_result["is_strange"] = True

  def process_to_pandas(self):
    output_dfs = {}
    for key, result in self.results.items():
      if not isinstance(result, dict):
        output_dfs[key] = result
        continue

      columns = result["columns"]
      if result["is_loop"]:
        data = result["data"]
        if isinstance(data[0], int):
          data = self.lines[data[0]:data[-1] + 1]
        else:
          data = [" ".join(row) for row in data]
        if len(data) == 1 and data[0].startswith(";"):
          data = self.process_multiline_loop_data(data[0])
        df = pd.DataFrame([data], columns=columns)
        if result["prefix"]:
          df.columns = [col.replace(result["prefix"] + ".", "") for col in df.columns]
        output_dfs[key] = df
      elif result["is_block"]:
        data = result["data"]
        if len(data) == len(columns):
          data = [data]
        df = pd.DataFrame(data, columns=columns)
        if result["prefix"]:
          df.columns = [col.replace(result["prefix"] + ".", "") for col in df.columns]
        output_dfs[key] = df
      elif result["is_single_value"]:
        output_dfs[key] = result["value"]

    return output_dfs

  def process_multiline_loop_data(self, data):
    data = data.strip().strip(";").split("\n")
    processed_data = []
    for line in data:
      line = line.strip()
      if line:
        processed_data.append(line)
    return [" ".join(processed_data)]

  def process_loop_data(self, line):
    if line.startswith("'") or line.startswith('"'):
      return [line.strip()]
    return line.split()

QUOTES_REGEX = re.compile(r'"[^"]*"|\'[^\']*\'|\S+')

def split_with_quotes(line):
  parts = QUOTES_REGEX.findall(line)
  return [part[1:-1] if part[0] in ("'", '"') else part for part in parts]

def handle_multiline(lines, current_index):
  multiline_content = []
  for j in range(current_index, len(lines)):
    line = lines[j]
    multiline_content.append(line)
    if line.strip() == ";":
      break
  return "\n".join(multiline_content), j




## Writing
###########

# Format and write cif text directly from dataframes

def infer_type(element):
  try:
    float_val = float(element)
    if float_val.is_integer():
      return int
    return float
  except ValueError:
    return str

def convert_column_types_based_on_first(df):
  sample_row = df.iloc[0]
  inferred_types = {col: infer_type(sample_row[col]) for col in df.columns}
  for col, dtype in inferred_types.items():
    if dtype in {float, int}:
      df[col] = pd.to_numeric(df[col], errors='coerce')
  return df

def quote_strings_with_spaces(df):
  for col in df.columns:
    if df[col].dtype == 'object':
      df[col] = df[col].apply(lambda x: f'"{x}"' if isinstance(x, str) and ' ' in x else x)
  return df

def format_dataframe_for_cif(df):
  for col in df.columns:
    if df[col].dtype == float or df[col].apply(lambda x: isinstance(x, float)).any():
      max_left_len = df[col].apply(lambda x: len(str(int(x))) if isinstance(x, float) and not pd.isna(x) else len(str(x))).max()
      max_right_len = df[col].apply(lambda x: len(str(x).split('.')[1]) if isinstance(x, float) and '.' in str(x) else 0).max()
      total_width = max_left_len + max_right_len + 1
      df[col] = df[col].apply(lambda x: f"{x:>{total_width}.{max_right_len}f}" if isinstance(x, float) else str(x).ljust(total_width))
    else:
      max_len = df[col].astype(str).str.len().max()
      df[col] = df[col].astype(str).str.ljust(max_len)
  formatted_df = df.apply(lambda x: ' '.join(x), axis=1)
  return '\n'.join(formatted_df)

def df_to_cif_lines(df, column_prefix=None, data_name=None):
  df.replace(to_replace=[None, "", " "], value='.', inplace=True)
  if len(df) > 1:
    lines_header = ["loop_"]
    for column in df.columns:
      column_name = column_prefix + "." + column if column_prefix else column
      if not column_name.startswith("_"):
        column_name = "_" + column_name
      lines_header.append(column_name)
    lines_body = format_dataframe_for_cif(df).splitlines()
    lines = lines_header + lines_body
  else:
    kv_list = list(df.to_dict("list").items())
    column_prefix_keys = [".".join([column_prefix, key]) for key, _ in kv_list]
    max_len = max(len(key) for key in column_prefix_keys)
    lines = []
    for i, (key, values) in enumerate(kv_list):
      assert len(values) == 1
      value = values[0]
      column_prefix_key = column_prefix_keys[i]
      line = column_prefix_key.ljust(max_len + 2) + str(value)
      lines.append(line)
  if data_name:
    if not data_name.startswith("data_"):
      data_name = "data_"+ data_name 
    lines = [data_name] + lines
  return lines

def write_dataframes_to_cif_file(dataframe_dict, filename):
  lines_out = []
  for data_key, data_value in dataframe_dict.items():
    if not data_key.startswith("data_"):
      data_key = "data_" + data_key
    lines_out.append(data_key)
    for group_key, group_value in data_value.items():
      if isinstance(group_value, pd.DataFrame):
        df = convert_column_types_based_on_first(group_value)
        df = quote_strings_with_spaces(df)
        lines = df_to_cif_lines(df, column_prefix=group_key)
        lines_out += lines
        lines_out.append("#")
      else:
        raise ValueError(f"Encountered non-dataframe value: {type(group_value)}")
  s = "\n".join(lines_out)
  with open(filename, "w") as fh:
    fh.write(s)