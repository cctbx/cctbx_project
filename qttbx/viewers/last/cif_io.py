from functools import reduce
from itertools import zip_longest
from collections import OrderedDict
import re

import pandas as pd
# Try import iotbx cif
try:
  from iotbx import cif
except:
  pass
# Try import pdbecif (pip install)
try:
  from pdbecif.mmcif_io import CifFileWriter, MMCIF2Dict
except:
  pass

"""
Tools for working with cif files in the context of pandas dataframes

The functions read_cif_file and write_cif_file do cif file io, and by default return nested
dictionaries with pandas dataframes as values. There are three options for the cif backend, iotbx cif, 
pdbecif, and a custom parser implemented here named "cifpd" which has only pandas for external dependencies.
"""
  
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

# A custom pure python cif parser with no external dependencies. Named cifpd

QUOTES_REGEX = re.compile(r'"[^"]*"|\'[^\']*\'|\S+')

def split_with_quotes_regex(line):
    parts = QUOTES_REGEX.findall(line)
    return [part[1:-1] if (part[0]=="'" or part[0]=='"') else part for part in parts]

split_with_quotes = split_with_quotes_regex


def check_for_semicolon(current_line_index, lines, used_semicolons):
    start = None
    end = None
    text_between_semicolons = []

    for i in range(current_line_index - 1, len(lines)):
        for j, char in enumerate(lines[i]):
            if char == ';' and (i, j) not in used_semicolons:
                if start is None:
                    start = (i, j)
                else:
                    end = (i, j)
                    break
        if start is not None:
            text_between_semicolons.append(lines[i])
        if end is not None:
            break

    if start and end:
        used_semicolons.extend([start, end])
        joined_text = ''.join(text_between_semicolons)
        # Remove the semicolons and anything before/after them in the start/end lines
        joined_text = joined_text[joined_text.find(';') + 1:joined_text.rfind(';')]
        return joined_text, used_semicolons
    else:
        return None, used_semicolons

      
def open_loop(lines, line_index):
    loop_keys = []
    loop_data = []
    in_loop = True
    line_index += 1  # Move to the next line after "loop_"
    used_semicolons = []
    
    while line_index < len(lines) and in_loop:
        line = lines[line_index]
        if len(line)>0:
          if line[0] =='_' and len(loop_data)==0:
              loop_keys.append(line)
          elif line and not line[0]=='_' and not line[0]=='#':
              parts = split_with_quotes(line)
              missing_parts = False
              if len(parts)>len(loop_keys): 
                assert False, f"Failed parsing line: {line}, to match number of keys: {len(loop_keys)}"
              elif len(parts)<len(loop_keys): 
                  # attempt to find missing parts
                  missing_parts = True
                  print("\nMissing parts at line:",line_index)
                  print(line)
                  print(len(parts),len(loop_keys))
                  print("\n")
            
                  s, used_semicolons = check_for_semicolon(line_index,lines,used_semicolons)
                  if s is not None:
                    parts.append(s)
                    line_index = used_semicolons[-1][0]
                  if len(parts)==len(loop_keys):
                    missing_parts = False
                    print("Found missing parts with semicolon search")
                  elif len(parts)>len(loop_keys): 
                    assert False, f"Failed parsing line: {line}, to match number of keys: {len(loop_keys)}"
                  elif len(parts)<len(loop_keys):  
                    # still not enough,  now check for non-semicolon next line continuation (only check one)
                    print("Checking next line continuation")
                    line = lines[line_index]
                    line_next = lines[line_index+1]
                    parts_next = split_with_quotes_regex(line_next)
                    if len(parts)+len(parts_next)==len(loop_keys):
                      # accept that the next line is a continuation
                      missing_parts = False
                      print("Found missing parts with next-line continuation")
                      parts = parts+parts_next
                      line_index += 1
                    else:
                      # A failure to assign values to keys
                      assert len(parts)==len(loop_keys), f"Unable to fix mismatch between number of values in loop row: {parts},  and number of keys({len(loop_keys)}): {loop_keys}"
              loop_data.append(parts)
          else:  # End of loop data
              in_loop = False
              line_index -= 1
        line_index += 1

    return loop_keys, loop_data, line_index

def add_single_entry(line,current_data_key,result):
  parts = split_with_quotes(line)
  group_key, column_key = parts[0].split('.')[0], parts[0].split('.')[-1]
  value = ' '.join(parts[1:])
  if group_key not in result[current_data_key]:
      result[current_data_key][group_key] = {}
  result[current_data_key][group_key][column_key] = value
  last_key = (current_data_key,group_key,column_key)
  return last_key
  
def close_loop(loop_keys, loop_data):
    result = {}
    group_key = loop_keys[0].split('.')[0] if loop_keys else ""
    result[group_key] = {}
    for i, key in enumerate(loop_keys):
        column_key = key.split('.')[-1]
        column_data = [data[i] for data in loop_data if i < len(data)]
        result[group_key][column_key] = column_data
    return result

def parse_cifpd(content):
    """
    The entry function to parse cif string contents. 
    """
    lines = [line.strip() for line in content.splitlines()]
    result = {}
    current_data_key = None
    last_key = None
    line_index = 0
    multi_line_single_start = None # flag for semicolon, denoting multi line text

    while line_index < len(lines):
        line = lines[line_index]
        if not line or line[0]=='#':
            line_index += 1
            continue  # Skip empty and comment lines

        if multi_line_single_start is not None:
          if line[0]!=";":
            line_index += 1
            continue
          else:
            # close multi line text
            multi_line_single_end = line_index
            current_data_key, group_key, column_key = last_key
            result[current_data_key][group_key][column_key] = "".join(lines[multi_line_single_start:multi_line_single_end]).replace(";","")
            multi_line_single_start = None
            line_index += 1
            continue
        if line[0]==";":
          multi_line_single_start = line_index
          line_index += 1
          continue
      
        if line[:5] == "data_":
            current_data_key = line.replace("data_","")
            result[current_data_key] = {}
        elif line == "loop_":
            loop_keys, loop_data, line_index = open_loop(lines, line_index)
            loop_result = close_loop(loop_keys, loop_data)
            result[current_data_key].update(loop_result)
            continue  # open_loop and close_loop handle line_index increment
        elif line[0]=="_" and current_data_key:
            last_key = add_single_entry(line,current_data_key,result)
        else:
          # maybe line is value for previous key on a new line
          current_data_key, group_key, column_key = last_key
          
          last_value= result[current_data_key][group_key][column_key]
          if last_value != "":
            print("IGNORED LINE: ",line)
          else:
            # assume value for previous key on a new line
            result[current_data_key][group_key][column_key] = split_with_quotes_regex(line)[0]
          
        line_index += 1

    return result

# Format and write cif text directly from dataframes

def infer_type(element):
    try:
        # Convert to float first
        float_val = float(element)
        # If the float is actually an integer, convert to int
        if float_val.is_integer():
            return int(float_val)
        return float
    except ValueError:
        # If conversion fails, it's a string
        return str

def convert_column_types_based_on_first(df):
    # Get the first row of the DataFrame
    sample_row = df.iloc[0]

    # Infer the type for each column based on the first row
    inferred_types = {col: infer_type(sample_row[col]) for col in df.columns}

    # Convert entire columns to the inferred types
    for col, dtype in inferred_types.items():
        if dtype is float or dtype is int:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        # Else, it's a string, no need to convert

    return df


  
def quote_strings_with_spaces(df):
  # Iterate over columns
  for col in df.columns:
      if df[col].dtype == 'object':  # Check if column is string type
          # Add quotes to strings with spaces
          df[col] = df[col].apply(lambda x: f'"{x}"' if isinstance(x, str) and ' ' in x else x)
  return df
  
def format_dataframe_for_cif(df):
    # Determine padding for each column
    col_format = {}
    for col in df.columns:
        if df[col].dtype == float or df[col].apply(lambda x: isinstance(x, float)).any():
            max_left_len = df[col].apply(lambda x: len(str(int(x)) if isinstance(x, float) and not pd.isna(x) else str(x))).max()
            max_right_len = df[col].apply(lambda x: len(str(x).split('.')[1]) if isinstance(x, float) and '.' in str(x) else 0).max()
            total_width = max_left_len + max_right_len + 1  # +1 for decimal
            df[col] = df[col].apply(lambda x: f"{x:>{total_width}.{max_right_len}f}" if isinstance(x, float) else str(x).ljust(total_width))
        else:
            max_len = df[col].astype(str).str.len().max()
            df[col] = df[col].astype(str).str.ljust(max_len)

    # Concatenate the formatted columns
    formatted_df = df.apply(lambda x: ' '.join(x), axis=1)

    return '\n'.join(formatted_df)

# Example DataFrame (replace with your actual Da


def df_to_cif_lines(df,decimal_places=3,integer_padding=4,decimal_padding=4,column_prefix=None,data_name=None):
  
  # replace python uncertain values with cif ones
  df.replace(to_replace=[None,""," "], value='.', inplace=True)
  
  if len(df)>1: #a loop
    lines_header = ["loop_"]
    for column in df.columns:
      if column_prefix is not None:
        column = column_prefix+"."+column
      if not column.startswith("_"):
        column = "_"+column
      lines_header.append(column)

    lines_body = format_dataframe_for_cif(df).splitlines()
    lines = lines_header+lines_body
    
  else: # not a loop
    kv_list = list(df.to_dict("list").items())
    column_prefix_keys = [".".join([column_prefix,key]) for key,values in kv_list]
    max_len = max([len(key) for key in column_prefix_keys])
    
    lines = []
    for i,(key,values) in enumerate(kv_list):
      assert len(values)==1
      value = values[0]
      column_prefix_key = column_prefix_keys[i]
      line = column_prefix_key.ljust(max_len+2)+str(value).ljust(max_len)
      lines.append(line)

  if data_name is not None:
    lines = ["data_"+data_name]+lines
  return lines

def write_dataframes_to_cif_file(dataframe_dict,filename):
  lines_out = []
  for data_key,data_value in dataframe_dict.items():
    lines_out.append("data_"+data_key)
    for group_key,group_value in data_value.items():
      if isinstance(group_value,pd.DataFrame):
        df = convert_column_types_based_on_first(group_value)
        df = quote_strings_with_spaces(df)
        #df = group_value
  
        #df = group_value.applymap(str)
        lines = df_to_cif_lines(df,column_prefix=group_key)
        lines_out+=lines
        lines_out.append("#")
      else:
        assert False, f"Encountered non-dataframe value: {type(group_value)}"
  s = "\n".join(lines_out)
  with open(filename,"w") as fh:
    fh.write(s)


# Finally, Actual IO functions

def read_cif_file(filename,method="pdbe",return_as="pandas"):
  assert method in ["pdbe","iotbx","cifpd"]
  assert return_as in ["dict","pandas","iotbx"]
  if method == "iotbx":
    if isinstance(filename,str):
      reader = cif.reader(file_path=filename)
    else:
      reader = cif.reader(file_obj=filename)
    model = reader.model()
    if return_as == "iotbx":
      return model

    d = convert_iotbx_cif_to_dict(model)
    if return_as == "dict":
      return d
      
    if return_as == "pandas":
      dfs = convert_dict_to_dataframes(d)
    return dfs
    
  elif method == "pdbe":
    d = MMCIF2Dict().parse(str(filename))
    d = clean_nested_dict(d)

    if return_as == "iotbx":
      raise NotImplementedError
    if return_as == "dict":
      return d

    if return_as == "pandas":
      dfs = convert_dict_to_dataframes(d)
      return dfs

  elif method == "cifpd":
    with open(filename,"r") as fh:
      content = fh.read()
      
    d = parse_cifpd(content)
    if return_as == "dict":
      return d
    if return_as == "pandas":
      dfs = convert_dict_to_dataframes(d)
      return dfs
    if return_as == "iotbx":
      raise NotImplementedError


def write_cif_file(inp,filename,inp_type="pandas",method="pdbe"):
  assert method in ["pdbe","iotbx", "cifpd"]
  assert inp_type in ["dict","pandas","iotbx"]
  if method == "pdbe":
    if inp_type== "pandas":
      d = convert_dataframes_to_dict(inp)
    elif inp_type == "dict":
      d = object
    elif inp_type == "iotbx":
      raise NotImplementedError
      
    CifFileWriter(filename).write(d)
    
  elif method == "iotbx":
    if inp_type == "iotbx":
      model = inp
    elif inp_type == "dict":
      dfs = convert_dict_to_dataframes(inp)
      model = convert_dataframes_to_iotbx_cif(dfs)
    elif inp_type == "pandas":
      model = convert_dataframes_to_iotbx_cif(inp)

    with open(filename,"w") as fh:
      model.show(out=fh)

  elif method == "cifpd":
    if inp_type == "pandas":
      dfs = inp
    elif inp_type == "dict":
      dfs = convert_dict_to_dataframes(inp)

    elif inp_type == "iotbx":
      d = convert_iotbx_cif_to_dict(inp)
      dfs = convert_dict_to_dataframes(d)

    write_dataframes_to_cif_file(dfs,filename)
      
# Testing tools

def compare_nested_dataframe_dicts(dict1, dict2, parent_key=''):
  """
  Take two nested dicts with pd.DataFrame leaves, compare the contents.
  """
  if dict1.keys() != dict2.keys():
      raise ValueError(f"Keys mismatch at {parent_key}: {dict1.keys()} vs {dict2.keys()}")

  for key in dict1.keys():
      current_key = f"{parent_key}.{key}" if parent_key else key
      if isinstance(dict1[key], pd.DataFrame) and isinstance(dict2[key], pd.DataFrame):
          if not dict1[key].equals(dict2[key]):
              raise ValueError(f"DataFrames mismatch at {current_key}")
      elif isinstance(dict1[key], dict) and isinstance(dict2[key], dict):
          compare_nested_dataframe_dicts(dict1[key], dict2[key], current_key)
      else:
          if dict1[key] != dict2[key]:
              raise ValueError(f"Values mismatch at {current_key}: {dict1[key]} vs {dict2[key]}")

  return True