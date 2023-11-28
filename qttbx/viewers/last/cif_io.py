import sys
import io
from pathlib import Path
from functools import reduce
from itertools import zip_longest

import numpy as np
import pandas as pd

from iotbx import cif
try:
  from pdbecif.mmcif_io import CifFileWriter, MMCIF2Dict
except:
  pass


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




def nest_dict(flat_dict):
    # take a flat dictionary with implicit nesting 
    # using dot syntax (like iotbx.cif) and make nesting explicit
    # It assumes a single head key
    head_keys = list(flat_dict.keys())
    assert len(head_keys)==1
    head_key = head_keys[0]
    flat_dict = flat_dict[head_key]
    nested_dict = {}
    for key, value in flat_dict.items():
        parts = key.split('.')
        d = nested_dict
        for part in parts[:-1]:
            if part not in d:
                d[part] = {}
            d = d[part]
        if not isinstance(value,str):
            value = list(value)
        d[parts[-1]] = value
    return {head_key:nested_dict}

def convert_to_dataframes(d):
    if isinstance(d, dict):
        # Check if all values are lists of the same length
        if all(isinstance(v, list) for v in d.values()) and len(set(len(v) for v in d.values())) == 1:
            return pd.DataFrame(d)
        elif all(not isinstance(child, dict) for child in d.values()):
            # Handle single-row dataframes
            return pd.DataFrame({k: [v] for k, v in d.items()})
        else:
            # Recursively apply the conversion to each child dictionary
            return {k: convert_to_dataframes(v) for k, v in d.items()}
    else:
        return d  # In case the value is not a dictionary

def max_dict_depth(d):
    if isinstance(d, dict):
        return 1 + (max(map(max_dict_depth, d.values())) if d else 0)
    return 0


def find_paths(d, path=None):
    paths = []
    for k, v in d.items():
        new_path = path + [k] if path else [k]
        if isinstance(v, dict):
            paths += find_paths(v, new_path)
        else:
            paths.append(new_path)
    return paths
  
def resolve_value(d, keys):
    for key in keys:
        if isinstance(d, dict) and key in d:
            d = d[key]
        else:
            raise KeyError('Invalid key: {}'.format(key))
    return d
  
  
def merge_dicts(dict1, dict2):
    # Recursively merge the dictionaries
    merged = {}
    for key in set(dict1.keys()).union(dict2.keys()):
        if key in dict1 and key in dict2:
            if isinstance(dict1[key], dict) and isinstance(dict2[key], dict):
                merged[key] = merge_dicts(dict1[key], dict2[key])  # Recursive call for nested dictionaries
            else:
                merged[key] = dict2[key]  # Override value in dict1 with value from dict2
        elif key in dict1:
            merged[key] = dict1[key]
        else:
            merged[key] = dict2[key]
    return merged
  
def set_value_at_nested_key(dictionary, keys, value):
    # Define the helper function to update the nested value
    def update_nested_value(d, key):
        if key not in d:
            d[key] = {}
        return d[key]

    # Use reduce to traverse the dictionary hierarchy and set the value
    reduce(update_nested_value, keys[:-1], dictionary)[keys[-1]] = value
    
    
def read_cif_file(fileinp,method="pdbe"):
  assert method in ["lastcif","pdbe","iotbx"]
  if method == "iotbx":
    if isinstance(fileinp,str):
      reader = cif.reader(file_path=fileinp)
    else:
      reader = cif.reader(file_obj=fileinp)
    model = reader.model()
    d = nest_dict(model)
    return d
  elif method == "pdbe":
    return MMCIF2Dict().parse(str(fileinp))
  elif method =="lastcif":
    if isinstance(fileinp,io.StringIO):
      lines = fileinp.readlines()

      return read_cif_lines(lines)
    filepath = Path(fileinp)
    if ".gz" in filepath.suffixes:
      with gzip.open(str(filepath),"rt") as fh:
        lines = fh.readlines()
        return read_cif_lines(lines)
    else:
      with open(str(fileinp),"r") as fh:
        lines = fh.readlines()

        return read_cif_lines(lines)
    
def split_and_quote(line):
  """
  This is a method to preserve quotes when splitting string by whitespace.
  It uses shlex and is slower than the other manual implementation
  """
  if "'" in line or '"' in line:
    split_line = shlex.split(line)
    for i, item in enumerate(split_line):
        if "'" + item + "'" in line or '"' + item + '"' in line:
            split_line[i] = "'" + item + "'"
    return split_line
  else:
    return shlex.split(line)
  
def split_with_quotes(line):
  """
  This is a method to preserve quotes when splitting string by whitespace.
  If no quotes are present in the line, it just returns line.split()
  """
  s = line
  if "'" in line or '"' in line:
    parts = []
    current_word = ""
    inside_quotes = False

    for c in s:
        if c == ' ' and not inside_quotes:
            if current_word:
                parts.append(current_word)
                current_word = ""
        elif c == "'":
            inside_quotes = not inside_quotes
        else:
            current_word += c

    if current_word:
        parts.append(current_word)
    return parts
  else:
    
    return line.split()

def read_cif_lines(lines):


  current_data_key = None
  loops = []
  failures = []
  in_loop = False
  loop_start = None
  loop_data_start = None
  loop_keys = []
  lines_len = len(lines)-1
  building_d = {}
  debug_stop = False
  for i,line in enumerate(lines):
      if debug_stop:
        continue
      
      
      line = line.replace("\n","").lstrip()
      whitespace_index = line.find(" ")
      if len(line)>0 and line[:4]=="data":
        data_key = line.split("data_")[-1].strip()
        building_d[data_key] = {}
        current_data_key = data_key
      elif len(line)>0 and line[0] == "_":
        if whitespace_index>-1:
          keys = line[:whitespace_index]
          value = line[whitespace_index:]
        else:
          keys = line
          value = ""
         
        value = value.strip()
        # if len(value)>2:
        #   if value[0]=="'" and value[-1]=="'":
        #     value = value[1:-1]
        key_split = keys.split(".")
        if not in_loop:
          nested_d = reduce(lambda x, y: {y: x}, key_split[::-1], value)
          building_d[current_data_key] = merge_dicts(building_d[current_data_key],nested_d)
        else:
          loop_keys.append(key_split)
      
      elif len(line)>0 and line[:4]== "loop":
        loop_start = i
        in_loop = True
      else:
        if in_loop:
          strip = line.strip()
          if loop_data_start == None:
            loop_data_start = i
          elif (len(strip) == 0 or strip == "#" or i==lines_len):
            if i==lines_len:
              i = i+1            
            loop_data = list(zip_longest(*[split_with_quotes(line) for line in lines[loop_data_start:i]]))
        
      
            loop_data = [list(e) for e in loop_data]
            loop_d = {"line_start":loop_start,
                         "line_data_start":loop_data_start,
                         "line_end": i,
                         "loop_keys":loop_keys,
                         "loop_lines":lines[loop_data_start:i],
                         "loop_data":loop_data}
            if len(loop_data) != len(loop_keys):
              failures.append(loop_d)
            else:
              loops.append(loop_d)
            in_loop = False
            loop_start = None
            loop_data_start = None
            loop_keys = []
          else:
            pass
  #print([loop["loop_keys"] for loop in failures])
  # load data into dicts
  loop_dicts = []
  for loop in loops:
  
      keys = loop["loop_keys"]
      loop_d = {}
      for i,keys in enumerate(loop["loop_keys"]):
        set_value_at_nested_key(loop_d,keys,loop["loop_data"][i])
      loop_dicts.append(loop_d)
      
      # add to building dict
      for d in loop_dicts:
          keys = list(d.keys())
          assert len(keys)==1
          key = keys[0]
          #assert key not in building_d[current_data_key] # ??
          building_d[current_data_key][key] = d[key]
  return building_d


# Define a function to pad each cell in the DataFrame
def pad_with_spaces(cell,decimal_places=3,integer_padding=2,decimal_padding=2):
    if isinstance(cell,(int)):
      cell = str(cell)
      return cell.ljust(integer_padding + decimal_padding + 1)  # +1 for the decimal point

    if isinstance(cell, (float)):
        # Split the number into integer and decimal parts
        integer_part, decimal_part = f"{round(cell, decimal_places):.0{decimal_places}f}".split('.')
        # Pad the integer part to the left and the decimal part to the right
        return f"{integer_part.rjust(integer_padding)}.{decimal_part.ljust(decimal_padding)}"
    elif isinstance(cell, str):
        # Pad the string to the right
        return cell.ljust(integer_padding + decimal_padding + 1)  # +1 for the decimal point
    else:
        # Convert the cell to a string and pad it to the right
        return str(cell).ljust(integer_padding + decimal_padding + 1)  # +1 for the decimal point
      
def df_to_cif_string(df,decimal_places=3,integer_padding=4,decimal_padding=4,column_prefix=None,data_name=None):
  lines = df_to_cif_lines(df,
                          decimal_places=decimal_places,
                          integer_padding=integer_padding,
                          decimal_padding=decimal_padding,
                          column_prefix=column_prefix,
                          data_name=data_name)
  return "\n".join(lines)

def df_to_cif_file(df,filename,decimal_places=3,integer_padding=4,decimal_padding=4,column_prefix=None,data_name=None,check=True,suppress=False):
  s = df_to_cif_string(df,
                       decimal_places=decimal_places,
                       integer_padding=integer_padding,
                       decimal_padding=decimal_padding,
                       column_prefix=column_prefix,
                       data_name=data_name)
  with open(filename,"w") as fh:
    fh.write(s)

  # roundtrip
  if check:
    df_in,df_out = df_to_cif_roundtrip_check(df,suppress=suppress)
    if suppress:
      return df_in, df_out
def df_to_cif_lines(df,decimal_places=3,integer_padding=4,decimal_padding=4,column_prefix=None,data_name=None):
  
  # replace python uncertain values with cif ones
  df.replace(to_replace=[None,""," "], value='.', inplace=True)
  
  # convert to string and pad
  df = df.applymap(lambda x: pad_with_spaces(x,
                                             decimal_places=decimal_places,
                                             integer_padding=integer_padding,
                                             decimal_padding=decimal_padding))

  if len(df)>1: #a loop
    lines_header = ["#","loop_"]
    for column in df.columns:
      if column_prefix is not None:
        column = column_prefix+"."+column
      if not column.startswith("_"):
        column = "_"+column
      lines_header.append(column)

    lines_body = [' '.join(map(str, row)) for row in df.values]
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
      line = column_prefix_key.ljust(max_len+2)+value.ljust(max_len)
      lines.append(line)

  if data_name is not None:
    lines = ["data_"+data_name]+lines
  return lines

def dict_to_cif_file(dict,filename,method="lastcif"):
  assert method in ["pdbe","lastcif"]
  if method=="pdbe":
    CifFileWriter("filename").write(dict)
  else:
    raise NotImplementedError("Only pdbe cif is fully functional")
  

def df_to_cif_roundtrip_check(df,filename,suppress=False):
    df = df.copy()
    df_to_cif_file(df,filename)
    d = parse_cif_file(filename)
    assert len(d.keys())==1, "Multiple head keys (data blocks) not supported"
    head_key = list(d.keys())[0]
    data = d[head_key]
    
    assert len(data.keys())==1, "Multiple tables not supported for this function"
    table_key = list(data.keys())[0]
    d = {table_key+"."+key:value for key,value in d[head_key][table_key].items()}
    df_backin = pd.DataFrame(d)
    for column in df_backin.columns:
        df_backin[column] = pd.to_numeric(df_backin[column],errors="raise")
    min_int = sys.maxsize * -1 - 1
    df_backin.fillna(min_int,inplace=True)

    #df = df.reset_index(drop=True)
    #df_backin = df_backin.reset_index(drop=True)
    if not suppress:
        assert all(df==df_backin), "Round trip df -> cif -> df failure. Use suppress keyword to debug"
    return df,df_backin

def write_cif_file(d,filename,method="pdbe"):
  assert method in ["pdbe"]
  if method == "pdbe":
    CifFileWriter(str(filename)).write(d)



  # # this uses functions from lastcif.py that need to be moved
  # # go from cif dict to list of dfs
  # assert len(list(td.keys()))==1, "Multiple data blocks not supported"
  # head_key = list(td.keys())[0]
  
  # # find second to last keys
  # all_paths = find_paths(td.data)
  
  # unique_paths = set()
  # for path in all_paths:
  #     if len(path)>1:
  #       p = path[:-1]
  #     else:
  #       p = path
  #     unique_paths.add(tuple(p))
  
  # dicts = {p[-1]:resolve_value(td.data,p) for p in unique_paths}
  # dfs = {}
  # # check for scaler values and make dfs
  # for k,d in dicts.items():
  #     keys = list(d.keys())
  #     is_scaler = False
  #     for key in keys:
  #         v_check = d[key]
  #         if isinstance(v_check,str) or len(v_check)==1:
  #             is_scaler = True
  #         if is_scaler:
  #             df = pd.DataFrame(d,index=[0])
  #         else:
  #             df = pd.DataFrame(d)
  #         dfs[k] = df
  
  # # # write to disk
  # cif_lines_all = ["data_"+head_key]
  # for key,df in dfs.items():
  
  #     cif_lines = df_to_cif_lines(df,column_prefix=key,
  #                                 decimal_places=3,
  #                                 integer_padding=2,
  #                                 decimal_padding=2)
  #     cif_lines_all+=cif_lines
  # with outfile.open("w") as fh:
  #   fh.write("\n".join(cif_lines))