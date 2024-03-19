from __future__ import division, print_function
import re
import json
from collections import defaultdict
import pandas as pd
import numpy as np

def tryfloat(value):
  # Try to make float
  try:
    return float(value)
  except Exception:
    return value

# Pattern to plit on either space or id_str
pattern = re.compile(r'pdb="[^"]*"|\S+')


# Define classes for restraints with fixed number of atoms
class BondRestraint:
  restraint_label = 'bond'
  @classmethod
  def from_geo_lines(cls, lines, settings):
    # start make restraint
    n_atoms = settings["n_atoms"]
    restraint_key = cls.restraint_label
    key_index = n_atoms
    value_index = n_atoms + 1

    # Extract atom sequences or pdb strings
    atoms_i_seqs = []
    for i in range(n_atoms):
        match = pattern.findall(lines[i].replace(cls.restraint_label,""))
        parts = [part.strip() for part in match if part.strip() != ""]
        if len(parts)>0:
            # Assuming the first match is the desired one; adjust as necessary
            atoms_i_seqs.append(parts[0])
        else:
          atoms_i_seqs.append(lines[i].replace(cls.restraint_label,"").strip())

    keys = pattern.findall(lines[key_index])
    values = [tryfloat(value) for value in pattern.findall(lines[value_index])]
    d = {"restraint_type": restraint_key}
    for i, i_seq in enumerate(atoms_i_seqs):
        if "pdb" in i_seq:
          k = f"id_str_{i}"
        else:
          k = f"i_seq_{i}"
        d[k] = i_seq

    d.update({key: value for key, value in zip(keys, values)})
    return cls(d)

  def __init__(self, data_dict):
    self.data_dict = data_dict
class AngleRestraint(BondRestraint):
  restraint_label = 'angle'
class DihedralRestraint(BondRestraint):
  restraint_label = 'dihedral'
class NonBondedRestraint(BondRestraint):
  restraint_label = 'nonbonded'
class ChiralityRestraint(BondRestraint):
  restraint_label = 'chirality'
class CBetaRestraint(BondRestraint):
  restraint_label = 'c-beta'

# Planes have variable n_atoms
class PlanarityRestraint:
  restraint_label = 'plane'

  @classmethod
  def from_geo_lines(cls, lines, settings):
    planes_data = []
    headers = None
    current_plane_data = {}

    for i,line in enumerate(lines):
      if 'delta' in line:  # Header line detected
        headers = line.split()
        continue  # Skip to the next iteration after setting headers


      # Handling a new "plane" line or a continuation of plane data
      if line.startswith(cls.restraint_label):
        # New "plane" entry; save the previous one if it exists
        if current_plane_data:
          planes_data.append(current_plane_data)
        # Reset current_plane_data for a new "plane" entry
        current_plane_data = {header: [] for header in headers}
        current_plane_data['plane_indices'] = []  # Ensure 'plane_indices' is initialized

      # Split the line using the precompiled pattern, then strip and filter
      matches = [match.group() for match in pattern.finditer(line)]
      parts = [part.strip() for part in matches if part.strip()!=""]
      #print(line)
      #print(parts)

      start_idx = 1 if line.startswith(cls.restraint_label) else 0
      #print("Start idx:",start_idx)
      # Adjust parsing according to the specific needs, for simplicity assuming plane_index is first
      plane_index = parts[start_idx]
      # Append the plane_index after conversion to integer if necessary
      current_plane_data['plane_indices'].append(int(plane_index) if plane_index.isdigit() else plane_index)

      # Assign values to their corresponding headers
      value_idx = start_idx + 1
      #print("Value idx:",value_idx)
      for header, value in zip(headers, parts[value_idx:]):
        current_plane_data[header].append(tryfloat(value))



    # Add the last plane's data if it exists
    #print(current_plane_data)
    if current_plane_data:
      planes_data.append(current_plane_data)


    # another pass to ensure equal length data
    assert len(planes_data)==1
    plane_data = planes_data[0]
    lens = [len(value) for value in plane_data.values()]
    max_len = max(lens)

    for key,value in list(plane_data.items()):
      if len(value)==1:
        plane_data[key] = value*max_len
      if key == "plane_indices":
        if any([isinstance(v,str) and ("pdb" in v) for v in value]):
          plane_data["id_str"] = value
        else:
          plane_data["i_seq"] = value
        del plane_data["plane_indices"]

    return cls({"plane": plane_data})


  def __init__(self, data_dict):
    self.data_dict = data_dict

  @staticmethod
  def form_dataframe(planes):
    planes = [d["plane"] for d in planes]

    # Find max length
    max_len = 0
    for plane in planes:
      lens = [len(value) for value in plane.values()]
      assert len(set(lens))==1, "Error: mismatched data lens"
      if max(lens)>max_len:
        max_len = max(lens)

    # expand each dict
    expanded_planes= []
    for plane in planes:
      plane_len = len(list(plane.values())[0])
      expanded_plane = {}
      for key,value in plane.items():
        cols = []
        values = []
        for i in range(1,max_len+1):
          col_name = f"{key}_{i}"
          if i<=plane_len:
            col_value = value[i-1]
          else:
            col_value = pd.NA
          cols.append(col_name)
          values.append(col_value)
        expanded_plane.update(dict(zip(cols,values)))
      expanded_planes.append(expanded_plane)
    # make df
    df = pd.DataFrame.from_records(expanded_planes)
    return df

restraint_settings = {
 "nonbonded":{                              # name of restraint
  "header_label":"Nonbonded interactions", # search string
  "n_atoms":2,                             # n_atoms to look for
  "cif_key":"_phenix_restraint_nonbonded", # key for restraint in result
  "restraint_class":NonBondedRestraint,
 },
 "angle":{
  "header_label":"Bond angle restraints",
  "n_atoms":3,
   "cif_key":"_phenix_restraint_angle",
   "restraint_class":AngleRestraint,
 },
 "bond":{
  "header_label":"Bond restraints",
  "n_atoms":2,
  "cif_key":"_phenix_restraint_bond",
  "restraint_class":BondRestraint,
 },
 "dihedral":{
  "header_label":"Dihedral angle restraints",
  "n_atoms":4,
  "cif_key":"_phenix_restraint_dihedral",
  "restraint_class":DihedralRestraint,
 },
  "chirality":{
  "header_label":"Chirality restraints",
  "n_atoms":4,
  "cif_key":"_phenix_restraint_chirality",
  "restraint_class":ChiralityRestraint,
 },
  "c-beta":{
  "header_label":"C-Beta improper torsion angle restraints",
  "n_atoms":4,
  "cif_key":"_phenix_restraint_c-beta",
  "restraint_class":CBetaRestraint,
 },
  "plane":{
   "header_label":"Planarity restraints",
   "n_atoms":None,
   "cif_key":"_phenix_restraint_planarity",
   "restraint_class":PlanarityRestraint,
 },
}
def split_into_sections(lines):
  """
  Split a .geo file lines into sub-sections of lines for each restraint type
  """
  sections = {}
  current_section = None
  for idx,line in enumerate(lines):
    line = line.strip()
    if not line:
      continue
    # Check if the line is a section header
    if ':' in line and any(keyword in line for keyword in ['restraints', 'interactions']):
      current_section = line
      sections[current_section] = []
    elif current_section is not None:
      # if 'plane' in line:
      #     sections[current_section].append(lines[idx-1])
      sections[current_section].append(line)
  # Rename labels
  output = {}
  for section_label,section_value in sections.items():
    for setting_label, setting_value in restraint_settings.items():
      search_term = setting_value["header_label"]
      if search_term in section_label:
        output[setting_label] = section_value
  return output

def extract_restraint_objs_from_group(group,
                   settings,
                   restraint_key=None):
  """
  Given a 'group' of lines for a given restraint, extract the data to a restraint instance
  """
  n_atoms = settings["n_atoms"]
  restraints = []
  restraint_lines = []
  if restraint_key == "c-beta":
    new_entry_indicator = 'dihedral'
  else:
    new_entry_indicator = restraint_key

  restraint_label = restraint_key
  for idx,line in enumerate(group):
    # Check if the line is the start of a new entry within the restraint
    if new_entry_indicator in line:
      if restraint_lines:
        # If there are already lines collected for a previous entry, save them
        if restraint_key == "plane":
          restraint_lines = restraint_lines[:-1]
        restraints.append(restraint_lines)
        restraint_lines = []
        if restraint_key == "plane":
          if idx>1:
            restraint_lines.append(group[idx-1])
        restraint_lines.append(line)
      else:
        # This is the first entry in the restraint
        if restraint_key == "plane":
          if idx>1:
            restraint_lines.append(group[idx-1])
        restraint_lines.append(line)
    elif restraint_lines:
      # Continue adding lines to the current entry
      restraint_lines.append(line)
  # Add the last set of lines if not empty
  if restraint_lines:
    if restraint_key == "plane":
      restraint_lines = restraint_lines[:-1]
    restraints.append(restraint_lines)
  s = set([len(r) for r in restraints])
  if restraint_key != "plane":
    assert len(s)==1, "Failed .geo parsing\n"+"\n".join(restraint_lines)
    assert next(iter(s))==n_atoms+2, "Failed .geo parsing"+"\n".join(restraints)
  restraint_objs = []
  for restraint_lines in restraints:
    settings = restraint_settings[restraint_label]
    if "restraint_class" in settings:
      restraint_class = restraint_settings[restraint_label]["restraint_class"]
      restraint_instance = restraint_class.from_geo_lines(restraint_lines,settings)
      restraint_objs.append(restraint_instance)
  return restraint_objs

def parse_geo_file(filename,return_format='dict'):
  """
  Args:
    filename (str): the .geo filename
    return_format (str): One of 'dict', 'df' (pd.DataFrame), or 'json'
  Returns:
    restraint_dicts (dict): nested dictionary with final data in column format
  """
  with open(filename,"r") as fh:
    lines = fh.readlines()
  restraint_dicts_all = {}
  restraint_sections = split_into_sections(lines)
  for restraint_label,setting_dict in restraint_settings.items():
    if restraint_label in restraint_sections:
      group = restraint_sections[restraint_label]
      restraint_objs = extract_restraint_objs_from_group(group,
                            setting_dict,
                            restraint_key=restraint_label)
      restraint_dicts_all[restraint_label] = [obj.data_dict for obj in restraint_objs]
  # convert dicts to dataframes and fill nan with None
  dfs = {key:pd.DataFrame(value).replace({np.nan: None}) for key,value in restraint_dicts_all.items() if key != "plane"}

  # add planes
  if "plane" in restraint_dicts_all:
    df = PlanarityRestraint.form_dataframe(restraint_dicts_all["plane"])
    dfs["plane"] = df
  if return_format == "df":
    return dfs
  elif return_format == "dict":
    # convert to 'records' format dict
    return {key:df.to_dict(orient='records') for key,df in dfs.items()}
  elif return_format == "json":
    # same format as dict, but as json string
    d = {key:df.to_json(orient='records') for key,df in dfs.items()}
    d = {key:json.loads(v) for key,v in d.items()}
    return json.dumps(d,indent=2)



def add_i_seq_columns_from_id_str(restraint_dfs,model):
  """
  Given a dict of pandas DataFrames, each containing one type of restraint data,
  add i_seq columns from a model by translating the existing id_str columns
  """
  restraint_dfs = d
  mapping_dict = {atom.id_str():atom.i_seq for atom in model.get_atoms()}
  for restraint_name,df in restraint_dfs.items():
    id_str_cols = [col for col in df.columns if "id_str" in col]
    i_seq_cols = [col.replace("id_str","i_seq") for col in id_str_cols]
    for i_seq_col,id_str_col in zip(i_seq_cols,id_str_cols):
      df[i_seq_col] = df[id_str_col].map(mapping_dict)
  return restraint_dfs