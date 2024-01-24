"""
Parse a '.geo' text file to structured restraint data.
"""
import re
import pandas as pd
import numpy as np
def tryfloat(value):
  try:
    return float(value)
  except:
    return value
class BondRestraint:
  restraint_label = 'bond'
  @classmethod
  def from_geo_lines(cls,lines,settings):
    # start make restraint
    n_atoms = settings["n_atoms"]
    restraint_key = cls.restraint_label
    key_index = n_atoms
    value_index = n_atoms+1
    atom_indices = ["i","j","k","l","m","n"]
    atoms_i_seqs = [int(re.search(r'\d+', lines[i]).group()) for i in range(0,n_atoms)]
    keys = [e for e in lines[key_index].split()]
    values = [tryfloat(e) for e in lines[value_index].split()]
    d = {"restraint_type":restraint_key}
    for i, i_seq in enumerate(atoms_i_seqs):
     k = "{}_seq".format(atom_indices[i])
     v = i_seq
     d[k] = v
    d.update({key:value for key,value in zip(keys,values)})
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
 #  "plane":{
 #   "header_label":"Planarity restraints",
 #   "n_atoms":None,
 #   "cif_key":"_phenix_restraint_planarity",
 #   "restraint_class":PlanarityRestraint,
 # },
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
        restraints.append(restraint_lines)
        restraint_lines = [line]
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
  dfs = {key:pd.DataFrame(value).replace({np.nan: None}) for key,value in restraint_dicts_all.items()}
  if return_format == "df":
    return dfs
  elif return_format == "dict":
    # convert to 'records' format dict
    return {key:df.to_dict(orient='records') for key,df in dfs.items()}
  elif return_format == "json":
    # same format as dict, but as json string
    return {key:df.to_json(orient='records') for key,df in dfs.items()}
