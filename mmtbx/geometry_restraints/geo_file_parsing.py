"""
Parse a '.geo' text file to structured restraint data.
"""
from __future__ import division
import re
import pandas as pd
import numpy as np



restraint_settings = {
  "nonbonded":{                              # name of restraint
    "header_label":"Nonbonded interactions", # search string
    "n_atoms":2,                             # n_atoms to look for
    "cif_key":"_phenix_restraint_nonbonded"  # key for restraint in result
  },
  "angle":{
    "header_label":"Bond angle restraints",
    "n_atoms":3,
      "cif_key":"_phenix_restraint_angle"
  },
  "bond":{
    "header_label":"Bond restraints",
    "n_atoms":2,
      "cif_key":"_phenix_restraint_bond"
  },
  "dihedral":{
    "header_label":"Dihedral angle restraints",
    "n_atoms":4,
      "cif_key":"_phenix_restraint_dihedral"
  },

}


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
  groups = geo_lines_to_groups(lines)
  restraint_dicts_all = {}
  for restraint_key,settings in restraint_settings.items():
    header_label = settings["header_label"]
    n_atoms = settings["n_atoms"]
    restraint_dicts = extract_restraint_dicts_from_groups(groups,
      restraint_key=restraint_key,
        header_label=header_label,
          n_atoms=n_atoms)
    reframed_dict = pd.DataFrame.from_records(restraint_dicts).to_dict("list")
    restraint_dicts_all[restraint_key] = reframed_dict
  for key,value in list(restraint_dicts_all.items()):
    cif_key = restraint_settings[key]["cif_key"]
    restraint_dicts_all[cif_key] = restraint_dicts_all.pop(key)

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





def geo_lines_to_groups(lines):
  """
  Parse raw text lines into 'groups' of lines for each restraint type
  """

  # group by newline into restraints groups (bonds, angles, etc)
  groups = []
  group_open = False
  group_lines = []
  for line in lines:
    if line == "\n":
      if group_open == False: # first group
        group_open = True
      else: # close open group
        groups.append(group_lines)
        group_lines = []
    elif group_open:
      group_lines.append(line)
  return groups

def extract_restraint_dicts_from_groups(groups,
                                        restraint_key=None,
                                        header_label=None,
                                        n_atoms=None,
                                        ignore_sym=True):
  """
  Given a 'group' of lines for a given restraint, extract the data to a dictionary

  ignore_sym (bool): Whether to include restraints that cross symetry boundary

  TODO: This only works up to dihedral. Need to implement restraints with uncertain n_atoms
  """
  assert [restraint_key,header_label,n_atoms].count(None)==0, "provide each parameter"
  # split groups into nonbonded restraints
  group = [group for group in groups if header_label in group[0]][0]

  restraints = []
  restraint_open = False
  restraint_lines = []
  for line in group:
    if restraint_key in line and header_label not in line:
      if restraint_open == False:
        restraint_open = True
      else:
        restraints.append(restraint_lines)
        restraint_lines = []
    if restraint_open:
      restraint_lines.append(line)

  s = set([len(r) for r in restraints])
  assert len(s)==1, "Failed .geo parsing"
  assert next(iter(s))==n_atoms+2, "Failed .geo parsing"
  key_index = n_atoms
  value_index = n_atoms+1

  restraint_dicts = []
  atom_indices = ["i","j","k","l","m","n"]
  for restraint in restraints:
    if "sym.op." in restraint[key_index] and ignore_sym:
      continue
    atoms_i_seqs = [int(re.search(r'\d+', restraint[i]).group()) for i in range(0,n_atoms)]

    keys = [e for e in restraint[key_index].split()]
    values = [float(e) for e in restraint[value_index].split()]
    d = {"restraint_type":restraint_key}
    for i, i_seq in enumerate(atoms_i_seqs):
      k = "{}_seq".format(atom_indices[i])
      v = i_seq
      d[k] = v
    d.update({key:value for key,value in zip(keys,values)})

    restraint_dicts.append(d)
  return restraint_dicts
