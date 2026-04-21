from __future__ import absolute_import, division, print_function
import os, sys

cif_keyword_dictionary = {
  "_chem_comp" : {"id" : str,
                  "name" : str,
                  "type" : str,
                  "pdbx_type" : str,
                  "formula" : str,
                  "mon_nstd_parent_comp_id" : str,
                  "pdbx_synonyms" : str,
                  "pdbx_formal_charge" : int,
                  "pdbx_initial_date" : str,
                  "pdbx_modified_date" : str,
                  "pdbx_ambiguous_flag" : str,
                  "pdbx_release_status" : str,
                  "pdbx_replaced_by" : str,
                  "pdbx_replaces" : str,
                  "formula_weight" : float,
                  "one_letter_code" : str,
                  "three_letter_code" : str,
                  "pdbx_model_coordinates_details" : str,
                  "pdbx_model_coordinates_missing_flag" : str,
                  "pdbx_ideal_coordinates_details" : str,
                  "pdbx_ideal_coordinates_missing_flag" : str,
                  "pdbx_model_coordinates_db_code" : str,
                  "pdbx_processing_site" : str,
                  # added 11/2008
                  "pdbx_subcomponent_list" : str,
                  },
  "_chem_comp_atom" : {"comp_id": str,
                       "atom_id": str,
                       "alt_atom_id": str,
                       "type_symbol": str,
                       "charge" : int,
                       "pdbx_align": int,
                       "pdbx_aromatic_flag": str,
                       "pdbx_leaving_atom_flag": str,
                       "pdbx_stereo_config": str,
                       "model_Cartn_x": float,
                       "model_Cartn_y": float,
                       "model_Cartn_z": float,
                       "pdbx_model_Cartn_x_ideal": float,
                       "pdbx_model_Cartn_y_ideal": float,
                       "pdbx_model_Cartn_z_ideal": float,
                       "pdbx_ordinal" : int,
                       # added 11/2008
                       "pdbx_component_atom_id" : str,
                       "pdbx_component_comp_id" : str,
                       },
  "_chem_comp_bond" : {"comp_id": str,
                       "atom_id_1": str,
                       "atom_id_2": str,
                       "value_order": str,
                       "pdbx_aromatic_flag": str,
                       "pdbx_stereo_config": str,
                       "pdbx_ordinal": int,
                       },
  "_pdbx_chem_comp_descriptor" : {"comp_id" : str,
                                  "type" : str,
                                  "program" : str,
                                  "program_version" : str,
                                  "descriptor" : str,
                                  },
  "_pdbx_chem_comp_identifier" : {"comp_id" : str,
                                  "type" : str,
                                  "program" : str,
                                  "program_version" : str,
                                  "identifier" : str,
                                  },
  }

class empty(object):
  def __repr__(self):
    outl = "\nObject"
    for attr in self.__dict__:
      outl += "\n  %s : %s %s" % (attr, getattr(self, attr), type(getattr(self, attr)))
    return outl

  def __len__(self):
    return len(self.__dict__)

def get_func(cif_key, sk):
  func = None
  if cif_key in cif_keyword_dictionary:
    func = cif_keyword_dictionary[cif_key].get(sk, None)
  return func

def get_typed_field(cif_key, sk, field, func=None):
  if func is None:
    func = get_func(cif_key, sk)
  if field not in ["?", "."] and func not in [None, str]:
    field = func(field)
  return field

def run(filename):
  from iotbx import cif
  if not os.path.exists(filename): return None
  complete_cif_data = {}
  cm = cif.reader(filename, strict=False).model()
  for code, rc in cm.items():
    for key, item in rc.items():
      cif_key = key.split('.')[0]
      sk = key.split(".")[1].strip()
      complete_cif_data.setdefault(cif_key, [])
      if type(item)==type(''):
        item=[item]
      else:
        continue
      for i, row in enumerate(item):
        if len(complete_cif_data[cif_key])<i+1:
          complete_cif_data[cif_key].append(empty())
        func = get_func(cif_key, sk)
        row = get_typed_field(cif_key, sk, row, func=func)
        setattr(complete_cif_data[cif_key][i], sk, row)
    for i, loop in enumerate(rc.iterloops()):
      if not loop: continue
      for j, (key, item) in enumerate(loop.items()):
        if not j:
          objs=[]
          for k in range(len(item)):
            objs.append(empty())
          cif_key = key.split('.')[0]
        sk = key.split(".")[1].strip()
        func = get_func(cif_key, sk)
        for k in range(len(item)):
          setattr(objs[k], sk, get_typed_field(cif_key, sk, item[k], func))
        complete_cif_data[cif_key] = objs
  return complete_cif_data
