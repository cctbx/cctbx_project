
import numpy as np
import pandas as pd

from iotbx.pdb.amino_acid_codes import three_letter_given_one_letter
from libtbx import group_args

"""
Parameters for viewers module

Concepts:
    1. core attributes: attributes that are unambiguous and must be fully present to define an atom

"""

core_map_to_mmcif = {
# Maps a 'core' attribute key to a default 'real' attribute key
'asym_id':"auth_asym_id",
'seq_id':"auth_seq_id",
'comp_id':"label_comp_id",
'atom_id':"label_atom_id",
'type_symbol':"type_symbol",
'alt_id':'label_alt_id',
'id':'id',
}
core_map_to_core = {v:k for k,v in core_map_to_mmcif.items()}

# Maps mmcif attribute name to phenix name
attrs_map_to_phenix = {
'auth_asym_id': 'chain',
'label_asym_id': 'chain',
'auth_seq_id': 'resseq',
'label_seq_id': 'resseq',
'label_comp_id': 'resname',
'type_symbol': 'element',
'label_alt_id': 'altloc',
'pdbx_PDB_model_num': 'model',
'pdbx_PDB_ins_code': 'icode',
'Cartn_x': 'x',
'Cartn_y': 'y',
'Cartn_z': 'z',
'B_iso_or_equiv': 'bfactor',
'occupancy': 'occupancy',
'pdbx_formal_charge': 'charge',
'group_PDB': 'group',
'id': 'id',
'label_atom_id': 'name',
'label_entity_id': 'segid'
}
# The reverse, 
attrs_map_to_mmcif = {}
for k,v in attrs_map_to_phenix.items():
    core_key = k.replace("auth_","").replace("label_","")
    if core_key in core_map_to_mmcif:
        mmcif_key = core_map_to_mmcif[core_key]
        attrs_map_to_mmcif[v] = mmcif_key # set the core:real map
        attrs_map_to_mmcif["auth_"+core_key] = mmcif_key # alias
        attrs_map_to_mmcif["label_"+core_key] = mmcif_key # alias


    else:
        attrs_map_to_mmcif[v] = k

#Create list of all keywords for parsing
keywords_all = list((set(list(attrs_map_to_mmcif.keys())) | 
                     set(list(attrs_map_to_mmcif.values())) | 
                     set(list(core_map_to_mmcif.keys())) | 
                     set(list(core_map_to_mmcif.values()))))

# Convert logic operators between pandas query and phenix selection
logic_map_to_pandas ={ 
    'or': '|',
    'and': '&',
    'not': '~',
    "and":"&",
    "or":"|",
    "not":"~",
    " ": "==",
    ">=": ">=",
    ">": ">",
    "<=": "<=",
    "<": "<",
}
# The reverse
logic_map_to_phenix = {v:k for k,v in logic_map_to_pandas.items() if v != "=="}

logic_map_molstar = {
    "==": "eq",
    ">=": "gre",
    ">": "gr",
    "<=": "lte",
    "<": "lt",
    "!=": "neq",
}


dtypes_mmcif = {
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
}

# Define functions to fill missing columns.
# 'None' indicates an absolutely required initial column
fill_functions = {
#'auth_asym_id': lambda sites: sites["label_asym_id"],
#'auth_asym_id': None,#lambda sites: sites["label_asym_id"],
#'label_asym_id': None,

#'auth_seq_id': lambda sites: sites["label_seq_id"],
#'auth_seq_id': None, #lambda sites: sites["label_seq_id"],
#'label_seq_id': None,

#'auth_comp_id': lambda sites: sites["label_comp_id"],
#'auth_comp_id':None, #lambda sites: sites["label_comp_id"],
#'label_comp_id': None,

'type_symbol': None,
'label_alt_id': None,
'pdbx_PDB_model_num': None,
'pdbx_PDB_ins_code': None,
'Cartn_x': None,
'Cartn_y': None,
'Cartn_z': None,
'B_iso_or_equiv': None,
'occupancy': None,
'pdbx_formal_charge': None,
'group_PDB': None,
'id': None,
'label_atom_id': None,
'label_entity_id': None,
}

# Regularize names to a consistent default name (no auth,label prefix)
attr_aliases = core_map_to_mmcif # core:real
#attr_aliases_reverse = {v:k for k,v in core_keys_to_mmcif_keys_default.items()}, # real:core

attrs_core = [ # attributes assumed to exists
'asym_id',
'seq_id',
'comp_id',
'atom_id',
'alt_id',
'id',
'type_symbol',
'occupancy',
'B_iso_or_equiv'
]
attrs_core_compositional = [
'asym_id',
'comp_id',
'atom_id',
'alt_id',
'type_symbol',
'seq_id',
]

attrs_compositional = []
for attr in attrs_core_compositional:
    if attr in core_map_to_mmcif:
        attr = core_map_to_mmcif[attr]
    attrs_compositional.append(attr)

attr_fill_values = {
# required
'asym_id':"?",
'seq_id':"?",
'comp_id':"?",
'atom_id':"?",
'id':"?",
'type_symbol':"?",

# other
'x':np.nan,
'y':np.nan,
'z':np.nan,
'group_id':"ATOM",
'entity_id':"",
'model_id':"1",
'occupancy':np.nan,
'bfactor':np.nan,
'icode':" ",
'charge':"",
'alt_id':"",
}
rounding = {
"Cartn_x":3,
"Cartn_y":3,
"Cartn_Z":3,
"B_iso_or_equiv":2,
"occupancy":2,

}
blanks = set([""," ",".","?",np.nan,pd.NA])


protein_comp_ids = [v.upper() for v in list(three_letter_given_one_letter.values())]


# Solvent/Ligands

solv_group_1 = ["HOH"]

solv_group_2 = [
  "NA",  # Sodium
  "K",   # Potassium
  "MG",  # Magnesium
  "CA",  # Calcium
  "ZN",  # Zinc
  "FE",  # Iron
  "MN",  # Manganese
  "CO",  # Cobalt
  "CU",  # Copper
  "NI",  # Nickel
  "AG",  # Silver
  "CD",  # Cadmium
  "HG"   # Mercury
]

solv_group_3 = [
  "CL",  # Chloride
  "BR",  # Bromide
  "I",   # Iodide
  "F"    # Fluoride
]

solv_group_4 = [
  "SO4",  # Sulfate
  "PO4",  # Phosphate
  "NO3",  # Nitrate
  "NH4"   # Ammonium
]

solvent_comp_ids = solv_group_1 + solv_group_2 + solv_group_3 + solv_group_4


params = group_args(**locals())