from __future__ import absolute_import, division, print_function
from mmtbx.conformation_dependent_library.bond_angle_registry import \
  bond_angle_registry

not_before_pro_groups = {
  "NonPGIV_nonxpro" : ["ALA",
                       "ARG",
                       "ASN",
                       "ASP",
                       "CYS",
                       "GLN",
                       "GLU",
                       "HIS",
                       "LEU",
                       "LYS",
                       "MET",
                       "PHE",
                       "SER",
                       "THR",
                       "TRP",
                       "TYR",
                       ],
  "IleVal_nonxpro" : ["ILE",
                      "VAL",
                      ],
  "Gly_nonxpro" : ["GLY"],
  "Pro_nonxpro" : ["PRO"],
}
before_pro_groups = {
  "NonPGIV_xpro" : not_before_pro_groups["NonPGIV_nonxpro"],
  "IleVal_xpro"  : not_before_pro_groups["IleVal_nonxpro"],
  "Gly_xpro"     : not_before_pro_groups["Gly_nonxpro"],
  "Pro_xpro"     : not_before_pro_groups["Pro_nonxpro"],
}
columns = [
  "",
  "",
  "mCNA", # C(-1) - N(0)  - Ca(0)
  "sCNA",
  "mNAB", # NAB   N(0)  - Ca(0) - Cb(0)
  "sNAB",
  "mNAC", # NAC   N(0)  - Ca(0) - C(0)
  "sNAC",
  "mBAC", # BAC   Cb(0) - Ca(0) - C(0)
  "sBAC",
  "mACO", # ACO   Ca(0) - C(0)  - O(0)
  "sACO",
  "mACN", # ACN   Ca(0) - C(0)  - N(+1)
  "sACN",
  "mOCN", # OCN   O(0)  - C(0)  - N(+1)
  "sOCN",
  "mCN",  # CN    C(-1) - N(0)
  "sCN",
  "mNA",  # NA    N(0)  - Ca(0)
  "sNA",
  "mAB",  # AB    Ca(0) - Cb(0)
  "sAB",
  "mAC",  # AC    Ca(0) - C(0)
  "sAC",
  "mCO",  # CO    C(0)  - O(0)
  "sCO",
  # needed for cis_127
  'mCND', # C(-1) - N(0) - Cd(0)
  'sCND',
  'mAND', # Ca(0) - N(0) - Cd(0)
  'sAND',
  'mNDG', # N(0) - Cd(0) - Cg(0)
  'sNDG',
  'mABG', # Ca(0) - Cb(0) - Cg(0)
  'sABG',
  'mBGD', # Cb(0) - Cg(0) - Cd(0)
  'sBGD',
  'mBG',
  'sBG',
  'mGD',
  'sGD',
  'mND',
  'sND',
  ]
headers = [
  "statistical type",       # 0
  "number",                 # 1
  "C(-1) - N(0)  - Ca(0)",  # 2
  "",
  "N(0)  - Ca(0) - Cb(0)",  # 4
  "",
  "N(0)  - Ca(0) - C(0)",   # 6
  "",
  "Cb(0) - Ca(0) - C(0)",   # 8
  "",
  "Ca(0) - C(0)  - O(0)",
  "",
  "Ca(0) - C(0)  - N(+1)",
  "",
  "O(0)  - C(0)  - N(+1)",
  "",
  "C(-1) - N(0)",
  "",
  "N(0)  - Ca(0)",
  "",
  "Ca(0) - Cb(0)",
  "",
  "Ca(0) - C(0)",
  "",
  "C(0)  - O(0)",
  "",
  ]

def setup_restraints(geometry, # restraints_manager
                     verbose=False,
                     ):
  ba_registry = bond_angle_registry()
  for angle in geometry.angle_proxies:
    ba_registry[angle.i_seqs]=angle
  return ba_registry

