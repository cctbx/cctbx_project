from __future__ import absolute_import, division, print_function
import sys

hpdl_database = {
    "Only ND1 protonated" : {
        ("CG",  "ND1") : [1.371, 0.009],
        ("ND1", "CE1") : [1.338, 0.008],
        ("CE1", "NE2") : [1.316, 0.009],
        ("NE2", "CD2") : [1.377, 0.007],
        ("CD2", "CG")  : [1.356, 0.009],
        ("CB",  "CG")  : [1.477, 0.019],
        ("ND1", "CG", "CD2") : [105.1, 0.5],
        ("CG", "ND1", "CE1") : [108.1, 0.5],
        ("ND1", "CE1", "NE2"): [111.3, 0.6],
        ("CE1", "NE2", "CD2"): [105.3, 0.6],
        ("NE2", "CD2", "CG") : [110.3, 0.6],
        ("CB", "CG", "ND1")  : [123.4, 1.6],
        ("CB", "CG", "CD2")  : [131.5, 1.6],
    },
    "Only NE2 protonated" : {
        ("CG",  "ND1") : [1.382, 0.009],
        ("ND1", "CE1") : [1.320, 0.008],
        ("CE1", "NE2") : [1.335, 0.014],
        ("NE2", "CD2") : [1.367, 0.009],
        ("CD2", "CG")  : [1.359, 0.010],
        ("CB",  "CG")  : [1.481, 0.022],
        ("ND1", "CG", "CD2") : [109.1, 0.5],
        ("CG", "ND1", "CE1") : [105.3, 0.9],
        ("ND1", "CE1", "NE2"): [111.9, 1.0],
        ("CE1", "NE2", "CD2"): [107.2, 0.8],
        ("NE2", "CD2", "CG") : [106.5, 0.7],
        ("CB", "CG", "ND1")  : [121.6, 1.3],
        ("CB", "CG", "CD2")  : [129.3, 1.4],
    },
    "ND1 and NE2 protonated" : {
        ("CG",  "ND1") : [1.379, 0.007],
        ("ND1", "CE1") : [1.325, 0.008],
        ("CE1", "NE2") : [1.316, 0.009],
        ("NE2", "CD2") : [1.373, 0.007],
        ("CD2", "CG")  : [1.353, 0.007],
        ("CB",  "CG")  : [1.490, 0.011],
        ("ND1", "CG", "CD2") : [105.9, 0.5],
        ("CG", "ND1", "CE1") : [109.3, 0.5],
        ("ND1", "CE1", "NE2"): [108.4, 0.5],
        ("CE1", "NE2", "CD2"): [108.9, 0.6],
        ("NE2", "CD2", "CG") : [107.4, 0.5],
        ("CB", "CG", "ND1")  : [122.7, 1.2],
        ("CB", "CG", "CD2")  : [131.3, 1.3],
    },
  }
"""
 This value seems to be in error. The corrected value for the EH99 NE2-CD2 bond length is 1.372  .
}
"""

def geometric_hydrogens():
  angles = [
    ("HD1", "ND1", "CG"),
    ("HD1", "ND1", "CE1"),
    ("HE1", "CE1", "ND1"),
    ("HE1", "CE1", "NE2"),
    ("HE2", "NE2", "CE1"),
    ("HE2", "NE2", "CD2"),
    ("HD2", "CD2", "NE2"),
    ("HD2", "CD2", "CG"),
    ]
  def _geometric_hydrogens(protonation):
    tmp = {}
    for angle in angles:
      for key, item in hpdl_database[protonation].items():
        if len(key)!=3: continue
        has_h=False
        for ta in angles:
          if ta[0] in key:
            has_h=True
            break
        if has_h: continue
        if angle[1]==key[1]:
          esd = item[1]
          na = (360-item[0])/2
          tmp[angle]=[na, esd]
    #for key in sorted(tmp):
    #  print key, tmp[key]
    return tmp

  hpdl_h_database = {}
  for key in hpdl_database:
    hpdl_database[key].update(_geometric_hydrogens(key))

def get_hpdl_database(include_hydrogens=True,
                      reasonable_esds=True,
                      ):
  if include_hydrogens: geometric_hydrogens()
  if reasonable_esds:
    for key, item in hpdl_database.items():
      for ic, values in item.items():
        if len(ic)==2:
          limit=0.01
          factor=2
        elif len(ic)==3:
          limit=1.5
          factor=2
        if values[1]<limit:
          values[1] = limit
        #else:
        #  values[1] *= factor
  return hpdl_database

def run(args):
  assert len(args) == 0
  print(hpdl_database["Only ND1 protonated"])
  for res_type in sorted(hpdl_database):
    print(res_type, len(hpdl_database[res_type]))
  geometric_hydrogens()
  print(hpdl_database["Only ND1 protonated"])
  for res_type in sorted(hpdl_database):
    print(res_type, len(hpdl_database[res_type]))

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
