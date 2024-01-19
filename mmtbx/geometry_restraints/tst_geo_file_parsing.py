from __future__ import division, print_function
import os
import json


import libtbx.load_env
from iotbx.data_manager import DataManager

from geo_file_parsing import parse_geo_file


"""
Test that parsing a geo file yields previous results
"""
test_data_json = """
{
  "_phenix_restraint_nonbonded": [
    {
      "restraint_type": "nonbonded",
      "i_seq": 0,
      "j_seq": 3,
      "model": 2.665,
      "vdw": 2.496
    },
    {
      "restraint_type": "nonbonded",
      "i_seq": 19,
      "j_seq": 36,
      "model": 2.726,
      "vdw": 2.52
    },
    {
      "restraint_type": "nonbonded",
      "i_seq": 41,
      "j_seq": 42,
      "model": 2.803,
      "vdw": 2.752
    }
  ],
  "_phenix_restraint_angle": [
    {
      "restraint_type": "angle",
      "i_seq": 12,
      "j_seq": 13,
      "k_seq": 14,
      "ideal": 108.9,
      "model": 113.46,
      "delta": -4.56,
      "sigma": 1.63,
      "weight": 0.376,
      "residual": 7.83
    },
    {
      "restraint_type": "angle",
      "i_seq": 21,
      "j_seq": 22,
      "k_seq": 23,
      "ideal": 120.33,
      "model": 122.26,
      "delta": -1.93,
      "sigma": 1.08,
      "weight": 0.857,
      "residual": 3.18
    },
    {
      "restraint_type": "angle",
      "i_seq": 29,
      "j_seq": 30,
      "k_seq": 31,
      "ideal": 108.9,
      "model": 111.01,
      "delta": -2.11,
      "sigma": 1.63,
      "weight": 0.376,
      "residual": 1.68
    }
  ],
  "_phenix_restraint_bond": [
    {
      "restraint_type": "bond",
      "i_seq": 0,
      "j_seq": 1,
      "ideal": 1.451,
      "model": 1.507,
      "delta": -0.056,
      "sigma": 0.016,
      "weight": 3910.0,
      "residual": 12.1
    },
    {
      "restraint_type": "bond",
      "i_seq": 5,
      "j_seq": 6,
      "ideal": 1.524,
      "model": 1.498,
      "delta": 0.026,
      "sigma": 0.0126,
      "weight": 6300.0,
      "residual": 4.15
    },
    {
      "restraint_type": "bond",
      "i_seq": 35,
      "j_seq": 37,
      "ideal": 1.328,
      "model": 1.348,
      "delta": -0.02,
      "sigma": 0.021,
      "weight": 2270.0,
      "residual": 0.875
    }
  ],
  "_phenix_restraint_dihedral": [
    {
      "restraint_type": "dihedral",
      "i_seq": 13,
      "j_seq": 14,
      "k_seq": 20,
      "l_seq": 21,
      "ideal": 180.0,
      "model": 166.19,
      "delta": 13.81,
      "harmonic": 0.0,
      "sigma": 5.0,
      "weight": 0.04,
      "residual": 7.63,
      "sinusoidal": null
    },
    {
      "restraint_type": "dihedral",
      "i_seq": 1,
      "j_seq": 2,
      "k_seq": 4,
      "l_seq": 5,
      "ideal": 180.0,
      "model": 171.51,
      "delta": 8.49,
      "harmonic": 0.0,
      "sigma": 5.0,
      "weight": 0.04,
      "residual": 2.88,
      "sinusoidal": null
    },
    {
      "restraint_type": "dihedral",
      "i_seq": 20,
      "j_seq": 21,
      "k_seq": 24,
      "l_seq": 25,
      "ideal": -180.0,
      "model": -164.64,
      "delta": -15.36,
      "harmonic": null,
      "sigma": 15.0,
      "weight": 0.00444,
      "residual": 1.45,
      "sinusoidal": 3.0
    }
  ]
}
"""

test_data = json.loads(test_data_json)

# load file
pdb_file = libtbx.env.find_in_repositories(
  relative_path="phenix_regression/pdb/1yjp_h.pdb",
  test=os.path.isfile)

dm = DataManager()
dm.process_model_file(pdb_file)
model = dm.get_model()

# write geo
model.add_crystal_symmetry_if_necessary()
model.process(make_restraints=True)
grm = model.get_restraints_manager().geometry
grm.write_geo_file(model.get_sites_cart(),file_name="test.geo")

# read back in
geo_dict = parse_geo_file("test.geo",return_format='dict')

# define comparison functions

def get_type(value):
  try:
    f = float(value)
    if f.is_integer():
      return int
    else:
      return float
  except Exception:
    return str

def compare_custom(d1,d2):
  """
  Compare two dictionaries by types, and 'ideal' floats
  with a large tolerance to enable future changes to restraints.
  """
  print(d1)
  print(d2)
  for key1,value1 in d1.items():
    assert key1 in d2
    value2 = d2[key1]
    print(value1,value2)
    t1,t2 = get_type(value1), get_type(value2)
    assert t1==t2, '%s!=%s' % (t1, t2)
    if isinstance(t1,float):
      if key1 == 'ideal':
        diff = abs(float(value1)-float(value2))
        assert diff < 5, "difference in expected and generated restraint out of range" # 5 ok?
  return True


# test
def test():
  for key,d_list in geo_dict.items():
    spot_indices = [0,3,11]
    for idx,i in enumerate(spot_indices):
      d = d_list[i]
      d_ref = test_data[key][idx]
      same = compare_custom(d,d_ref)
      print(same)

if __name__ == '__main__':
  test()
