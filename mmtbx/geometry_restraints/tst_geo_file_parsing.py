from __future__ import division, print_function
import json

from geo_file_parsing import parse_geo_file


"""
Test that parsing a geo file yields previous results
"""
sample_geo_text = """
# Geometry restraints

Bond restraints: 2
bond 38
     39
  ideal  model  delta    sigma   weight residual
  1.522  1.553 -0.030 1.18e-02 7.18e+03 6.53e+00
bond 37
     38
  ideal  model  delta    sigma   weight residual
  1.460  1.485 -0.025 1.17e-02 7.31e+03 4.40e+00

Bond angle restraints: 2
Sorted by residual:
angle 23
      24
      25
    ideal   model   delta    sigma   weight residual
   108.90  113.48   -4.58 1.63e+00 3.76e-01 7.90e+00
angle 37
      38
      39
    ideal   model   delta    sigma   weight residual
   108.02  111.93   -3.91 1.78e+00 3.16e-01 4.84e+00


Dihedral angle restraints: 2
  sinusoidal: 1
    harmonic: 1
Sorted by residual:
dihedral 24
         25
         37
         38
    ideal   model   delta  harmonic     sigma   weight residual
   180.00  166.21   13.79     0      5.00e+00 4.00e-02 7.60e+00
dihedral 58
         59
         60
         61
    ideal   model   delta sinusoidal    sigma   weight residual
     0.00  -72.39   72.39     2      3.00e+01 1.11e-03 4.85e+00


C-Beta improper torsion angle restraints: 2
Sorted by residual:
dihedral 37
         39
         38
         41
    ideal   model   delta  harmonic     sigma   weight residual
   122.80  126.95   -4.15     0      2.50e+00 1.60e-01 2.75e+00
dihedral 56
         54
         55
         58
    ideal   model   delta  harmonic     sigma   weight residual
  -122.60 -126.01    3.41     0      2.50e+00 1.60e-01 1.86e+00


Chirality restraints: 2
Sorted by residual:
chirality 55
          54
          56
          58
  both_signs  ideal   model   delta    sigma   weight residual
    False      2.51    2.39    0.12 2.00e-01 2.50e+01 3.48e-01
chirality 72
          71
          73
          75
  both_signs  ideal   model   delta    sigma   weight residual
    False      2.51    2.62   -0.11 2.00e-01 2.50e+01 2.86e-01

Planarity restraints: 2
Sorted by residual:
             delta    sigma   weight rms_deltas residual
plane 89     0.007 2.00e-02 2.50e+03   1.43e-02 6.15e+00
      90     0.029 2.00e-02 2.50e+03
      91    -0.003 2.00e-02 2.50e+03
      92     0.001 2.00e-02 2.50e+03
      93     0.003 2.00e-02 2.50e+03
      94    -0.001 2.00e-02 2.50e+03
      95    -0.013 2.00e-02 2.50e+03
      96    -0.001 2.00e-02 2.50e+03
      102   -0.029 2.00e-02 2.50e+03
      103   -0.016 2.00e-02 2.50e+03
      104    0.017 2.00e-02 2.50e+03
      105    0.005 2.00e-02 2.50e+03
            delta    sigma   weight rms_deltas residual
plane 13   -0.003 2.00e-02 2.50e+03   1.55e-02 3.62e+00
      14   -0.022 2.00e-02 2.50e+03
      15    0.019 2.00e-02 2.50e+03
      16   -0.002 2.00e-02 2.50e+03
      19   -0.012 2.00e-02 2.50e+03
      20    0.020 2.00e-02 2.50e+03


Nonbonded interactions: 3
Sorted by model distance:
nonbonded 106
          110
   model   vdw sym.op.
   1.719 1.850 -x+1,y-1/2,-z+1
nonbonded 29
          34
   model   vdw sym.op.
   1.859 1.850 x,y+1,z

nonbonded 33
          61
   model   vdw
   1.876 1.850

"""

test_data_json = """
{
  "nonbonded": [
    {
      "restraint_type": "nonbonded",
      "i_seq": 106,
      "j_seq": 110,
      "model": 1.719,
      "vdw": 1.85,
      "sym.op.": "-x+1,y-1/2,-z+1"
    },
    {
      "restraint_type": "nonbonded",
      "i_seq": 29,
      "j_seq": 34,
      "model": 1.859,
      "vdw": 1.85,
      "sym.op.": "x,y+1,z"
    },
    {
      "restraint_type": "nonbonded",
      "i_seq": 33,
      "j_seq": 61,
      "model": 1.876,
      "vdw": 1.85,
      "sym.op.": null
    }
  ],
  "angle": [
    {
      "restraint_type": "angle",
      "i_seq": 23,
      "j_seq": 24,
      "k_seq": 25,
      "ideal": 108.9,
      "model": 113.48,
      "delta": -4.58,
      "sigma": 1.63,
      "weight": 0.376,
      "residual": 7.9
    },
    {
      "restraint_type": "angle",
      "i_seq": 37,
      "j_seq": 38,
      "k_seq": 39,
      "ideal": 108.02,
      "model": 111.93,
      "delta": -3.91,
      "sigma": 1.78,
      "weight": 0.316,
      "residual": 4.84
    }
  ],
  "bond": [
    {
      "restraint_type": "bond",
      "i_seq": 38,
      "j_seq": 39,
      "ideal": 1.522,
      "model": 1.553,
      "delta": -0.03,
      "sigma": 0.0118,
      "weight": 7180.0,
      "residual": 6.53
    },
    {
      "restraint_type": "bond",
      "i_seq": 37,
      "j_seq": 38,
      "ideal": 1.46,
      "model": 1.485,
      "delta": -0.025,
      "sigma": 0.0117,
      "weight": 7310.0,
      "residual": 4.4
    }
  ],
  "dihedral": [
    {
      "restraint_type": "dihedral",
      "i_seq": 24,
      "j_seq": 25,
      "k_seq": 37,
      "l_seq": 38,
      "ideal": 180.0,
      "model": 166.21,
      "delta": 13.79,
      "harmonic": 0.0,
      "sigma": 5.0,
      "weight": 0.04,
      "residual": 7.6,
      "sinusoidal": null
    },
    {
      "restraint_type": "dihedral",
      "i_seq": 58,
      "j_seq": 59,
      "k_seq": 60,
      "l_seq": 61,
      "ideal": 0.0,
      "model": -72.39,
      "delta": 72.39,
      "harmonic": null,
      "sigma": 30.0,
      "weight": 0.00111,
      "residual": 4.85,
      "sinusoidal": 2.0
    }
  ],
  "chirality": [
    {
      "restraint_type": "chirality",
      "i_seq": 55,
      "j_seq": 54,
      "k_seq": 56,
      "l_seq": 58,
      "both_signs": "False",
      "ideal": 2.51,
      "model": 2.39,
      "delta": 0.12,
      "sigma": 0.2,
      "weight": 25.0,
      "residual": 0.348
    },
    {
      "restraint_type": "chirality",
      "i_seq": 72,
      "j_seq": 71,
      "k_seq": 73,
      "l_seq": 75,
      "both_signs": "False",
      "ideal": 2.51,
      "model": 2.62,
      "delta": -0.11,
      "sigma": 0.2,
      "weight": 25.0,
      "residual": 0.286
    }
  ],
  "c-beta": [
    {
      "restraint_type": "c-beta",
      "i_seq": 37,
      "j_seq": 39,
      "k_seq": 38,
      "l_seq": 41,
      "ideal": 122.8,
      "model": 126.95,
      "delta": -4.15,
      "harmonic": 0.0,
      "sigma": 2.5,
      "weight": 0.16,
      "residual": 2.75
    },
    {
      "restraint_type": "c-beta",
      "i_seq": 56,
      "j_seq": 54,
      "k_seq": 55,
      "l_seq": 58,
      "ideal": -122.6,
      "model": -126.01,
      "delta": 3.41,
      "harmonic": 0.0,
      "sigma": 2.5,
      "weight": 0.16,
      "residual": 1.86
    }
  ]
}

"""

test_data = json.loads(test_data_json)

# load file
with open("tst_geo_parsing.geo","w") as fh:
  fh.write(sample_geo_text)

# read back in
geo_dict = parse_geo_file("tst_geo_parsing.geo",return_format='dict')

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
  #print(d1)
  #print(d2)
  for key1,value1 in d1.items():
    assert key1 in d2
    value2 = d2[key1]
    #print(value1,value2)
    t1,t2 = get_type(value1), get_type(value2)
    assert t1==t2, '%s!=%s' % (t1, t2)
    assert str(value1)==str(value2), "Invalid comparisons: "+str(value1)+" and "+str(value2)
  return True



# test
def test():
  for key,d_list in geo_dict.items():
    for j,d in enumerate(d_list):
      d_ref = test_data[key][j]
      same = compare_custom(d,d_ref)
      print(same)
if __name__ == '__main__':
  test()
  print("Finished")
