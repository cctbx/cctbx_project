from __future__ import division, print_function
import json
import pandas as pd
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
      "i_seq_0": "106",
      "i_seq_1": "110",
      "model": 1.719,
      "vdw": 1.85,
      "sym.op.": "-x+1,y-1/2,-z+1"
    },
    {
      "restraint_type": "nonbonded",
      "i_seq_0": "29",
      "i_seq_1": "34",
      "model": 1.859,
      "vdw": 1.85,
      "sym.op.": "x,y+1,z"
    },
    {
      "restraint_type": "nonbonded",
      "i_seq_0": "33",
      "i_seq_1": "61",
      "model": 1.876,
      "vdw": 1.85,
      "sym.op.": null
    }
  ],
  "angle": [
    {
      "restraint_type": "angle",
      "i_seq_0": "23",
      "i_seq_1": "24",
      "i_seq_2": "25",
      "ideal": 108.9,
      "model": 113.48,
      "delta": -4.58,
      "sigma": 1.63,
      "weight": 0.376,
      "residual": 7.9
    },
    {
      "restraint_type": "angle",
      "i_seq_0": "37",
      "i_seq_1": "38",
      "i_seq_2": "39",
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
      "i_seq_0": "38",
      "i_seq_1": "39",
      "ideal": 1.522,
      "model": 1.553,
      "delta": -0.03,
      "sigma": 0.0118,
      "weight": 7180.0,
      "residual": 6.53
    },
    {
      "restraint_type": "bond",
      "i_seq_0": "37",
      "i_seq_1": "38",
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
      "i_seq_0": "24",
      "i_seq_1": "25",
      "i_seq_2": "37",
      "i_seq_3": "38",
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
      "i_seq_0": "58",
      "i_seq_1": "59",
      "i_seq_2": "60",
      "i_seq_3": "61",
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
      "i_seq_0": "55",
      "i_seq_1": "54",
      "i_seq_2": "56",
      "i_seq_3": "58",
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
      "i_seq_0": "72",
      "i_seq_1": "71",
      "i_seq_2": "73",
      "i_seq_3": "75",
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
      "i_seq_0": "dihedral",
      "i_seq_1": "39",
      "i_seq_2": "38",
      "i_seq_3": "41",
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
      "i_seq_0": "dihedral",
      "i_seq_1": "54",
      "i_seq_2": "55",
      "i_seq_3": "58",
      "ideal": -122.6,
      "model": -126.01,
      "delta": 3.41,
      "harmonic": 0.0,
      "sigma": 2.5,
      "weight": 0.16,
      "residual": 1.86
    }
  ],
  "plane": [
    {
      "delta_1": 0.007,
      "delta_2": 0.029,
      "delta_3": -0.003,
      "delta_4": 0.001,
      "delta_5": 0.003,
      "delta_6": -0.001,
      "delta_7": -0.013,
      "delta_8": -0.001,
      "delta_9": -0.029,
      "delta_10": -0.016,
      "delta_11": 0.017,
      "delta_12": 0.005,
      "sigma_1": 0.02,
      "sigma_2": 0.02,
      "sigma_3": 0.02,
      "sigma_4": 0.02,
      "sigma_5": 0.02,
      "sigma_6": 0.02,
      "sigma_7": 0.02,
      "sigma_8": 0.02,
      "sigma_9": 0.02,
      "sigma_10": 0.02,
      "sigma_11": 0.02,
      "sigma_12": 0.02,
      "weight_1": 2500.0,
      "weight_2": 2500.0,
      "weight_3": 2500.0,
      "weight_4": 2500.0,
      "weight_5": 2500.0,
      "weight_6": 2500.0,
      "weight_7": 2500.0,
      "weight_8": 2500.0,
      "weight_9": 2500.0,
      "weight_10": 2500.0,
      "weight_11": 2500.0,
      "weight_12": 2500.0,
      "rms_deltas_1": 0.0143,
      "rms_deltas_2": 0.0143,
      "rms_deltas_3": 0.0143,
      "rms_deltas_4": 0.0143,
      "rms_deltas_5": 0.0143,
      "rms_deltas_6": 0.0143,
      "rms_deltas_7": 0.0143,
      "rms_deltas_8": 0.0143,
      "rms_deltas_9": 0.0143,
      "rms_deltas_10": 0.0143,
      "rms_deltas_11": 0.0143,
      "rms_deltas_12": 0.0143,
      "residual_1": 6.15,
      "residual_2": 6.15,
      "residual_3": 6.15,
      "residual_4": 6.15,
      "residual_5": 6.15,
      "residual_6": 6.15,
      "residual_7": 6.15,
      "residual_8": 6.15,
      "residual_9": 6.15,
      "residual_10": 6.15,
      "residual_11": 6.15,
      "residual_12": 6.15,
      "i_seq_1": 89,
      "i_seq_2": 90,
      "i_seq_3": 91,
      "i_seq_4": 92,
      "i_seq_5": 93,
      "i_seq_6": 94,
      "i_seq_7": 95,
      "i_seq_8": 96,
      "i_seq_9": 102,
      "i_seq_10": 103,
      "i_seq_11": 104,
      "i_seq_12": 105
    },
    {
      "delta_1": -0.003,
      "delta_2": -0.022,
      "delta_3": 0.019,
      "delta_4": -0.002,
      "delta_5": -0.012,
      "delta_6": null,
      "delta_7": null,
      "delta_8": null,
      "delta_9": null,
      "delta_10": null,
      "delta_11": null,
      "delta_12": null,
      "sigma_1": 0.02,
      "sigma_2": 0.02,
      "sigma_3": 0.02,
      "sigma_4": 0.02,
      "sigma_5": 0.02,
      "sigma_6": null,
      "sigma_7": null,
      "sigma_8": null,
      "sigma_9": null,
      "sigma_10": null,
      "sigma_11": null,
      "sigma_12": null,
      "weight_1": 2500.0,
      "weight_2": 2500.0,
      "weight_3": 2500.0,
      "weight_4": 2500.0,
      "weight_5": 2500.0,
      "weight_6": null,
      "weight_7": null,
      "weight_8": null,
      "weight_9": null,
      "weight_10": null,
      "weight_11": null,
      "weight_12": null,
      "rms_deltas_1": 0.0155,
      "rms_deltas_2": 0.0155,
      "rms_deltas_3": 0.0155,
      "rms_deltas_4": 0.0155,
      "rms_deltas_5": 0.0155,
      "rms_deltas_6": null,
      "rms_deltas_7": null,
      "rms_deltas_8": null,
      "rms_deltas_9": null,
      "rms_deltas_10": null,
      "rms_deltas_11": null,
      "rms_deltas_12": null,
      "residual_1": 3.62,
      "residual_2": 3.62,
      "residual_3": 3.62,
      "residual_4": 3.62,
      "residual_5": 3.62,
      "residual_6": null,
      "residual_7": null,
      "residual_8": null,
      "residual_9": null,
      "residual_10": null,
      "residual_11": null,
      "residual_12": null,
      "i_seq_1": 13,
      "i_seq_2": 14,
      "i_seq_3": 15,
      "i_seq_4": 16,
      "i_seq_5": 19,
      "i_seq_6": null,
      "i_seq_7": null,
      "i_seq_8": null,
      "i_seq_9": null,
      "i_seq_10": null,
      "i_seq_11": null,
      "i_seq_12": null
    }
  ]
}
"""
# Now with labels rather than i_seqs
sample_geo_text_2 = """
Bond restraints: 3
Sorted by residual:
bond pdb=" N   GLY A   1 "
     pdb=" CA  GLY A   1 "
  ideal  model  delta    sigma   weight residual
  1.451  1.507 -0.056 1.60e-02 3.91e+03 1.23e+01
bond pdb=" CA  GLN A   4 "
     pdb=" C   GLN A   4 "
  ideal  model  delta    sigma   weight residual
  1.522  1.553 -0.030 1.18e-02 7.18e+03 6.53e+00
bond pdb=" N   GLN A   4 "
     pdb=" CA  GLN A   4 "
  ideal  model  delta    sigma   weight residual
  1.460  1.485 -0.025 1.17e-02 7.31e+03 4.40e+00

Bond angle restraints: 2
Sorted by residual:
angle pdb=" N   ASN A   3 "
      pdb=" CA  ASN A   3 "
      pdb=" C   ASN A   3 "
    ideal   model   delta    sigma   weight residual
   108.90  113.48   -4.58 1.63e+00 3.76e-01 7.90e+00
angle pdb=" N   GLN A   4 "
      pdb=" CA  GLN A   4 "
      pdb=" C   GLN A   4 "
    ideal   model   delta    sigma   weight residual
   108.02  111.93   -3.91 1.78e+00 3.16e-01 4.84e+00

Dihedral angle restraints: 3
  sinusoidal: 15
    harmonic: 7
Sorted by residual:
dihedral pdb=" CA  ASN A   3 "
         pdb=" C   ASN A   3 "
         pdb=" N   GLN A   4 "
         pdb=" CA  GLN A   4 "
    ideal   model   delta  harmonic     sigma   weight residual
   180.00  166.21   13.79     0      5.00e+00 4.00e-02 7.60e+00
dihedral pdb=" CB  GLN A   5 "
         pdb=" CG  GLN A   5 "
         pdb=" CD  GLN A   5 "
         pdb=" OE1 GLN A   5 "
    ideal   model   delta sinusoidal    sigma   weight residual
     0.00  -72.39   72.39     2      3.00e+01 1.11e-03 4.85e+00
dihedral pdb=" CB  GLN A   4 "
         pdb=" CG  GLN A   4 "
         pdb=" CD  GLN A   4 "
         pdb=" OE1 GLN A   4 "
    ideal   model   delta sinusoidal    sigma   weight residual
     0.00   54.08  -54.08     2      3.00e+01 1.11e-03 3.50e+00
Planarity restraints: 2
Sorted by residual:
                               delta    sigma   weight rms_deltas residual
plane pdb=" CB  TYR A   7 "   -0.006 2.00e-02 2.50e+03   9.66e-03 1.87e+00
      pdb=" CG  TYR A   7 "    0.022 2.00e-02 2.50e+03
      pdb=" CD1 TYR A   7 "   -0.008 2.00e-02 2.50e+03
      pdb=" CD2 TYR A   7 "   -0.004 2.00e-02 2.50e+03
      pdb=" CE1 TYR A   7 "    0.002 2.00e-02 2.50e+03
      pdb=" CE2 TYR A   7 "   -0.001 2.00e-02 2.50e+03
      pdb=" CZ  TYR A   7 "   -0.011 2.00e-02 2.50e+03
      pdb=" OH  TYR A   7 "    0.006 2.00e-02 2.50e+03
                               delta    sigma   weight rms_deltas residual
plane pdb=" CB  ASN A   2 "    0.006 2.00e-02 2.50e+03   1.19e-02 1.42e+00
      pdb=" CG  ASN A   2 "   -0.021 2.00e-02 2.50e+03
      pdb=" OD1 ASN A   2 "    0.008 2.00e-02 2.50e+03
      pdb=" ND2 ASN A   2 "    0.007 2.00e-02 2.50e+03
"""
test_data_json_2 = """
{
  "angle": [
    {
      "restraint_type": "angle",
      "id_str_0": "pdb=\\\" N   ASN A   3 \\\"",
      "id_str_1": "pdb=\\\" CA  ASN A   3 \\\"",
      "id_str_2": "pdb=\\\" C   ASN A   3 \\\"",
      "ideal": 108.9,
      "model": 113.48,
      "delta": -4.58,
      "sigma": 1.63,
      "weight": 0.376,
      "residual": 7.9
    },
    {
      "restraint_type": "angle",
      "id_str_0": "pdb=\\\" N   GLN A   4 \\\"",
      "id_str_1": "pdb=\\\" CA  GLN A   4 \\\"",
      "id_str_2": "pdb=\\\" C   GLN A   4 \\\"",
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
      "id_str_0": "pdb=\\\" N   GLY A   1 \\\"",
      "id_str_1": "pdb=\\\" CA  GLY A   1 \\\"",
      "ideal": 1.451,
      "model": 1.507,
      "delta": -0.056,
      "sigma": 0.016,
      "weight": 3910.0,
      "residual": 12.3
    },
    {
      "restraint_type": "bond",
      "id_str_0": "pdb=\\\" CA  GLN A   4 \\\"",
      "id_str_1": "pdb=\\\" C   GLN A   4 \\\"",
      "ideal": 1.522,
      "model": 1.553,
      "delta": -0.03,
      "sigma": 0.0118,
      "weight": 7180.0,
      "residual": 6.53
    },
    {
      "restraint_type": "bond",
      "id_str_0": "pdb=\\\" N   GLN A   4 \\\"",
      "id_str_1": "pdb=\\\" CA  GLN A   4 \\\"",
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
      "id_str_0": "pdb=\\\" CA  ASN A   3 \\\"",
      "id_str_1": "pdb=\\\" C   ASN A   3 \\\"",
      "id_str_2": "pdb=\\\" N   GLN A   4 \\\"",
      "id_str_3": "pdb=\\\" CA  GLN A   4 \\\"",
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
      "id_str_0": "pdb=\\\" CB  GLN A   5 \\\"",
      "id_str_1": "pdb=\\\" CG  GLN A   5 \\\"",
      "id_str_2": "pdb=\\\" CD  GLN A   5 \\\"",
      "id_str_3": "pdb=\\\" OE1 GLN A   5 \\\"",
      "ideal": 0.0,
      "model": -72.39,
      "delta": 72.39,
      "harmonic": null,
      "sigma": 30.0,
      "weight": 0.00111,
      "residual": 4.85,
      "sinusoidal": 2.0
    },
    {
      "restraint_type": "dihedral",
      "id_str_0": "pdb=\\\" CB  GLN A   4 \\\"",
      "id_str_1": "pdb=\\\" CG  GLN A   4 \\\"",
      "id_str_2": "pdb=\\\" CD  GLN A   4 \\\"",
      "id_str_3": "pdb=\\\" OE1 GLN A   4 \\\"",
      "ideal": 0.0,
      "model": 54.08,
      "delta": -54.08,
      "harmonic": null,
      "sigma": 30.0,
      "weight": 0.00111,
      "residual": 3.5,
      "sinusoidal": 2.0
    }
  ],
  "plane": [
    {
      "delta_1": -0.006,
      "delta_2": 0.022,
      "delta_3": -0.008,
      "delta_4": -0.004,
      "delta_5": 0.002,
      "delta_6": -0.001,
      "delta_7": -0.011,
      "delta_8": 0.006,
      "sigma_1": 0.02,
      "sigma_2": 0.02,
      "sigma_3": 0.02,
      "sigma_4": 0.02,
      "sigma_5": 0.02,
      "sigma_6": 0.02,
      "sigma_7": 0.02,
      "sigma_8": 0.02,
      "weight_1": 2500.0,
      "weight_2": 2500.0,
      "weight_3": 2500.0,
      "weight_4": 2500.0,
      "weight_5": 2500.0,
      "weight_6": 2500.0,
      "weight_7": 2500.0,
      "weight_8": 2500.0,
      "rms_deltas_1": 0.00966,
      "rms_deltas_2": 0.00966,
      "rms_deltas_3": 0.00966,
      "rms_deltas_4": 0.00966,
      "rms_deltas_5": 0.00966,
      "rms_deltas_6": 0.00966,
      "rms_deltas_7": 0.00966,
      "rms_deltas_8": 0.00966,
      "residual_1": 1.87,
      "residual_2": 1.87,
      "residual_3": 1.87,
      "residual_4": 1.87,
      "residual_5": 1.87,
      "residual_6": 1.87,
      "residual_7": 1.87,
      "residual_8": 1.87,
      "id_str_1": "pdb=\\\" CB  TYR A   7 \\\"",
      "id_str_2": "pdb=\\\" CG  TYR A   7 \\\"",
      "id_str_3": "pdb=\\\" CD1 TYR A   7 \\\"",
      "id_str_4": "pdb=\\\" CD2 TYR A   7 \\\"",
      "id_str_5": "pdb=\\\" CE1 TYR A   7 \\\"",
      "id_str_6": "pdb=\\\" CE2 TYR A   7 \\\"",
      "id_str_7": "pdb=\\\" CZ  TYR A   7 \\\"",
      "id_str_8": "pdb=\\\" OH  TYR A   7 \\\""
    },
    {
      "delta_1": 0.006,
      "delta_2": -0.021,
      "delta_3": 0.008,
      "delta_4": null,
      "delta_5": null,
      "delta_6": null,
      "delta_7": null,
      "delta_8": null,
      "sigma_1": 0.02,
      "sigma_2": 0.02,
      "sigma_3": 0.02,
      "sigma_4": null,
      "sigma_5": null,
      "sigma_6": null,
      "sigma_7": null,
      "sigma_8": null,
      "weight_1": 2500.0,
      "weight_2": 2500.0,
      "weight_3": 2500.0,
      "weight_4": null,
      "weight_5": null,
      "weight_6": null,
      "weight_7": null,
      "weight_8": null,
      "rms_deltas_1": 0.0119,
      "rms_deltas_2": 0.0119,
      "rms_deltas_3": 0.0119,
      "rms_deltas_4": null,
      "rms_deltas_5": null,
      "rms_deltas_6": null,
      "rms_deltas_7": null,
      "rms_deltas_8": null,
      "residual_1": 1.42,
      "residual_2": 1.42,
      "residual_3": 1.42,
      "residual_4": null,
      "residual_5": null,
      "residual_6": null,
      "residual_7": null,
      "residual_8": null,
      "id_str_1": "pdb=\\\" CB  ASN A   2 \\\"",
      "id_str_2": "pdb=\\\" CG  ASN A   2 \\\"",
      "id_str_3": "pdb=\\\" OD1 ASN A   2 \\\"",
      "id_str_4": null,
      "id_str_5": null,
      "id_str_6": null,
      "id_str_7": null,
      "id_str_8": null
    }
  ]
}
"""


test_data = json.loads(test_data_json)
test_data_2 = json.loads(test_data_json_2)

# load file
with open("tst_geo_parsing.geo","w") as fh:
  fh.write(sample_geo_text)

# load file
with open("tst_geo_parsing_2.geo","w") as fh:
  fh.write(sample_geo_text_2)

# read back in
geo_dict = parse_geo_file("tst_geo_parsing.geo",return_format='dict')
geo_dict_2 = parse_geo_file("tst_geo_parsing_2.geo",return_format='dict')


# define comparison functions

def get_type(value):
  if pd.isna(value):
    return None
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

    # get all none-like values to be None
    if t1 is None:
      value1 = None
    if t2 is None:
      value2 = None
    assert str(value1)==str(value2), "Invalid comparisons: "+str(value1)+" and "+str(value2)
  return True



# test
def test1():
  print("Test1")
  for key,d_list in geo_dict.items():
    for j,d in enumerate(d_list):
      d_ref = test_data[key][j]
      same = compare_custom(d,d_ref)
      print(same)

def test2():
  print("Test2")
  for key,d_list in geo_dict_2.items():
    for j,d in enumerate(d_list):
      d_ref = test_data_2[key][j]
      same = compare_custom(d,d_ref)
      print(same)

if __name__ == '__main__':
  test1()
  test2()
  print("Finished")

