from __future__ import division
import json
import re
from collections import defaultdict
from io import StringIO
from itertools import chain
from iotbx.data_manager import DataManager
from cctbx.crystal.tst_super_cell import pdb_str_1yjp
from mmtbx.geometry_restraints.geo_file_parsing import GeoParser
from cctbx import geometry_restraints
from libtbx.test_utils import show_diff

"""
Example usage:

# Have a processed model
grm = model.restraints_manager.geometry

# Define the usual atom labels
site_labels = model.get_xray_structure().scatterers().extract_labels()

# write to geo_string
buffer = StringIO()
grm.write_geo_file(model.get_sites_cart(),site_labels=site_labels,file_descriptor=buffer)
geo_str = buffer.getvalue()

# Parse geo str
geo_lines = geo_str.split("\n")
geo_container = GeoParser(geo_lines,model=model)


# access entries
entry = geo_container.entries["Dihedral angle"][0]

# access entry as a pure dict
entry.record

{'i_seqs': [13, 14, 20, 21],
 'atom_labels': ['pdb=" CA  ASN A   3 "',
  'pdb=" C   ASN A   3 "',
  'pdb=" N   GLN A   4 "',
  'pdb=" CA  GLN A   4 "'],
 'ideal': 180.0,
 'model': 166.21,
 'delta': 13.79,
 'harmonic': 0,
 'sigma': 5.0,
 'weight': 0.04,
 'residual': 7.6,
 'origin_id': 0}


"""
tst_1_geo = """
# Geometry restraints


Bond | covalent geometry | restraints: 1
Sorted by residual:
bond pdb=" C1' C     643 " segid="D16S"
     pdb=" N1  C     643 " segid="D16S"
  ideal  model  delta    sigma   weight residual
  1.480  1.481 -0.001 1.50e-02 4.44e+03 8.32e-03

Bond | Misc. | restraints: 1
Sorted by residual:
bond pdb=" OP1 C     643 " segid="D16S"
     pdb=" N   SER   115 " segid="DS08"
  ideal  model  delta    sigma   weight residual
  1.297  1.297 -0.000 1.00e-02 1.00e+04 2.18e-14

Bond | link_ALPHA2-6 | restraints: 1
Sorted by residual:
bond pdb=" O6  BMA B1954 "
     pdb=" C2  MAN B1957 "
  ideal  model  delta    sigma   weight residual
  1.439  1.439 -0.000 2.00e-02 2.50e+03 2.41e-15

Bond | Disulphide bridge | restraints: 1
Sorted by residual:
bond pdb=" SG  CYS A  27 "
     pdb=" SG  CYS A 123 "
  ideal  model  delta    sigma   weight residual
  2.031  2.031  0.000 2.00e-02 2.50e+03 8.59e-17

Bond | link_NAG-ASN | restraints: 1
Sorted by residual:
bond pdb=" ND2 ASN A 303 "
     pdb=" C1  NAG B1349 "
  ideal  model  delta    sigma   weight residual
  1.439  1.441 -0.002 2.00e-02 2.50e+03 8.66e-03

Bond | Bond-like | restraints: 1
Sorted by residual:
bond pdb=" O2   DC A  10 "
     pdb=" N2   DG B  11 "
  ideal  model  delta    sigma   weight residual
  2.780  2.873 -0.093 1.00e-01 1.00e+02 8.58e-01

Bond | Metal coordination | restraints: 1
Sorted by residual:
bond pdb="MG   MG  A   1 " segid="A   "
     pdb=" O   HOH A   2 " segid="A   "
  ideal  model  delta    sigma   weight residual sym.op.
  2.070  2.070  0.000 5.00e-02 4.00e+02 1.18e-17 -x-1/2,y+1/2,-z+3/4

Bond | link_TRANS | restraints: 1
Sorted by residual:
bond pdb=" N   SER A   1 "
     pdb=" C   GLY A  34 "
  ideal  model  delta    sigma   weight residual
  1.329  1.329 -0.000 1.40e-02 5.10e+03 3.00e-05

Bond | Custom Glycosidic | restraints: 1
Sorted by residual:
bond pdb=" O6  NAG D1584 "
     pdb=" C1  FU4 D1588 "
  ideal  model  delta    sigma   weight residual
  1.439  1.439 -0.000 2.00e-02 2.50e+03 1.04e-10

Bond | link_BETA1-6 | restraints: 1
Sorted by residual:
bond pdb=" O6  NAG A 294 "
     pdb=" C1  FUC A 299 "
  ideal  model  delta    sigma   weight residual
  1.439  1.443 -0.004 2.00e-02 2.50e+03 3.32e-02

Bond | link_BETA1-3 | restraints: 1
Sorted by residual:
bond pdb=" O3  BGC B 401 "
     pdb=" C1  BGC B 402 "
  ideal  model  delta    sigma   weight residual
  1.439  1.439  0.000 2.00e-02 2.50e+03 3.54e-16

Bond | User supplied | restraints: 1
Sorted by residual:
bond pdb=" SG  CYS A 470 "
     pdb="FE   HEM A 601 "
  ideal  model  delta    sigma   weight residual
  2.330  2.330  0.000 2.00e-02 2.50e+03 3.62e-12

Bond | link_BETA1-4 | restraints: 1
Sorted by residual:
bond pdb=" O4  NAG A 361 "
     pdb=" C1  NAG A 362 "
  ideal  model  delta    sigma   weight residual
  1.439  1.442 -0.003 2.00e-02 2.50e+03 2.84e-02

Bond angle | covalent geometry | restraints: 1
Sorted by residual:
angle pdb=" C2' C     643 " segid="D16S"
      pdb=" C1' C     643 " segid="D16S"
      pdb=" N1  C     643 " segid="D16S"
    ideal   model   delta    sigma   weight residual
   112.00  112.26   -0.26 1.50e+00 4.44e-01 3.04e-02

Bond angle | link_ALPHA2-6 | restraints: 1
Sorted by residual:
angle pdb=" O6  BMA B1954 "
      pdb=" C2  MAN B1957 "
      pdb=" C1  MAN B1957 "
    ideal   model   delta    sigma   weight residual
   112.30  112.30    0.00 3.00e+00 1.11e-01 8.81e-13

Bond angle | Disulphide bridge | restraints: 1
Sorted by residual:
angle pdb=" SG  CYS A  27 "
      pdb=" SG  CYS A 123 "
      pdb=" CB  CYS A 123 "
    ideal   model   delta    sigma   weight residual
   104.20  104.20    0.00 2.10e+00 2.27e-01 1.59e-15

Bond angle | link_NAG-ASN | restraints: 1
Sorted by residual:
angle pdb=" ND2 ASN A 303 "
      pdb=" C1  NAG B1349 "
      pdb=" O5  NAG B1349 "
    ideal   model   delta    sigma   weight residual
   112.30  111.42    0.88 3.00e+00 1.11e-01 8.70e-02

Bond angle | Secondary Structure restraints around h-bond | restraints: 1
Sorted by residual:
angle pdb=" C4   DC A  10 "
      pdb=" N4   DC A  10 "
      pdb=" O6   DG B  11 "
    ideal   model   delta    sigma   weight residual
   117.30  120.06   -2.76 2.86e+00 1.22e-01 9.34e-01

Bond angle | link_TRANS | restraints: 1
Sorted by residual:
angle pdb=" CA  SER A   1 "
      pdb=" N   SER A   1 "
      pdb=" C   GLY A  34 "
    ideal   model   delta    sigma   weight residual
   121.70  121.74   -0.04 1.80e+00 3.09e-01 5.23e-04

Bond angle | Custom Glycosidic | restraints: 1
Sorted by residual:
angle pdb=" O6  NAG D1584 "
      pdb=" C1  FU4 D1588 "
      pdb=" C2  FU4 D1588 "
    ideal   model   delta    sigma   weight residual
   109.47  109.47   -0.00 3.00e+00 1.11e-01 2.59e-10

Bond angle | link_BETA1-6 | restraints: 1
Sorted by residual:
angle pdb=" O6  NAG A 294 "
      pdb=" C1  FUC A 299 "
      pdb=" O5  FUC A 299 "
    ideal   model   delta    sigma   weight residual
   112.30  110.82    1.48 3.00e+00 1.11e-01 2.43e-01

Bond angle | link_BETA1-3 | restraints: 1
Sorted by residual:
angle pdb=" O3  BGC B 401 "
      pdb=" C1  BGC B 402 "
      pdb=" O5  BGC B 402 "
    ideal   model   delta    sigma   weight residual
   112.30  112.30    0.00 3.00e+00 1.11e-01 2.31e-14

Bond angle | link_BETA1-4 | restraints: 1
Sorted by residual:
angle pdb=" O4  NAG A 361 "
      pdb=" C1  NAG A 362 "
      pdb=" O5  NAG A 362 "
    ideal   model   delta    sigma   weight residual
   112.30  111.17    1.13 3.00e+00 1.11e-01 1.41e-01

Dihedral angle | covalent geometry | restraints: 1
  sinusoidal: 13
    harmonic: 0
Sorted by residual:
dihedral pdb=" C1' C     643 " segid="D16S"
         pdb=" N1  C     643 " segid="D16S"
         pdb=" C6  C     643 " segid="D16S"
         pdb=" C5  C     643 " segid="D16S"
    ideal   model   delta sinusoidal    sigma   weight residual
  -106.54 -179.87   73.33     2      3.00e+01 1.11e-03 4.89e+00

Dihedral angle | C-Beta improper | restraints: 1
Sorted by residual:
dihedral pdb=" C   SER   115 " segid="DS08"
         pdb=" N   SER   115 " segid="DS08"
         pdb=" CA  SER   115 " segid="DS08"
         pdb=" CB  SER   115 " segid="DS08"
    ideal   model   delta  harmonic     sigma   weight residual
  -122.60 -122.55   -0.05     0      2.50e+00 1.60e-01 4.55e-04

Dihedral angle | Side chain | restraints: 1
Sorted by residual:
dihedral pdb=" N   SER   115 " segid="DS08"
         pdb=" CA  SER   115 " segid="DS08"
         pdb=" CB  SER   115 " segid="DS08"
         pdb=" OG  SER   115 " segid="DS08"
    ideal   model   delta  harmonic     sigma   weight residual
   -63.99  -63.99   -0.00     0      1.00e+01 1.00e-02 4.15e-11

Dihedral angle | link_TRANS | restraints: 1
Sorted by residual:
dihedral pdb=" N   SER A   1 "
         pdb=" CA  GLY A  34 "
         pdb=" C   GLY A  34 "
         pdb=" N   GLY A  34 "
    ideal   model   delta sinusoidal    sigma   weight residual
  -160.00 -145.03  -14.97     2      3.00e+01 1.11e-03 3.56e-01

Chirality | covalent geometry | restraints: 1
Sorted by residual:
chirality pdb=" C1' C     643 " segid="D16S"
          pdb=" O4' C     643 " segid="D16S"
          pdb=" C2' C     643 " segid="D16S"
          pdb=" N1  C     643 " segid="D16S"
  both_signs  ideal   model   delta    sigma   weight residual
    False      2.47    2.45    0.01 2.00e-01 2.50e+01 4.38e-03

Chirality | link_ALPHA2-6 | restraints: 1
Sorted by residual:
chirality pdb=" C2  MAN B1957 "
          pdb=" O6  BMA B1954 "
          pdb=" C1  MAN B1957 "
          pdb=" C1  MAN B1957 "
  both_signs  ideal   model   delta    sigma   weight residual
    False     -2.40   -0.00   -2.40 2.00e-02 2.50e+03 1.44e+04

Chirality | link_NAG-ASN | restraints: 1
Sorted by residual:
chirality pdb=" C1  NAG B1349 "
          pdb=" ND2 ASN A 303 "
          pdb=" C2  NAG B1349 "
          pdb=" O5  NAG B1349 "
  both_signs  ideal   model   delta    sigma   weight residual
    False     -2.40   -2.28   -0.12 2.00e-01 2.50e+01 3.45e-01

Chirality | link_BETA1-6 | restraints: 1
Sorted by residual:
chirality pdb=" C1  FUC A 299 "
          pdb=" O6  NAG A 294 "
          pdb=" C2  FUC A 299 "
          pdb=" O5  FUC A 299 "
  both_signs  ideal   model   delta    sigma   weight residual
    False     -2.40   -2.40   -0.00 2.00e-02 2.50e+03 1.20e-02

Chirality | link_BETA1-3 | restraints: 1
Sorted by residual:
chirality pdb=" C1  BGC B 402 "
          pdb=" O3  BGC B 401 "
          pdb=" C2  BGC B 402 "
          pdb=" O5  BGC B 402 "
  both_signs  ideal   model   delta    sigma   weight residual
    False     -2.40   -2.40   -0.00 2.00e-02 2.50e+03 2.55e-18

Chirality | link_BETA1-4 | restraints: 1
Sorted by residual:
chirality pdb=" C1  NAG A 362 "
          pdb=" O4  NAG A 361 "
          pdb=" C2  NAG A 362 "
          pdb=" O5  NAG A 362 "
  both_signs  ideal   model   delta    sigma   weight residual
    False     -2.40   -2.40   -0.00 2.00e-02 2.50e+03 7.20e-03

Planarity | covalent geometry | restraints: 1
Sorted by residual:
                                            delta    sigma   weight rms_deltas residual
plane pdb=" C1' C     643 " segid="D16S"   -0.001 2.00e-02 2.50e+03   6.70e-04 1.01e-02
      pdb=" N1  C     643 " segid="D16S"    0.002 2.00e-02 2.50e+03
      pdb=" C2  C     643 " segid="D16S"    0.000 2.00e-02 2.50e+03
      pdb=" O2  C     643 " segid="D16S"   -0.000 2.00e-02 2.50e+03
      pdb=" N3  C     643 " segid="D16S"    0.000 2.00e-02 2.50e+03
      pdb=" C4  C     643 " segid="D16S"    0.000 2.00e-02 2.50e+03
      pdb=" N4  C     643 " segid="D16S"    0.000 2.00e-02 2.50e+03
      pdb=" C5  C     643 " segid="D16S"   -0.001 2.00e-02 2.50e+03
      pdb=" C6  C     643 " segid="D16S"   -0.000 2.00e-02 2.50e+03

Planarity | link_NAG-ASN | restraints: 1
Sorted by residual:
                               delta    sigma   weight rms_deltas residual
plane pdb=" CB  ASN A 303 "    0.000 2.00e-02 2.50e+03   4.37e-10 2.39e-15
      pdb=" CG  ASN A 303 "   -0.000 2.00e-02 2.50e+03
      pdb=" OD1 ASN A 303 "   -0.000 2.00e-02 2.50e+03
      pdb=" ND2 ASN A 303 "   -0.000 2.00e-02 2.50e+03
      pdb=" C1  NAG B1349 "    0.000 2.00e-02 2.50e+03

Planarity | link_TRANS | restraints: 1
Sorted by residual:
                               delta    sigma   weight rms_deltas residual
plane pdb=" N   SER A   1 "   -0.000 2.00e-02 2.50e+03   6.43e-04 4.14e-03
      pdb=" CA  GLY A  34 "   -0.000 2.00e-02 2.50e+03
      pdb=" C   GLY A  34 "    0.001 2.00e-02 2.50e+03
      pdb=" O   GLY A  34 "   -0.000 2.00e-02 2.50e+03

Nonbonded | unspecified | interactions: 1
Sorted by model distance:
nonbonded pdb=" O4' C     643 " segid="D16S"
          pdb=" C6  C     643 " segid="D16S"
   model   vdw
   2.683 2.672

Parallelity | Stacking parallelity | restraints: 1
Sorted by residual:
    plane 1                plane 2                residual  delta(deg) sigma
    pdb=" C1'  DG B  11 "  pdb=" C1'  DC B  12 "  6.47e+00   5.5671    0.0270
    pdb=" N9   DG B  11 "  pdb=" N1   DC B  12 "
    pdb=" C8   DG B  11 "  pdb=" C2   DC B  12 "
    pdb=" N7   DG B  11 "  pdb=" O2   DC B  12 "
    pdb=" C5   DG B  11 "  pdb=" N3   DC B  12 "
    pdb=" C6   DG B  11 "  pdb=" C4   DC B  12 "
    pdb=" O6   DG B  11 "  pdb=" N4   DC B  12 "
    pdb=" N1   DG B  11 "  pdb=" C5   DC B  12 "
    pdb=" C2   DG B  11 "  pdb=" C6   DC B  12 "
    pdb=" N2   DG B  11 "
    pdb=" N3   DG B  11 "
    pdb=" C4   DG B  11 "

Parallelity | Basepair parallelity | restraints: 1
Sorted by residual:
    plane 1                plane 2                residual  delta(deg) sigma
    pdb=" C1'  DG A   9 "  pdb=" C1'  DC B  12 "  2.26e+01  12.9202    0.0335
    pdb=" N9   DG A   9 "  pdb=" N1   DC B  12 "
    pdb=" C8   DG A   9 "  pdb=" C2   DC B  12 "
    pdb=" N7   DG A   9 "  pdb=" O2   DC B  12 "
    pdb=" C5   DG A   9 "  pdb=" N3   DC B  12 "
    pdb=" C6   DG A   9 "  pdb=" C4   DC B  12 "
    pdb=" O6   DG A   9 "  pdb=" N4   DC B  12 "
    pdb=" N1   DG A   9 "  pdb=" C5   DC B  12 "
    pdb=" C2   DG A   9 "  pdb=" C6   DC B  12 "
    pdb=" N2   DG A   9 "
    pdb=" N3   DG A   9 "
    pdb=" C4   DG A   9 "

"""

tst_2_geo = """
# Geometry restraints

Bond | covalent geometry | restraints: 4
Sorted by residual:
bond 0
     1
  ideal  model  delta    sigma   weight residual
  1.451  1.507 -0.056 1.60e-02 3.91e+03 1.23e+01
bond 21
     22
  ideal  model  delta    sigma   weight residual
  1.522  1.553 -0.030 1.18e-02 7.18e+03 6.53e+00
bond 20
     21
  ideal  model  delta    sigma   weight residual
  1.460  1.485 -0.025 1.17e-02 7.31e+03 4.40e+00
bond 5
     6
  ideal  model  delta    sigma   weight residual
  1.524  1.498  0.025 1.26e-02 6.30e+03 4.00e+00
"""
# RESULT STRINGS
# Json strings that represent entry records.
# This is the structured data we expect to parse from the .geo text

# 1yjp .geo text (with labels), write the first entry of each type to json
result_json_1 = """
{
  "Bond": [
    {
      "i_seqs": [
        0,
        1
      ],
      "atom_labels": [
        "pdb=\\\" N   GLY A   1 \\\"",
        "pdb=\\\" CA  GLY A   1 \\\""
      ],
      "ideal": 1.451,
      "model": 1.507,
      "delta": -0.056,
      "sigma": 0.016,
      "weight": 3910.0,
      "residual": 12.3,
      "origin_id": 0
    }
  ],
  "Bond angle": [
    {
      "i_seqs": [
        12,
        13,
        14
      ],
      "atom_labels": [
        "pdb=\\\" N   ASN A   3 \\\"",
        "pdb=\\\" CA  ASN A   3 \\\"",
        "pdb=\\\" C   ASN A   3 \\\""
      ],
      "ideal": 108.9,
      "model": 113.48,
      "delta": -4.58,
      "sigma": 1.63,
      "weight": 0.376,
      "residual": 7.9,
      "origin_id": 0
    }
  ],
  "Dihedral angle": [
    {
      "i_seqs": [
        13,
        14,
        20,
        21
      ],
      "atom_labels": [
        "pdb=\\\" CA  ASN A   3 \\\"",
        "pdb=\\\" C   ASN A   3 \\\"",
        "pdb=\\\" N   GLN A   4 \\\"",
        "pdb=\\\" CA  GLN A   4 \\\""
      ],
      "ideal": 180.0,
      "model": 166.21,
      "delta": 13.79,
      "harmonic": 0,
      "sigma": 5.0,
      "weight": 0.04,
      "residual": 7.6,
      "origin_id": 0
    }
  ],
  "Chirality": [
    {
      "i_seqs": [
        30,
        29,
        31,
        33
      ],
      "atom_labels": [
        "pdb=\\\" CA  GLN A   5 \\\"",
        "pdb=\\\" N   GLN A   5 \\\"",
        "pdb=\\\" C   GLN A   5 \\\"",
        "pdb=\\\" CB  GLN A   5 \\\""
      ],
      "both_signs": "False",
      "ideal": 2.51,
      "model": 2.39,
      "delta": 0.12,
      "sigma": 0.2,
      "weight": 25.0,
      "residual": 0.348,
      "origin_id": 0
    }
  ],
  "Planarity": [
    {
      "i_seqs": [
        50,
        51,
        52,
        53,
        54,
        55,
        56,
        57
      ],
      "atom_labels": [
        "pdb=\\\" CB  TYR A   7 \\\"",
        "pdb=\\\" CG  TYR A   7 \\\"",
        "pdb=\\\" CD1 TYR A   7 \\\"",
        "pdb=\\\" CD2 TYR A   7 \\\"",
        "pdb=\\\" CE1 TYR A   7 \\\"",
        "pdb=\\\" CE2 TYR A   7 \\\"",
        "pdb=\\\" CZ  TYR A   7 \\\"",
        "pdb=\\\" OH  TYR A   7 \\\""
      ],
      "delta": [
        -0.006,
        0.022,
        -0.008,
        -0.004,
        0.002,
        -0.001,
        -0.011,
        0.006
      ],
      "sigma": [
        0.02,
        0.02,
        0.02,
        0.02,
        0.02,
        0.02,
        0.02,
        0.02
      ],
      "weight": [
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0
      ],
      "rms_deltas": [
        0.00966,
        0.00966,
        0.00966,
        0.00966,
        0.00966,
        0.00966,
        0.00966,
        0.00966
      ],
      "residual": [
        1.87,
        1.87,
        1.87,
        1.87,
        1.87,
        1.87,
        1.87,
        1.87
      ],
      "origin_id": 0
    }
  ],
  "Nonbonded": [
    {
      "i_seqs": [
        57,
        62
      ],
      "atom_labels": [
        "pdb=\\\" OH  TYR A   7 \\\"",
        "pdb=\\\" O   HOH A  11 \\\""
      ],
      "model": 2.525,
      "vdw": 3.04,
      "sym.op.": "-x+1,y-1/2,-z+1",
      "origin_id": 0
    }
  ]
}
"""
# 1yjp .geo text (without labels), write the first entry of each type to json
result_json_2 = """
{
  "Bond": [
    {
      "i_seqs": [
        0,
        1
      ],
      "atom_labels": [],
      "ideal": 1.451,
      "model": 1.507,
      "delta": -0.056,
      "sigma": 0.016,
      "weight": 3910.0,
      "residual": 12.3,
      "origin_id": 0
    }
  ],
  "Bond angle": [
    {
      "i_seqs": [
        12,
        13,
        14
      ],
      "atom_labels": [],
      "ideal": 108.9,
      "model": 113.48,
      "delta": -4.58,
      "sigma": 1.63,
      "weight": 0.376,
      "residual": 7.9,
      "origin_id": 0
    }
  ],
  "Dihedral angle": [
    {
      "i_seqs": [
        13,
        14,
        20,
        21
      ],
      "atom_labels": [],
      "ideal": 180.0,
      "model": 166.21,
      "delta": 13.79,
      "harmonic": 0,
      "sigma": 5.0,
      "weight": 0.04,
      "residual": 7.6,
      "origin_id": 0
    }
  ],
  "Chirality": [
    {
      "i_seqs": [
        30,
        29,
        31,
        33
      ],
      "atom_labels": [],
      "both_signs": "False",
      "ideal": 2.51,
      "model": 2.39,
      "delta": 0.12,
      "sigma": 0.2,
      "weight": 25.0,
      "residual": 0.348,
      "origin_id": 0
    }
  ],
  "Planarity": [
    {
      "i_seqs": [
        50,
        51,
        52,
        53,
        54,
        55,
        56,
        57
      ],
      "atom_labels": [],
      "delta": [
        -0.006,
        0.022,
        -0.008,
        -0.004,
        0.002,
        -0.001,
        -0.011,
        0.006
      ],
      "sigma": [
        0.02,
        0.02,
        0.02,
        0.02,
        0.02,
        0.02,
        0.02,
        0.02
      ],
      "weight": [
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0
      ],
      "rms_deltas": [
        0.00966,
        0.00966,
        0.00966,
        0.00966,
        0.00966,
        0.00966,
        0.00966,
        0.00966
      ],
      "residual": [
        1.87,
        1.87,
        1.87,
        1.87,
        1.87,
        1.87,
        1.87,
        1.87
      ],
      "origin_id": 0
    }
  ],
  "Nonbonded": [
    {
      "i_seqs": [
        57,
        62
      ],
      "atom_labels": [],
      "model": 2.525,
      "vdw": 3.04,
      "sym.op.": "-x+1,y-1/2,-z+1",
      "origin_id": 0
    }
  ]
}
"""
# 1yjp .geo text (with labels but no model), So no i_seqs possible
result_json_3 = """
{
  "Bond": [
    {
      "i_seqs": [],
      "atom_labels": [
        "pdb=\\\" N   GLY A   1 \\\"",
        "pdb=\\\" CA  GLY A   1 \\\""
      ],
      "ideal": 1.451,
      "model": 1.507,
      "delta": -0.056,
      "sigma": 0.016,
      "weight": 3910.0,
      "residual": 12.3,
      "origin_id": 0
    }
  ],
  "Bond angle": [
    {
      "i_seqs": [],
      "atom_labels": [
        "pdb=\\\" N   ASN A   3 \\\"",
        "pdb=\\\" CA  ASN A   3 \\\"",
        "pdb=\\\" C   ASN A   3 \\\""
      ],
      "ideal": 108.9,
      "model": 113.48,
      "delta": -4.58,
      "sigma": 1.63,
      "weight": 0.376,
      "residual": 7.9,
      "origin_id": 0
    }
  ],
  "Dihedral angle": [
    {
      "i_seqs": [],
      "atom_labels": [
        "pdb=\\\" CA  ASN A   3 \\\"",
        "pdb=\\\" C   ASN A   3 \\\"",
        "pdb=\\\" N   GLN A   4 \\\"",
        "pdb=\\\" CA  GLN A   4 \\\""
      ],
      "ideal": 180.0,
      "model": 166.21,
      "delta": 13.79,
      "harmonic": 0,
      "sigma": 5.0,
      "weight": 0.04,
      "residual": 7.6,
      "origin_id": 0
    }
  ],
  "Chirality": [
    {
      "i_seqs": [],
      "atom_labels": [
        "pdb=\\\" CA  GLN A   5 \\\"",
        "pdb=\\\" N   GLN A   5 \\\"",
        "pdb=\\\" C   GLN A   5 \\\"",
        "pdb=\\\" CB  GLN A   5 \\\""
      ],
      "both_signs": "False",
      "ideal": 2.51,
      "model": 2.39,
      "delta": 0.12,
      "sigma": 0.2,
      "weight": 25.0,
      "residual": 0.348,
      "origin_id": 0
    }
  ],
  "Planarity": [
    {
      "i_seqs": [],
      "atom_labels": [
        "pdb=\\\" CB  TYR A   7 \\\"",
        "pdb=\\\" CG  TYR A   7 \\\"",
        "pdb=\\\" CD1 TYR A   7 \\\"",
        "pdb=\\\" CD2 TYR A   7 \\\"",
        "pdb=\\\" CE1 TYR A   7 \\\"",
        "pdb=\\\" CE2 TYR A   7 \\\"",
        "pdb=\\\" CZ  TYR A   7 \\\"",
        "pdb=\\\" OH  TYR A   7 \\\""
      ],
      "delta": [
        -0.006,
        0.022,
        -0.008,
        -0.004,
        0.002,
        -0.001,
        -0.011,
        0.006
      ],
      "sigma": [
        0.02,
        0.02,
        0.02,
        0.02,
        0.02,
        0.02,
        0.02,
        0.02
      ],
      "weight": [
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0
      ],
      "rms_deltas": [
        0.00966,
        0.00966,
        0.00966,
        0.00966,
        0.00966,
        0.00966,
        0.00966,
        0.00966
      ],
      "residual": [
        1.87,
        1.87,
        1.87,
        1.87,
        1.87,
        1.87,
        1.87,
        1.87
      ],
      "origin_id": 0
    }
  ],
  "Nonbonded": [
    {
      "i_seqs": [],
      "atom_labels": [
        "pdb=\\\" OH  TYR A   7 \\\"",
        "pdb=\\\" O   HOH A  11 \\\""
      ],
      "model": 2.525,
      "vdw": 3.04,
      "sym.op.": "-x+1,y-1/2,-z+1",
      "origin_id": 0
    }
  ]
}
"""
# tst_1_geo parsed (testing multiple types of origin ids), the first entry for each entry type.
result_json_4 = """
{
  "Bond": [
    {
      "i_seqs": [],
      "atom_labels": [
        "pdb=\\\" C1' C     643 \\\" segid=\\\"D16S\\\"",
        "pdb=\\\" N1  C     643 \\\" segid=\\\"D16S\\\""
      ],
      "ideal": 1.48,
      "model": 1.481,
      "delta": -0.001,
      "sigma": 0.015,
      "weight": 4440.0,
      "residual": 0.00832,
      "origin_id": 0
    }
  ],
  "Bond angle": [
    {
      "i_seqs": [],
      "atom_labels": [
        "pdb=\\\" C2' C     643 \\\" segid=\\\"D16S\\\"",
        "pdb=\\\" C1' C     643 \\\" segid=\\\"D16S\\\"",
        "pdb=\\\" N1  C     643 \\\" segid=\\\"D16S\\\""
      ],
      "ideal": 112.0,
      "model": 112.26,
      "delta": -0.26,
      "sigma": 1.5,
      "weight": 0.444,
      "residual": 0.0304,
      "origin_id": 0
    }
  ],
  "Dihedral angle": [
    {
      "i_seqs": [],
      "atom_labels": [
        "pdb=\\\" C1' C     643 \\\" segid=\\\"D16S\\\"",
        "pdb=\\\" N1  C     643 \\\" segid=\\\"D16S\\\"",
        "pdb=\\\" C6  C     643 \\\" segid=\\\"D16S\\\"",
        "pdb=\\\" C5  C     643 \\\" segid=\\\"D16S\\\""
      ],
      "ideal": -106.54,
      "model": -179.87,
      "delta": 73.33,
      "sinusoidal": 2,
      "sigma": 30.0,
      "weight": 0.00111,
      "residual": 4.89,
      "origin_id": 0
    }
  ],
  "Chirality": [
    {
      "i_seqs": [],
      "atom_labels": [
        "pdb=\\\" C1' C     643 \\\" segid=\\\"D16S\\\"",
        "pdb=\\\" O4' C     643 \\\" segid=\\\"D16S\\\"",
        "pdb=\\\" C2' C     643 \\\" segid=\\\"D16S\\\"",
        "pdb=\\\" N1  C     643 \\\" segid=\\\"D16S\\\""
      ],
      "both_signs": "False",
      "ideal": 2.47,
      "model": 2.45,
      "delta": 0.01,
      "sigma": 0.2,
      "weight": 25.0,
      "residual": 0.00438,
      "origin_id": 0
    }
  ],
  "Planarity": [
    {
      "i_seqs": [],
      "atom_labels": [
        "pdb=\\\" C1' C     643 \\\" segid=D16S",
        "pdb=\\\" N1  C     643 \\\" segid=D16S",
        "pdb=\\\" C2  C     643 \\\" segid=D16S",
        "pdb=\\\" O2  C     643 \\\" segid=D16S",
        "pdb=\\\" N3  C     643 \\\" segid=D16S",
        "pdb=\\\" C4  C     643 \\\" segid=D16S",
        "pdb=\\\" N4  C     643 \\\" segid=D16S",
        "pdb=\\\" C5  C     643 \\\" segid=D16S",
        "pdb=\\\" C6  C     643 \\\" segid=D16S"
      ],
      "delta": [
        -0.001,
        0.002,
        0.0,
        -0.0,
        0.0,
        0.0,
        0.0,
        -0.001,
        -0.0
      ],
      "sigma": [
        0.02,
        0.02,
        0.02,
        0.02,
        0.02,
        0.02,
        0.02,
        0.02,
        0.02
      ],
      "weight": [
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0
      ],
      "rms_deltas": [
        0.00067,
        0.00067,
        0.00067,
        0.00067,
        0.00067,
        0.00067,
        0.00067,
        0.00067,
        0.00067
      ],
      "residual": [
        0.0101,
        0.0101,
        0.0101,
        0.0101,
        0.0101,
        0.0101,
        0.0101,
        0.0101,
        0.0101
      ],
      "origin_id": 0
    }
  ],
  "Nonbonded": [
    {
      "i_seqs": [],
      "atom_labels": [
        "pdb=\\\" O4' C     643 \\\" segid=\\\"D16S\\\"",
        "pdb=\\\" C6  C     643 \\\" segid=\\\"D16S\\\""
      ],
      "model": 2.683,
      "vdw": 2.672,
      "origin_id": 0
    }
  ],
  "Parallelity": [
    {
      "i_seqs": [],
      "j_seqs": [],
      "atom_labels_i": [
        "pdb=\\\" C1'  DG B  11 \\\"",
        "pdb=\\\" N9   DG B  11 \\\"",
        "pdb=\\\" C8   DG B  11 \\\"",
        "pdb=\\\" N7   DG B  11 \\\"",
        "pdb=\\\" C5   DG B  11 \\\"",
        "pdb=\\\" C6   DG B  11 \\\"",
        "pdb=\\\" O6   DG B  11 \\\"",
        "pdb=\\\" N1   DG B  11 \\\"",
        "pdb=\\\" C2   DG B  11 \\\"",
        "pdb=\\\" N2   DG B  11 \\\"",
        "pdb=\\\" N3   DG B  11 \\\"",
        "pdb=\\\" C4   DG B  11 \\\""
      ],
      "atom_labels_j": [
        "pdb=\\\" C1'  DC B  12 \\\"",
        "pdb=\\\" N1   DC B  12 \\\"",
        "pdb=\\\" C2   DC B  12 \\\"",
        "pdb=\\\" O2   DC B  12 \\\"",
        "pdb=\\\" N3   DC B  12 \\\"",
        "pdb=\\\" C4   DC B  12 \\\"",
        "pdb=\\\" N4   DC B  12 \\\"",
        "pdb=\\\" C5   DC B  12 \\\"",
        "pdb=\\\" C6   DC B  12 \\\""
      ],
      "residual": 6.47,
      "delta(deg)": 5.5671,
      "sigma": 0.027,
      "origin_id": 6
    }
  ]
}
"""
# tst_1_geo but the descriptive labels replaced with artifical i_seqs.
# Tests ability to read atom labels in varied forms (id_str, i_seq)
result_json_5 = """
{
  "Bond": [
    {
      "i_seqs": [
        0,
        1
      ],
      "atom_labels": [],
      "ideal": 1.48,
      "model": 1.481,
      "delta": -0.001,
      "sigma": 0.015,
      "weight": 4440.0,
      "residual": 0.00832,
      "origin_id": 0
    }
  ],
  "Bond angle": [
    {
      "i_seqs": [
        26,
        27,
        28
      ],
      "atom_labels": [],
      "ideal": 112.0,
      "model": 112.26,
      "delta": -0.26,
      "sigma": 1.5,
      "weight": 0.444,
      "residual": 0.0304,
      "origin_id": 0
    }
  ],
  "Dihedral angle": [
    {
      "i_seqs": [
        25,
        26,
        27,
        28
      ],
      "atom_labels": [],
      "ideal": -106.54,
      "model": -179.87,
      "delta": 73.33,
      "sinusoidal": 2,
      "sigma": 30.0,
      "weight": 0.00111,
      "residual": 4.89,
      "origin_id": 0
    }
  ],
  "Chirality": [
    {
      "i_seqs": [
        10,
        11,
        12,
        13
      ],
      "atom_labels": [],
      "both_signs": "False",
      "ideal": 2.47,
      "model": 2.45,
      "delta": 0.01,
      "sigma": 0.2,
      "weight": 25.0,
      "residual": 0.00438,
      "origin_id": 0
    }
  ],
  "Planarity": [
    {
      "i_seqs": [
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11
      ],
      "atom_labels": [],
      "delta": [
        -0.001,
        0.002,
        0.000,
        -0.0,
        0.0,
        0.0,
        0.0,
        -0.001,
        -0.0
      ],
      "sigma": [
        0.02,
        0.02,
        0.02,
        0.02,
        0.02,
        0.02,
        0.02,
        0.02,
        0.02
      ],
      "weight": [
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0,
        2500.0
      ],
      "rms_deltas": [
        0.00067,
        0.00067,
        0.00067,
        0.00067,
        0.00067,
        0.00067,
        0.00067,
        0.00067,
        0.00067
      ],
      "residual": [
        0.0101,
        0.0101,
        0.0101,
        0.0101,
        0.0101,
        0.0101,
        0.0101,
        0.0101,
        0.0101
      ],
      "origin_id": 0
    }
  ],
  "Nonbonded": [
    {
      "i_seqs": [
        21,
        22
      ],
      "atom_labels": [],
      "model": 2.683,
      "vdw": 2.672,
      "origin_id": 0
    }
  ],
  "Parallelity": [
    {
      "i_seqs": [
        23,
        25,
        27,
        29,
        0,
        2,
        4,
        6,
        8,
        10,
        11,
        12
      ],
      "j_seqs": [
        24,
        26,
        28,
        30,
        1,
        3,
        5,
        7,
        9
      ],
      "atom_labels_i": [],
      "atom_labels_j": [],
      "residual": 6.47,
      "delta(deg)": 5.5671,
      "sigma": 0.027,
      "origin_id": 6
    }
  ]
}
"""
def replace_idstr_with_int(text, max_int=100):
  """
  Replace id_strs in a geo_file str with integers
    For testing
  """
  index = [0]  # Use a list to allow modification inside the replacement function

  def replacement(match):
    original_length = len(match.group(0))
    replacement_text = '{:d}'.format(index[0]).ljust(original_length)
    index[0] += 1
    if index[0] > max_int:
      index[0] = 0
    return replacement_text

  new_text = re.sub(r'pdb=".*?"', replacement, text)
  new_text = re.sub(r'segid=".*?"', "", new_text)
  return new_text


def extract_results(container,print_result=False):
  results = defaultdict(list)
  for entry_name,entries in container.entries.items():
    entry = entries[0]
    results[entry_name].append(entry.record)

  #Print output to make tests
  js = json.dumps(results,indent=2)
  js = js.replace('\\\"','\\\\\\"')
  if print_result:
    print(js)
  return results, js

# Start tests
def init_model():
  dm = DataManager()
  dm.process_model_str("1yjp",pdb_str_1yjp)
  model= dm.get_model()
  model.process(make_restraints=True)
  return model

def tst_01(model,printing=False):
  # Test a 1yjp with YES labels and YES a model
  # (Can build proxies from label matching to model i_seqs)
  expected= json.loads(result_json_1)

  grm = model.restraints_manager.geometry

  buffer = StringIO()
  site_labels = model.get_xray_structure().scatterers().extract_labels()
  grm.write_geo_file(model.get_sites_cart(),site_labels=site_labels,file_descriptor=buffer)
  geo_str = buffer.getvalue()
  geo_lines = geo_str.split("\n")
  geo_container = GeoParser(geo_lines,model=model)

  if printing:
    print("\n\ntst_01")
  results, str_results = extract_results(geo_container,print_result=printing)

  # Check values
  if not expected==results:
    # show_diff can fail even if expected == results due to string formatting
    # Only show if the dict comparison fails
    assert not show_diff(str_results, result_json_1.replace('\\\"','\\\\\\"'))


  # Check numbers
  records = geo_container.records_list
  entries = geo_container.entries_list
  assert len(records) == len(entries)
  assert len(entries) ==   1369
  if not geo_container.has_proxies():
    geo_container.build_proxies()
  assert geo_container.has_proxies()
  assert len(geo_container.proxies_list) == len(entries)-len(geo_container.entries["Nonbonded"])

def tst_02(model,printing=False):
  # Test a 1yjp with NO labels and YES a model
  # (i_seqs present in .geo string because no labels, will build proxies)

  expected= json.loads(result_json_2)

  grm = model.restraints_manager.geometry

  buffer = StringIO()
  grm.write_geo_file(model.get_sites_cart(),site_labels=None,file_descriptor=buffer)
  geo_str = buffer.getvalue()
  geo_lines = geo_str.split("\n")
  geo_container = GeoParser(geo_lines,model=model)


  if printing:
    print("\n\ntst_02")
  results, str_results = extract_results(geo_container,print_result=printing)

  # Check values
  if not expected==results:
    assert not show_diff(str_results, result_json_2.replace('\\\"','\\\\\\"'))

  # Check numbers
  records = geo_container.records_list
  entries = geo_container.entries_list
  assert len(records) == len(entries)
  assert len(entries) ==   1369
  if not geo_container.has_proxies():
    geo_container.build_proxies()
  assert geo_container.has_proxies()
  assert len(geo_container.proxies_list) == len(entries)-len(geo_container.entries["Nonbonded"])

def tst_03(model,printing=False):
  # Test a 1yjp with NO labels and NO a model
  # (i_seqs present in .geo string because no labels, will build proxies)

  expected= json.loads(result_json_2)

  grm = model.restraints_manager.geometry

  buffer = StringIO()
  grm.write_geo_file(model.get_sites_cart(),site_labels=None,file_descriptor=buffer)
  geo_str = buffer.getvalue()
  geo_lines = geo_str.split("\n")
  geo_container = GeoParser(geo_lines,model=model)



  if printing:
    print("\n\ntst_03")
  results, str_results = extract_results(geo_container,print_result=printing)

  # Check values
  if not expected==results:
    assert not show_diff(str_results, result_json_2.replace('\\\"','\\\\\\"'))

  # Check numbers
  records = geo_container.records_list
  entries = geo_container.entries_list
  assert len(records) == len(entries)
  assert len(entries) ==   1369
  if not geo_container.has_proxies():
    geo_container.build_proxies()
  assert geo_container.has_proxies()
  assert len(geo_container.proxies_list) == len(entries)-len(geo_container.entries["Nonbonded"])

def tst_04(model,printing=False):
  # Test a 1yjp with YES labels and NO a model
  # (i_seqs not present in .geo string and not moel, cannot build proxies)

  expected= json.loads(result_json_3)

  grm = model.restraints_manager.geometry

  buffer = StringIO()
  site_labels = model.get_xray_structure().scatterers().extract_labels()
  grm.write_geo_file(model.get_sites_cart(),site_labels=site_labels,file_descriptor=buffer)
  geo_str = buffer.getvalue()
  geo_lines = geo_str.split("\n")
  geo_container = GeoParser(geo_lines,model=None)



  if printing:
    print("\n\ntst_04")
  results, str_results = extract_results(geo_container,print_result=printing)

  # Check values
  if not expected==results:
    assert not show_diff(str_results, result_json_3.replace('\\\"','\\\\\\"'))

  # Check numbers
  records = geo_container.records_list
  entries = geo_container.entries_list
  assert len(records) == len(entries)
  assert len(entries) ==   1369
  assert not geo_container.has_proxies()

def tst_05(model,printing=False):
  # Test reading complicated geo file
  # YES labels and NO model, cannot build proxies

  expected= json.loads(result_json_4)

  geo_lines = tst_1_geo.split("\n")
  geo_container = GeoParser(geo_lines,model=None)



  if printing:
    print("\n\ntst_05")
  results, str_results = extract_results(geo_container,print_result=printing)

  # Check values
  if not expected==results:
    assert not show_diff(str_results, result_json_4.replace('\\\"','\\\\\\"'))

  # Check numbers
  records = geo_container.records_list
  entries = geo_container.entries_list
  assert len(records) == len(entries)
  assert len(entries) ==   39
  assert not geo_container.has_proxies()
  assert "sym.op." in records[6] and records[6]["sym.op."].strip() == '-x-1/2,y+1/2,-z+3/4', "Failure to parse symmetry op, "+str(records[6]['sym.op.'].strip())

  origins = [entry.origin_label for entry in geo_container.entries_list]
  expected_origins = ['covalent geometry', 'Misc.', 'link_ALPHA2-6', 'Disulphide bridge', 'link_NAG-ASN', '-like', 'Metal coordination', 'link_TRANS', 'Custom Glycosidic', 'link_BETA1-6', 'link_BETA1-3', 'User supplied', 'link_BETA1-4', 'covalent geometry', 'link_ALPHA2-6', 'Disulphide bridge', 'link_NAG-ASN', 'Secondary Structure restraints around h-bond', 'link_TRANS', 'Custom Glycosidic', 'link_BETA1-6', 'link_BETA1-3', 'link_BETA1-4', 'covalent geometry', 'C-Beta improper', 'Side chain', 'link_TRANS', 'covalent geometry', 'link_ALPHA2-6', 'link_NAG-ASN', 'link_BETA1-6', 'link_BETA1-3', 'link_BETA1-4', 'covalent geometry', 'link_NAG-ASN', 'link_TRANS', '', 'Stacking parallelity', 'Basepair parallelity']
  assert origins ==expected_origins, "Got:" + str(origins)

def tst_06(model,printing=False):
  # Test reading complicated geo file
  # Use 'dummy' i_seqs and 1yjp to simulate a small model with complex a .geo file
  # NO labels (so i_seqs) and NO model, will build proxies

  expected= json.loads(result_json_5)


  tst_1_geo_iseqs = replace_idstr_with_int(tst_1_geo,max_int=30)
  geo_lines = tst_1_geo_iseqs.split("\n")
  geo_container = GeoParser(geo_lines,model=model)

  if printing:
    print("\n\ntst_06")
  results, str_results = extract_results(geo_container,print_result=printing)

  # Check values
  if not expected==results:
    assert not show_diff(str_results, result_json_5.replace('\\\"','\\\\\\"'))

  # Check numbers
  records = geo_container.records_list
  entries = geo_container.entries_list
  assert len(records) == len(entries)
  assert len(entries) ==   39
  if not geo_container.has_proxies():
    geo_container.build_proxies()
  assert geo_container.has_proxies()
  assert len(geo_container.proxies_list) == len(entries)-len(geo_container.entries["Nonbonded"])
  # Check for successful build of the single asu text proxy
  proxy_class_list = [type(p) for p in geo_container.proxies_list]
  assert geometry_restraints.bond_asu_proxy in proxy_class_list

  origins = [entry.origin_label for entry in geo_container.entries_list]
  expected_origins = ['covalent geometry', 'Misc.', 'link_ALPHA2-6', 'Disulphide bridge', 'link_NAG-ASN', '-like', 'Metal coordination', 'link_TRANS', 'Custom Glycosidic', 'link_BETA1-6', 'link_BETA1-3', 'User supplied', 'link_BETA1-4', 'covalent geometry', 'link_ALPHA2-6', 'Disulphide bridge', 'link_NAG-ASN', 'Secondary Structure restraints around h-bond', 'link_TRANS', 'Custom Glycosidic', 'link_BETA1-6', 'link_BETA1-3', 'link_BETA1-4', 'covalent geometry', 'C-Beta improper', 'Side chain', 'link_TRANS', 'covalent geometry', 'link_ALPHA2-6', 'link_NAG-ASN', 'link_BETA1-6', 'link_BETA1-3', 'link_BETA1-4', 'covalent geometry', 'link_NAG-ASN', 'link_TRANS', '', 'Stacking parallelity', 'Basepair parallelity']
  assert origins == expected_origins, "Got: "+ str(origins)


def tst_07(model,printing=False):
  # Test reading integer (i_seq) geo file.

  result_js = """
{
  "Bond": [
    {
      "i_seqs": [
        0,
        1
      ],
      "atom_labels": [
      ],
      "ideal": 1.451,
      "model": 1.507,
      "delta": -0.056,
      "sigma": 0.016,
      "weight": 3910.0,
      "residual": 12.3,
      "origin_id": 0
    }
  ]
}
  """
  expected = json.loads(result_js)

  geo_lines = tst_2_geo.split("\n")
  geo_container = GeoParser(geo_lines,model=model)


  if printing:
    print("\n\ntst_07")
  results, str_results = extract_results(geo_container,print_result=printing)


  # Check values
  if not expected==results:
    assert not show_diff(str_results, result_js.replace('\\\"','\\\\\\"'))


  # Check numbers
  records = geo_container.records_list
  entries = geo_container.entries_list
  assert len(records) == len(entries)
  assert len(entries) ==   4
  if not geo_container.has_proxies():
    geo_container.build_proxies()
  assert geo_container.has_proxies()
  assert len(geo_container.proxies_list) == len(entries)-len(geo_container.entries["Nonbonded"])



def main():
  printing = False # Print results
  model = init_model()
  tst_01(model,printing=printing)
  tst_02(model,printing=printing)
  tst_03(model,printing=printing)
  tst_04(model,printing=printing)
  tst_05(model,printing=printing)
  tst_06(model,printing=printing)
  tst_07(model,printing=printing)
  print('OK')

if __name__ == '__main__':
  main()
