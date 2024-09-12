import json
import re
from io import StringIO
from itertools import chain
from iotbx.data_manager import DataManager
from cctbx.crystal.tst_super_cell import pdb_str_1yjp
from mmtbx.geometry_restraints.geo_file_parsing import GeoParseContainer

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
geo_container = GeoParseContainer(geo_lines,model=model)
  

# access entries
entry = geo_container.entries["dihedral"][0]

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
 'harmonic': '0',
 'sigma': 5.0,
 'weight': 0.04,
 'residual': 7.6,
 'origin_id': 0}


"""
tst_1_geo = """
# Geometry restraints

Bond restraints: 637113
Sorted by residual:
bond pdb=" C   KBEDW   1 "
     pdb=" N   DPPDW   2 "
  ideal  model  delta    sigma   weight residual
  1.329  1.498 -0.169 1.40e-02 5.10e+03 1.45e+02
bond pdb=" C   KBEHW   1 "
     pdb=" N   DPPHW   2 "
  ideal  model  delta    sigma   weight residual
  1.329  1.493 -0.164 1.40e-02 5.10e+03 1.37e+02
bond pdb=" C   DPPBW   2 "
     pdb=" N   SERBW   3 "
  ideal  model  delta    sigma   weight residual
  1.329  1.492 -0.163 1.40e-02 5.10e+03 1.35e+02

Metal coordination restraints: 16
Sorted by residual:
bond pdb=" SG  CYSA4  11 "
     pdb="ZN    ZNA4 102 "
  ideal  model  delta    sigma   weight residual
  2.318  1.877  0.441 2.70e-02 1.37e+03 2.66e+02
bond pdb=" SG  CYSE4  14 "
     pdb="ZN    ZNE4 101 "
  ideal  model  delta    sigma   weight residual
  2.318  1.878  0.440 2.70e-02 1.37e+03 2.66e+02


Misc. restraints: 4
Sorted by residual:
bond pdb=" NG  DPPBW   2 "
     pdb=" C   5OHBW   6 "
  ideal  model  delta    sigma   weight residual
  1.430  1.512 -0.082 1.00e-02 1.00e+04 6.67e+01
bond pdb=" NG  DPPHW   2 "
     pdb=" C   5OHHW   6 "
  ideal  model  delta    sigma   weight residual
  1.430  1.502 -0.072 1.00e-02 1.00e+04 5.23e+01


Bond angle restraints: 949535
Sorted by residual:
angle pdb=" N   GLNC5 122 "
      pdb=" CA  GLNC5 122 "
      pdb=" C   GLNC5 122 "
    ideal   model   delta    sigma   weight residual
   108.78  121.25  -12.47 8.20e-01 1.49e+00 2.31e+02
angle pdb=" C   THRAD 151 "
      pdb=" N   PROAD 152 "
      pdb=" CA  PROAD 152 "
    ideal   model   delta    sigma   weight residual
   119.19  103.99   15.20 1.06e+00 8.90e-01 2.06e+02

Secondary Structure restraints around h-bond angle restraints: 36
Sorted by residual:
angle pdb=" C   PHE A  13 "
      pdb=" O   PHE A  13 "
      pdb=" N   TYR A  17 "
    ideal   model   delta    sigma   weight residual
   155.00  172.15  -17.15 1.00e+01 1.00e-02 2.94e+00
angle pdb=" C   SER G  14 "
      pdb=" O   SER G  14 "
      pdb=" N   CYS G  18 "
    ideal   model   delta    sigma   weight residual
   155.00  163.04   -8.04 5.00e+00 4.00e-02 2.58e+00
angle pdb=" C   PHE E  13 "
      pdb=" O   PHE E  13 "
      pdb=" N   TYR E  17 "
    ideal   model   delta    sigma   weight residual
   155.00  169.75  -14.75 1.00e+01 1.00e-02 2.17e+00

Dihedral angle restraints: 360110
  sinusoidal: 334059
    harmonic: 26051
Sorted by residual:
dihedral pdb=" CA  TRPBV 218 "
         pdb=" C   TRPBV 218 "
         pdb=" N   HISBV 219 "
         pdb=" CA  HISBV 219 "
    ideal   model   delta  harmonic     sigma   weight residual
  -180.00 -133.25  -46.75     0      5.00e+00 4.00e-02 8.74e+01
dihedral pdb=" O4'   UCA1174 "
         pdb=" C1'   UCA1174 "
         pdb=" N1    UCA1174 "
         pdb=" C2    UCA1174 "
    ideal   model   delta sinusoidal    sigma   weight residual
   200.00   22.28  177.72     1      1.50e+01 4.44e-03 8.53e+01
dihedral pdb=" O4'   UAA 546 "
         pdb=" C1'   UAA 546 "
         pdb=" N1    UAA 546 "
         pdb=" C2    UAA 546 "
    ideal   model   delta sinusoidal    sigma   weight residual
   200.00   28.25  171.75     1      1.50e+01 4.44e-03 8.49e+01


C-Beta improper torsion angle restraints: 47272
Sorted by residual:
dihedral pdb=" C   UALFW   5 "
         pdb=" N   UALFW   5 "
         pdb=" CA  UALFW   5 "
         pdb=" CB  UALFW   5 "
    ideal   model   delta  harmonic     sigma   weight residual
  -122.60 -179.97   57.37     0      2.50e+00 1.60e-01 5.27e+02
dihedral pdb=" N   UALFW   5 "
         pdb=" C   UALFW   5 "
         pdb=" CA  UALFW   5 "
         pdb=" CB  UALFW   5 "
    ideal   model   delta  harmonic     sigma   weight residual
   122.80  179.97  -57.17     0      2.50e+00 1.60e-01 5.23e+02

Chirality restraints: 120781
Sorted by residual:
chirality pdb=" CB  VALE5 108 "
          pdb=" CA  VALE5 108 "
          pdb=" CG1 VALE5 108 "
          pdb=" CG2 VALE5 108 "
  both_signs  ideal   model   delta    sigma   weight residual
    False     -2.63    0.29   -2.92 2.00e-01 2.50e+01 2.13e+02
chirality pdb=" P     AGA1866 "
          pdb=" OP1   AGA1866 "
          pdb=" OP2   AGA1866 "
          pdb=" O5'   AGA1866 "
  both_signs  ideal   model   delta    sigma   weight residual
    True       2.41    0.24    2.17 2.00e-01 2.50e+01 1.18e+02

Planarity restraints: 53425
Sorted by residual:
                               delta    sigma   weight rms_deltas residual
plane pdb=" N   UALFW   5 "   -0.055 2.00e-02 2.50e+03   1.19e-01 1.76e+02
      pdb=" CA  UALFW   5 "   -0.042 2.00e-02 2.50e+03
      pdb=" C   UALFW   5 "    0.115 2.00e-02 2.50e+03
      pdb=" CB  UALFW   5 "   -0.171 2.00e-02 2.50e+03
      pdb=" N1  UALFW   5 "    0.153 2.00e-02 2.50e+03
                               delta    sigma   weight rms_deltas residual
plane pdb=" N   UALDW   5 "   -0.024 2.00e-02 2.50e+03   6.39e-02 5.11e+01
      pdb=" CA  UALDW   5 "   -0.027 2.00e-02 2.50e+03
      pdb=" C   UALDW   5 "    0.063 2.00e-02 2.50e+03
      pdb=" CB  UALDW   5 "   -0.093 2.00e-02 2.50e+03
      pdb=" N1  UALDW   5 "    0.080 2.00e-02 2.50e+03
                               delta    sigma   weight rms_deltas residual
plane pdb=" C1'   AGA1095 "    0.035 2.00e-02 2.50e+03   2.83e-02 2.20e+01
      pdb=" N9    AGA1095 "   -0.078 2.00e-02 2.50e+03
      pdb=" C8    AGA1095 "    0.003 2.00e-02 2.50e+03
      pdb=" N7    AGA1095 "    0.001 2.00e-02 2.50e+03
      pdb=" C5    AGA1095 "    0.022 2.00e-02 2.50e+03
      pdb=" C6    AGA1095 "    0.015 2.00e-02 2.50e+03
      pdb=" N6    AGA1095 "   -0.008 2.00e-02 2.50e+03
      pdb=" N1    AGA1095 "   -0.009 2.00e-02 2.50e+03
      pdb=" C2    AGA1095 "   -0.004 2.00e-02 2.50e+03
      pdb=" N3    AGA1095 "   -0.001 2.00e-02 2.50e+03
      pdb=" C4    AGA1095 "    0.024 2.00e-02 2.50e+03

Stacking parallelity restraints: 2
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

    plane 1                plane 2                residual  delta(deg) sigma
    pdb=" C1'  DG A   9 "  pdb=" C1'  DC A  10 "  6.00e+00   5.3613    0.0270
    pdb=" N9   DG A   9 "  pdb=" N1   DC A  10 "
    pdb=" C8   DG A   9 "  pdb=" C2   DC A  10 "
    pdb=" N7   DG A   9 "  pdb=" O2   DC A  10 "
    pdb=" C5   DG A   9 "  pdb=" N3   DC A  10 "
    pdb=" C6   DG A   9 "  pdb=" C4   DC A  10 "
    pdb=" O6   DG A   9 "  pdb=" N4   DC A  10 "
    pdb=" N1   DG A   9 "  pdb=" C5   DC A  10 "
    pdb=" C2   DG A   9 "  pdb=" C6   DC A  10 "
    pdb=" N2   DG A   9 "  
    pdb=" N3   DG A   9 "  
    pdb=" C4   DG A   9 "  


Basepair parallelity restraints: 2
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

    plane 1                plane 2                residual  delta(deg) sigma
    pdb=" C1'  DC A  10 "  pdb=" C1'  DG B  11 "  5.40e+00   6.3101    0.0335
    pdb=" N1   DC A  10 "  pdb=" N9   DG B  11 "
    pdb=" C2   DC A  10 "  pdb=" C8   DG B  11 "
    pdb=" O2   DC A  10 "  pdb=" N7   DG B  11 "
    pdb=" N3   DC A  10 "  pdb=" C5   DG B  11 "
    pdb=" C4   DC A  10 "  pdb=" C6   DG B  11 "
    pdb=" N4   DC A  10 "  pdb=" O6   DG B  11 "
    pdb=" C5   DC A  10 "  pdb=" N1   DG B  11 "
    pdb=" C6   DC A  10 "  pdb=" C2   DG B  11 "
                           pdb=" N2   DG B  11 "
                           pdb=" N3   DG B  11 "
                           pdb=" C4   DG B  11 "



Nonbonded interactions: 5594317
Sorted by model distance:
nonbonded pdb="MG    MGAA3098 "
          pdb=" O   HOHAA3655 "
   model   vdw
   1.653 2.170
nonbonded pdb="MG    MGEA3013 "
          pdb=" O   HOHEA3263 "
   model   vdw
   1.659 2.170
nonbonded pdb="MG    MGAA3029 "
          pdb=" O   HOHAA3329 "
   model   vdw
   1.675 2.170
nonbonded pdb="MG    MGDA1642 "
          pdb=" O   HOHDA1879 "
   model   vdw
   1.699 2.170
"""

tst_2_geo = """
# Geometry restraints

Bond restraints: 59
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


def replace_idstr_with_int(text,max_int=100):
  """
  Replace id_strs in a geo_file str with integers
    For testing
  """
  index = 0
  def replacement(match):
    nonlocal index # allow accessing index
    original_length = len(match.group(0))
    replacement_text = f'{index}'.ljust(original_length)
    index += 1
    if index>max_int:
      index = 0
    return replacement_text

  new_text = re.sub(r'pdb=".*?"', replacement, text)
  return new_text


def extract_results(container,print_result=False):
  results = {}
    
  for entry_name,entries in container.entries.items():
    results[entry_name] = []
    idxs = [0,-1]
    for idx in idxs:
      entry = entries[idx]
      value_idxs = [0,1,2,3] # just take a few values to check
      for i,(key,value) in enumerate(entry.record.items()):
        if i in value_idxs:
          results[entry_name].append(value)

  if print_result:
    #Print output to make tests
    print("expected= {")
    for key,value in results.items():
      print(f"'{key}':",value,",")
    print("}")
  return results

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
  expected= {
  'bond': [[0, 1], ['pdb=" N   GLY A   1 "', 'pdb=" CA  GLY A   1 "'], 1.451, 1.507, [6, 12], ['pdb=" C   ASN A   2 "', 'pdb=" N   ASN A   3 "'], 1.331, 1.33] ,
  'angle': [[12, 13, 14], ['pdb=" N   ASN A   3 "', 'pdb=" CA  ASN A   3 "', 'pdb=" C   ASN A   3 "'], 108.9, 113.48, [47, 48, 49], ['pdb=" CA  TYR A   7 "', 'pdb=" C   TYR A   7 "', 'pdb=" O   TYR A   7 "'], 121.0, 120.98] ,
  'dihedral': [[13, 14, 20, 21], ['pdb=" CA  ASN A   3 "', 'pdb=" C   ASN A   3 "', 'pdb=" N   GLN A   4 "', 'pdb=" CA  GLN A   4 "'], 180.0, 166.21, [14, 12, 13, 16], ['pdb=" C   ASN A   3 "', 'pdb=" N   ASN A   3 "', 'pdb=" CA  ASN A   3 "', 'pdb=" CB  ASN A   3 "'], -122.6, -122.56] ,
  'chiral': [[30, 29, 31, 33], ['pdb=" CA  GLN A   5 "', 'pdb=" N   GLN A   5 "', 'pdb=" C   GLN A   5 "', 'pdb=" CB  GLN A   5 "'], 'False', 2.51, [13, 12, 14, 16], ['pdb=" CA  ASN A   3 "', 'pdb=" N   ASN A   3 "', 'pdb=" C   ASN A   3 "', 'pdb=" CB  ASN A   3 "'], 'False', 2.51] ,
  'plane': [[50, 51, 52, 53, 54, 55, 56, 57], ['pdb=" CB  TYR A   7 "', 'pdb=" CG  TYR A   7 "', 'pdb=" CD1 TYR A   7 "', 'pdb=" CD2 TYR A   7 "', 'pdb=" CE1 TYR A   7 "', 'pdb=" CE2 TYR A   7 "', 'pdb=" CZ  TYR A   7 "', 'pdb=" OH  TYR A   7 "'], [-0.006, 0.022, -0.008, -0.004, 0.002, -0.001, -0.011, 0.006], [0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02], [13, 14, 15], ['pdb=" CA  ASN A   3 "', 'pdb=" C   ASN A   3 "', 'pdb=" O   ASN A   3 "'], ['0.000', -0.002, 0.001], [0.02, 0.02, 0.02]] ,
  'nonbonded': [[57, 62], ['pdb=" OH  TYR A   7 "', 'pdb=" O   HOH A  11 "'], 2.525, 3.04, [38, 51], ['pdb=" N   ASN A   6 "', 'pdb=" CG  TYR A   7 "'], 4.9, 3.34] ,
  }

  grm = model.restraints_manager.geometry
  
  buffer = StringIO()
  site_labels = model.get_xray_structure().scatterers().extract_labels()
  grm.write_geo_file(model.get_sites_cart(),site_labels=site_labels,file_descriptor=buffer)
  geo_str = buffer.getvalue()
  geo_lines = geo_str.split("\n")
  geo_container = GeoParseContainer(geo_lines,model=model)
  
  if printing:
    print("\n\ntst_01")
  results = extract_results(geo_container,print_result=printing)
  
  # Check values
  assert expected==results
  
  # Check numbers
  records = geo_container.records_list
  entries = geo_container.entries_list
  assert len(records) == len(entries)
  assert len(entries) ==   1369
  assert geo_container.has_proxies
  assert len(geo_container.proxies_list) == len(entries)-len(geo_container.entries["nonbonded"])

def tst_02(model,printing=False):
  # Test a 1yjp with NO labels and YES a model
  # (i_seqs present in .geo string because no labels, will build proxies)
    
  expected= {
  'bond': [[0, 1], ['0', '1'], 1.451, 1.507, [6, 12], ['6', '12'], 1.331, 1.33] ,
  'angle': [[12, 13, 14], ['12', '13', '14'], 108.9, 113.48, [47, 48, 49], ['47', '48', '49'], 121.0, 120.98] ,
  'dihedral': [[13, 14, 20, 21], ['13', '14', '20', '21'], 180.0, 166.21, [14, 12, 13, 16], ['14', '12', '13', '16'], -122.6, -122.56] ,
  'chiral': [[30, 29, 31, 33], ['30', '29', '31', '33'], 'False', 2.51, [13, 12, 14, 16], ['13', '12', '14', '16'], 'False', 2.51] ,
  'plane': [[50, 51, 52, 53, 54, 55, 56, 57], ['50', '51', '52', '53', '54', '55', '56', '57'], [-0.006, 0.022, -0.008, -0.004, 0.002, -0.001, -0.011, 0.006], [0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02], [13, 14, 15], ['13', '14', '15'], ['0.000', -0.002, 0.001], [0.02, 0.02, 0.02]] ,
  'nonbonded': [[57, 62], ['57', '62'], 2.525, 3.04, [38, 51], ['38', '51'], 4.9, 3.34] ,
  }

  grm = model.restraints_manager.geometry
  
  buffer = StringIO()
  grm.write_geo_file(model.get_sites_cart(),site_labels=None,file_descriptor=buffer)
  geo_str = buffer.getvalue()
  geo_lines = geo_str.split("\n")
  geo_container = GeoParseContainer(geo_lines,model=model)
  
  
  if printing:
    print("\n\ntst_02")
  results = extract_results(geo_container,print_result=printing)
  
  # Check values
  assert expected==results
  
  # Check numbers
  records = geo_container.records_list
  entries = geo_container.entries_list
  assert len(records) == len(entries)
  assert len(entries) ==   1369
  assert geo_container.has_proxies
  assert len(geo_container.proxies_list) == len(entries)-len(geo_container.entries["nonbonded"])

def tst_03(model,printing=False):
  # Test a 1yjp with NO labels and NO a model
  # (i_seqs present in .geo string because no labels, will build proxies)
  
  expected= {
  'bond': [[0, 1], ['0', '1'], 1.451, 1.507, [6, 12], ['6', '12'], 1.331, 1.33] ,
  'angle': [[12, 13, 14], ['12', '13', '14'], 108.9, 113.48, [47, 48, 49], ['47', '48', '49'], 121.0, 120.98] ,
  'dihedral': [[13, 14, 20, 21], ['13', '14', '20', '21'], 180.0, 166.21, [14, 12, 13, 16], ['14', '12', '13', '16'], -122.6, -122.56] ,
  'chiral': [[30, 29, 31, 33], ['30', '29', '31', '33'], 'False', 2.51, [13, 12, 14, 16], ['13', '12', '14', '16'], 'False', 2.51] ,
  'plane': [[50, 51, 52, 53, 54, 55, 56, 57], ['50', '51', '52', '53', '54', '55', '56', '57'], [-0.006, 0.022, -0.008, -0.004, 0.002, -0.001, -0.011, 0.006], [0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02], [13, 14, 15], ['13', '14', '15'], ['0.000', -0.002, 0.001], [0.02, 0.02, 0.02]] ,
  'nonbonded': [[57, 62], ['57', '62'], 2.525, 3.04, [38, 51], ['38', '51'], 4.9, 3.34] ,
  }


  grm = model.restraints_manager.geometry
  
  buffer = StringIO()
  grm.write_geo_file(model.get_sites_cart(),site_labels=None,file_descriptor=buffer)
  geo_str = buffer.getvalue()
  geo_lines = geo_str.split("\n")
  geo_container = GeoParseContainer(geo_lines,model=model)
  
  
  
  if printing:
    print("\n\ntst_03")
  results = extract_results(geo_container,print_result=printing)
  
  # Check values
  assert expected==results
  
  # Check numbers
  records = geo_container.records_list
  entries = geo_container.entries_list
  assert len(records) == len(entries)
  assert len(entries) ==   1369
  assert geo_container.has_proxies
  assert len(geo_container.proxies_list) == len(entries)-len(geo_container.entries["nonbonded"])

def tst_04(model,printing=False):
  # Test a 1yjp with YES labels and NO a model
  # (i_seqs not present in .geo string and not moel, cannot build proxies)
    
  expected= {
  'bond': [[], ['pdb=" N   GLY A   1 "', 'pdb=" CA  GLY A   1 "'], 1.451, 1.507, [], ['pdb=" C   ASN A   2 "', 'pdb=" N   ASN A   3 "'], 1.331, 1.33] ,
  'angle': [[], ['pdb=" N   ASN A   3 "', 'pdb=" CA  ASN A   3 "', 'pdb=" C   ASN A   3 "'], 108.9, 113.48, [], ['pdb=" CA  TYR A   7 "', 'pdb=" C   TYR A   7 "', 'pdb=" O   TYR A   7 "'], 121.0, 120.98] ,
  'dihedral': [[], ['pdb=" CA  ASN A   3 "', 'pdb=" C   ASN A   3 "', 'pdb=" N   GLN A   4 "', 'pdb=" CA  GLN A   4 "'], 180.0, 166.21, [], ['pdb=" C   ASN A   3 "', 'pdb=" N   ASN A   3 "', 'pdb=" CA  ASN A   3 "', 'pdb=" CB  ASN A   3 "'], -122.6, -122.56] ,
  'chiral': [[], ['pdb=" CA  GLN A   5 "', 'pdb=" N   GLN A   5 "', 'pdb=" C   GLN A   5 "', 'pdb=" CB  GLN A   5 "'], 'False', 2.51, [], ['pdb=" CA  ASN A   3 "', 'pdb=" N   ASN A   3 "', 'pdb=" C   ASN A   3 "', 'pdb=" CB  ASN A   3 "'], 'False', 2.51] ,
  'plane': [[], ['pdb=" CB  TYR A   7 "', 'pdb=" CG  TYR A   7 "', 'pdb=" CD1 TYR A   7 "', 'pdb=" CD2 TYR A   7 "', 'pdb=" CE1 TYR A   7 "', 'pdb=" CE2 TYR A   7 "', 'pdb=" CZ  TYR A   7 "', 'pdb=" OH  TYR A   7 "'], [-0.006, 0.022, -0.008, -0.004, 0.002, -0.001, -0.011, 0.006], [0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02], [], ['pdb=" CA  ASN A   3 "', 'pdb=" C   ASN A   3 "', 'pdb=" O   ASN A   3 "'], ['0.000', -0.002, 0.001], [0.02, 0.02, 0.02]] ,
  'nonbonded': [[], ['pdb=" OH  TYR A   7 "', 'pdb=" O   HOH A  11 "'], 2.525, 3.04, [], ['pdb=" N   ASN A   6 "', 'pdb=" CG  TYR A   7 "'], 4.9, 3.34] ,
  }

  grm = model.restraints_manager.geometry
  
  buffer = StringIO()
  site_labels = model.get_xray_structure().scatterers().extract_labels()
  grm.write_geo_file(model.get_sites_cart(),site_labels=site_labels,file_descriptor=buffer)
  geo_str = buffer.getvalue()
  geo_lines = geo_str.split("\n")
  geo_container = GeoParseContainer(geo_lines,model=None)
  
  
  
  if printing:
    print("\n\ntst_04")
  results = extract_results(geo_container,print_result=printing)
  
  # Check values
  assert expected==results
  
  # Check numbers
  records = geo_container.records_list
  entries = geo_container.entries_list
  assert len(records) == len(entries)
  assert len(entries) ==   1369
  assert not geo_container.has_proxies

def tst_05(model,printing=False):
  # Test reading complicated geo file
  # YES labels and NO model, cannot build proxies
  
  expected= {
  'bond': [[], ['pdb=" C   KBEDW   1 "', 'pdb=" N   DPPDW   2 "'], 1.329, 1.498, [], ['pdb=" NG  DPPHW   2 "', 'pdb=" C   5OHHW   6 "'], 1.43, 1.502] ,
  'angle': [[], ['pdb=" N   GLNC5 122 "', 'pdb=" CA  GLNC5 122 "', 'pdb=" C   GLNC5 122 "'], 108.78, 121.25, [], ['pdb=" C   PHE E  13 "', 'pdb=" O   PHE E  13 "', 'pdb=" N   TYR E  17 "'], 155.0, 169.75] ,
  'dihedral': [[], ['pdb=" CA  TRPBV 218 "', 'pdb=" C   TRPBV 218 "', 'pdb=" N   HISBV 219 "', 'pdb=" CA  HISBV 219 "'], -180.0, -133.25, [], ['pdb=" N   UALFW   5 "', 'pdb=" C   UALFW   5 "', 'pdb=" CA  UALFW   5 "', 'pdb=" CB  UALFW   5 "'], 122.8, 179.97] ,
  'chiral': [[], ['pdb=" CB  VALE5 108 "', 'pdb=" CA  VALE5 108 "', 'pdb=" CG1 VALE5 108 "', 'pdb=" CG2 VALE5 108 "'], 'False', -2.63, [], ['pdb=" P     AGA1866 "', 'pdb=" OP1   AGA1866 "', 'pdb=" OP2   AGA1866 "', 'pdb=" O5\'   AGA1866 "'], 'True', 2.41] ,
  'plane': [[], ['pdb=" N   UALFW   5 "', 'pdb=" CA  UALFW   5 "', 'pdb=" C   UALFW   5 "', 'pdb=" CB  UALFW   5 "', 'pdb=" N1  UALFW   5 "'], [-0.055, -0.042, 0.115, -0.171, 0.153], [0.02, 0.02, 0.02, 0.02, 0.02], [], ['pdb=" C1\'   AGA1095 "', 'pdb=" N9    AGA1095 "', 'pdb=" C8    AGA1095 "', 'pdb=" N7    AGA1095 "', 'pdb=" C5    AGA1095 "', 'pdb=" C6    AGA1095 "', 'pdb=" N6    AGA1095 "', 'pdb=" N1    AGA1095 "', 'pdb=" C2    AGA1095 "', 'pdb=" N3    AGA1095 "'], [0.035, -0.078, 0.003, 0.001, 0.022, 0.015, -0.008, -0.009, -0.004, -0.001], [0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02]] ,
  'parallelity': [[], [], ['pdb=" C1\'  DG B  11 "', 'pdb=" N9   DG B  11 "', 'pdb=" C8   DG B  11 "', 'pdb=" N7   DG B  11 "', 'pdb=" C5   DG B  11 "', 'pdb=" C6   DG B  11 "', 'pdb=" O6   DG B  11 "', 'pdb=" N1   DG B  11 "', 'pdb=" C2   DG B  11 "', 'pdb=" N2   DG B  11 "', 'pdb=" N3   DG B  11 "', 'pdb=" C4   DG B  11 "'], ['pdb=" C1\'  DC B  12 "', 'pdb=" N1   DC B  12 "', 'pdb=" C2   DC B  12 "', 'pdb=" O2   DC B  12 "', 'pdb=" N3   DC B  12 "', 'pdb=" C4   DC B  12 "', 'pdb=" N4   DC B  12 "', 'pdb=" C5   DC B  12 "', 'pdb=" C6   DC B  12 "'], [], [], ['pdb=" C1\'  DC A  10 "', 'pdb=" N1   DC A  10 "', 'pdb=" C2   DC A  10 "', 'pdb=" O2   DC A  10 "', 'pdb=" N3   DC A  10 "', 'pdb=" C4   DC A  10 "', 'pdb=" N4   DC A  10 "', 'pdb=" C5   DC A  10 "', 'pdb=" C6   DC A  10 "', 'pdb=" N2   DG B  11 "', 'pdb=" N3   DG B  11 "', 'pdb=" C4   DG B  11 "'], ['pdb=" C1\'  DG B  11 "', 'pdb=" N9   DG B  11 "', 'pdb=" C8   DG B  11 "', 'pdb=" N7   DG B  11 "', 'pdb=" C5   DG B  11 "', 'pdb=" C6   DG B  11 "', 'pdb=" O6   DG B  11 "', 'pdb=" N1   DG B  11 "', 'pdb=" C2   DG B  11 "']] ,
  'nonbonded': [[], ['pdb="MG    MGAA3098 "', 'pdb=" O   HOHAA3655 "'], 1.653, 2.17, [], ['pdb="MG    MGDA1642 "', 'pdb=" O   HOHDA1879 "'], 1.699, 2.17] ,
  }


  
  geo_lines = tst_1_geo.split("\n")
  geo_container = GeoParseContainer(geo_lines,model=None)
  
  
  
  if printing:
    print("\n\ntst_05")
  results = extract_results(geo_container,print_result=printing)
  
  
  # Check values
  assert expected==results
  
  # Check numbers
  records = geo_container.records_list
  entries = geo_container.entries_list
  assert len(records) == len(entries)
  assert len(entries) ==   30
  assert not geo_container.has_proxies
  
  origin_ids = [entry.origin_id for entry in geo_container.entries_list]
  assert origin_ids == [0, 0, 0, 3, 3, 9, 9, 0, 0, 2, 2, 2, 0, 0, 0, 11, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

def tst_06(model,printing=False):
  # Test reading complicated geo file
  # Use 'dummy' i_seqs and 1yjp to simulate a small model with complex a .geo file
  # NO labels (so i_seqs) and NO model, will build proxies
  
  expected= {
  'bond': [[0, 1], ['0', '1'], 1.329, 1.498, [12, 13], ['12', '13'], 1.43, 1.502] ,
  'angle': [[14, 15, 16], ['14', '15', '16'], 108.78, 121.25, [26, 27, 28], ['26', '27', '28'], 155.0, 169.75] ,
  'dihedral': [[29, 30, 0, 1], ['29', '30', '0', '1'], -180.0, -133.25, [14, 15, 16, 17], ['14', '15', '16', '17'], 122.8, 179.97] ,
  'chiral': [[18, 19, 20, 21], ['18', '19', '20', '21'], 'False', -2.63, [22, 23, 24, 25], ['22', '23', '24', '25'], 'True', 2.41] ,
  'plane': [[26, 27, 28, 29, 30], ['26', '27', '28', '29', '30'], [-0.055, -0.042, 0.115, -0.171, 0.153], [0.02, 0.02, 0.02, 0.02, 0.02], [5, 6, 7, 8, 9, 10, 11, 12, 13, 14], ['5', '6', '7', '8', '9', '10', '11', '12', '13', '14'], [0.035, -0.078, 0.003, 0.001, 0.022, 0.015, -0.008, -0.009, -0.004, -0.001], [0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02]] ,
  'parallelity': [[16, 18, 20, 22, 24, 26, 28, 30, 1, 3, 4, 5], [17, 19, 21, 23, 25, 27, 29, 0, 2], ['16', '18', '20', '22', '24', '26', '28', '30', '1', '3', '4', '5'], ['17', '19', '21', '23', '25', '27', '29', '0', '2'], [17, 19, 21, 23, 25, 27, 29, 0, 2, 4, 5, 6], [18, 20, 22, 24, 26, 28, 30, 1, 3], ['17', '19', '21', '23', '25', '27', '29', '0', '2', '4', '5', '6'], ['18', '20', '22', '24', '26', '28', '30', '1', '3']] ,
  'nonbonded': [[7, 8], ['7', '8'], 1.653, 2.17, [13, 14], ['13', '14'], 1.699, 2.17] ,
  }


  
  tst_1_geo_iseqs = replace_idstr_with_int(tst_1_geo,max_int=30)
  geo_lines = tst_1_geo_iseqs.split("\n")
  geo_container = GeoParseContainer(geo_lines,model=model)
  
  if printing:
    print("\n\ntst_06")
  results = extract_results(geo_container,print_result=printing)
  
  # Check values
  assert expected==results
  
  # Check numbers
  records = geo_container.records_list
  entries = geo_container.entries_list
  assert len(records) == len(entries)
  assert len(entries) ==   30
  assert geo_container.has_proxies
  assert len(geo_container.proxies_list) == len(entries)-len(geo_container.entries["nonbonded"])
  
  origin_ids = [entry.origin_id for entry in geo_container.entries_list]
  assert origin_ids == [0, 0, 0, 3, 3, 9, 9, 0, 0, 2, 2, 2, 0, 0, 0, 11, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]


def tst_07(model,printing=False):
  # Test reading integer (i_seq) geo file
  
  expected= {
  'bond': [[0, 1], ['0', '1'], 1.451, 1.507, [5, 6], ['5', '6'], 1.524, 1.498] ,
  }

  
  geo_lines = tst_2_geo.split("\n")
  geo_container = GeoParseContainer(geo_lines,model=model)
  
  
  if printing:
    print("\n\ntst_01")
  results = extract_results(geo_container,print_result=printing)
  
  
  # Check values
  assert expected==results
  
  # Check numbers
  records = geo_container.records_list
  entries = geo_container.entries_list
  assert len(records) == len(entries)
  assert len(entries) ==   4
  assert geo_container.has_proxies
  assert len(geo_container.proxies_list) == len(entries)-len(geo_container.entries["nonbonded"])


def tst_08(model,printing=False):
  # Test reading without whitespace between sections


  expected= {
  'bond': [[], ['pdb=" C   KBEDW   1 "', 'pdb=" N   DPPDW   2 "'], 1.329, 1.498, [], ['pdb=" NG  DPPHW   2 "', 'pdb=" C   5OHHW   6 "'], 1.43, 1.502] ,
  'angle': [[], ['pdb=" N   GLNC5 122 "', 'pdb=" CA  GLNC5 122 "', 'pdb=" C   GLNC5 122 "'], 108.78, 121.25, [], ['pdb=" C   PHE E  13 "', 'pdb=" O   PHE E  13 "', 'pdb=" N   TYR E  17 "'], 155.0, 169.75] ,
  'dihedral': [[], ['pdb=" CA  TRPBV 218 "', 'pdb=" C   TRPBV 218 "', 'pdb=" N   HISBV 219 "', 'pdb=" CA  HISBV 219 "'], -180.0, -133.25, [], ['pdb=" N   UALFW   5 "', 'pdb=" C   UALFW   5 "', 'pdb=" CA  UALFW   5 "', 'pdb=" CB  UALFW   5 "'], 122.8, 179.97] ,
  'chiral': [[], ['pdb=" CB  VALE5 108 "', 'pdb=" CA  VALE5 108 "', 'pdb=" CG1 VALE5 108 "', 'pdb=" CG2 VALE5 108 "'], 'False', -2.63, [], ['pdb=" P     AGA1866 "', 'pdb=" OP1   AGA1866 "', 'pdb=" OP2   AGA1866 "', 'pdb=" O5\'   AGA1866 "'], 'True', 2.41] ,
  'plane': [[], ['pdb=" N   UALFW   5 "', 'pdb=" CA  UALFW   5 "', 'pdb=" C   UALFW   5 "', 'pdb=" CB  UALFW   5 "', 'pdb=" N1  UALFW   5 "'], [-0.055, -0.042, 0.115, -0.171, 0.153], [0.02, 0.02, 0.02, 0.02, 0.02], [], ['pdb=" C1\'   AGA1095 "', 'pdb=" N9    AGA1095 "', 'pdb=" C8    AGA1095 "', 'pdb=" N7    AGA1095 "', 'pdb=" C5    AGA1095 "', 'pdb=" C6    AGA1095 "', 'pdb=" N6    AGA1095 "', 'pdb=" N1    AGA1095 "', 'pdb=" C2    AGA1095 "', 'pdb=" N3    AGA1095 "'], [0.035, -0.078, 0.003, 0.001, 0.022, 0.015, -0.008, -0.009, -0.004, -0.001], [0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02]] ,
  'parallelity': [[], [], ['pdb=" C1\'  DG B  11 "', 'pdb=" N9   DG B  11 "', 'pdb=" C8   DG B  11 "', 'pdb=" N7   DG B  11 "', 'pdb=" C5   DG B  11 "', 'pdb=" C6   DG B  11 "', 'pdb=" O6   DG B  11 "', 'pdb=" N1   DG B  11 "', 'pdb=" C2   DG B  11 "', 'pdb=" N2   DG B  11 "', 'pdb=" N3   DG B  11 "', 'pdb=" C4   DG B  11 "'], ['pdb=" C1\'  DC B  12 "', 'pdb=" N1   DC B  12 "', 'pdb=" C2   DC B  12 "', 'pdb=" O2   DC B  12 "', 'pdb=" N3   DC B  12 "', 'pdb=" C4   DC B  12 "', 'pdb=" N4   DC B  12 "', 'pdb=" C5   DC B  12 "', 'pdb=" C6   DC B  12 "'], [], [], ['pdb=" C1\'  DC A  10 "', 'pdb=" N1   DC A  10 "', 'pdb=" C2   DC A  10 "', 'pdb=" O2   DC A  10 "', 'pdb=" N3   DC A  10 "', 'pdb=" C4   DC A  10 "', 'pdb=" N4   DC A  10 "', 'pdb=" C5   DC A  10 "', 'pdb=" C6   DC A  10 "', 'pdb=" N2   DG B  11 "', 'pdb=" N3   DG B  11 "', 'pdb=" C4   DG B  11 "'], ['pdb=" C1\'  DG B  11 "', 'pdb=" N9   DG B  11 "', 'pdb=" C8   DG B  11 "', 'pdb=" N7   DG B  11 "', 'pdb=" C5   DG B  11 "', 'pdb=" C6   DG B  11 "', 'pdb=" O6   DG B  11 "', 'pdb=" N1   DG B  11 "', 'pdb=" C2   DG B  11 "']] ,
  'nonbonded': [[], ['pdb="MG    MGAA3098 "', 'pdb=" O   HOHAA3655 "'], 1.653, 2.17, [], ['pdb="MG    MGDA1642 "', 'pdb=" O   HOHDA1879 "'], 1.699, 2.17] ,
  }



  geo_lines = tst_1_geo.split("\n")
  geo_lines = [line for line in geo_lines if len(line.strip().replace("\n",""))>0]
  geo_container = GeoParseContainer(geo_lines,model=model)

  
  
  if printing:
    print("\n\ntst_08")
  results = extract_results(geo_container,print_result=printing)
  

  assert expected==results
  
  # Check numbers
  records = geo_container.records_list
  entries = geo_container.entries_list
  assert len(records) == len(entries)
  assert len(entries) ==   30
  assert geo_container.has_proxies
  assert len(geo_container.proxies_list) == len(entries)-len(geo_container.entries["nonbonded"])
  
  origin_ids = [entry.origin_id for entry in geo_container.entries_list]
  assert origin_ids == [0, 0, 0, 3, 3, 9, 9, 0, 0, 2, 2, 2, 0, 0, 0, 11, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]



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
  tst_08(model,printing=printing)
  print('OK')

if __name__ == '__main__':
  main()