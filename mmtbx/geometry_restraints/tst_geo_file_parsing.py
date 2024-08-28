import json
import re
from io import StringIO
from itertools import chain
from iotbx.data_manager import DataManager
from cctbx.crystal.tst_super_cell import pdb_str_1yjp
from geo_file_parsing import GeoParseContainer


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


# some random noise text
$$jasfdlksadf4#23454wr4t4t43t
a43tew\"'d.
da43409q38798790980p98pfbond"
asdfdsaangle"ade44deChirality restraints:

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
      value_idxs = [0,1,3,5,7,len(entry.record)] # just take a few values to check
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
dm = DataManager()
dm.process_model_str("1yjp",pdb_str_1yjp)
model= dm.get_model()
model.process(make_restraints=True)

def tst_01():
  # Test a 1yjp with YES labels and YES a model
  # (Can build proxies from label matching to model i_seqs)
  expected= {
  'bond': [0, 1, ' CA  GLY A   1 ', 1.507, 0.016, 6, 12, ' N   ASN A   3 ', 1.33, 0.0134] ,
  'angle': [12, 13, ' N   ASN A   3 ', ' C   ASN A   3 ', 113.48, 47, 48, ' CA  TYR A   7 ', ' O   TYR A   7 ', 120.98] ,
  'dihedral': [13, 14, 21, ' C   ASN A   3 ', ' CA  GLN A   4 ', 14, 12, 16, ' N   ASN A   3 ', ' CB  ASN A   3 '] ,
  'chiral': [30, 29, 33, ' N   GLN A   5 ', ' CB  GLN A   5 ', 13, 12, 16, ' N   ASN A   3 ', ' CB  ASN A   3 '] ,
  'plane': [50, 51, 53, 55, 57, 13, 14, ' CA  ASN A   3 ', ' O   ASN A   3 ', [0.02, 0.02, 0.02]] ,
  'nonbonded': [57, 62, ' O   HOH A  11 ', 3.04, 0, 38, 51, ' CG  TYR A   7 ', 3.34] ,
  }
  grm = model.restraints_manager.geometry
  
  buffer = StringIO()
  site_labels = model.get_xray_structure().scatterers().extract_labels()
  grm.write_geo_file(model.get_sites_cart(),site_labels=site_labels,file_descriptor=buffer)
  geo_str = buffer.getvalue()
  geo_container = GeoParseContainer.from_geo_str(geo_str,model=model)
  
  
  results = extract_results(geo_container)
  
  # Check values
  assert expected==results
  
  # Check numbers
  records = geo_container.records_list
  entries = geo_container.entries_list
  assert len(records) == len(entries)
  assert len(entries) ==   1369
  assert geo_container.has_proxies
  assert len(geo_container.proxies_list) == len(entries)-len(geo_container.entries["nonbonded"])

def tst_02():
  # Test a 1yjp with NO labels and YES a model
  # (i_seqs present in .geo string because no labels, will build proxies)
  
  expected= {
  'bond': [0, 1, '1', 1.507, 0.016, 6, 12, '12', 1.33, 0.0134] ,
  'angle': [12, 13, '12', '14', 113.48, 47, 48, '47', '49', 120.98] ,
  'dihedral': [13, 14, 21, '14', '21', 14, 12, 16, '12', '16'] ,
  'chiral': [30, 29, 33, '29', '33', 13, 12, 16, '12', '16'] ,
  'plane': [50, 51, 53, 55, 57, 13, 14, '13', '15', [0.02, 0.02, 0.02]] ,
  'nonbonded': [57, 62, '62', 3.04, 0, 38, 51, '51', 3.34] ,
  }
  grm = model.restraints_manager.geometry
  
  buffer = StringIO()
  grm.write_geo_file(model.get_sites_cart(),site_labels=None,file_descriptor=buffer)
  geo_str = buffer.getvalue()
  geo_container = GeoParseContainer.from_geo_str(geo_str,model=model)
  
  
  results = extract_results(geo_container)
  
  # Check values
  assert expected==results
  
  # Check numbers
  records = geo_container.records_list
  entries = geo_container.entries_list
  assert len(records) == len(entries)
  assert len(entries) ==   1369
  assert geo_container.has_proxies
  assert len(geo_container.proxies_list) == len(entries)-len(geo_container.entries["nonbonded"])

def tst_03():
  # Test a 1yjp with NO labels and NO a model
  # (i_seqs present in .geo string because no labels, will build proxies)
  
  expected= {
  'bond': [0, 1, '1', 1.507, 0.016, 6, 12, '12', 1.33, 0.0134] ,
  'angle': [12, 13, '12', '14', 113.48, 47, 48, '47', '49', 120.98] ,
  'dihedral': [13, 14, 21, '14', '21', 14, 12, 16, '12', '16'] ,
  'chiral': [30, 29, 33, '29', '33', 13, 12, 16, '12', '16'] ,
  'plane': [50, 51, 53, 55, 57, 13, 14, '13', '15', [0.02, 0.02, 0.02]] ,
  'nonbonded': [57, 62, '62', 3.04, 0, 38, 51, '51', 3.34] ,
  }
  grm = model.restraints_manager.geometry
  
  buffer = StringIO()
  grm.write_geo_file(model.get_sites_cart(),site_labels=None,file_descriptor=buffer)
  geo_str = buffer.getvalue()
  geo_container = GeoParseContainer.from_geo_str(geo_str,model=None)
  
  
  results = extract_results(geo_container)
  
  # Check values
  assert expected==results
  
  # Check numbers
  records = geo_container.records_list
  entries = geo_container.entries_list
  assert len(records) == len(entries)
  assert len(entries) ==   1369
  assert geo_container.has_proxies
  assert len(geo_container.proxies_list) == len(entries)-len(geo_container.entries["nonbonded"])

def tst_04():
  # Test a 1yjp with YES labels and NO a model
  # (i_seqs not present in .geo string and not moel, cannot build proxies)
  
  expected= {
  'bond': [' N   GLY A   1 ', ' CA  GLY A   1 ', 1.507, 0.016, 12.3, ' C   ASN A   2 ', ' N   ASN A   3 ', 1.33, 0.0134, 1.29e-05] ,
  'angle': [' N   ASN A   3 ', ' CA  ASN A   3 ', 108.9, -4.58, 0.376, ' CA  TYR A   7 ', ' C   TYR A   7 ', 121.0, 0.02, 0.111] ,
  'dihedral': [' CA  ASN A   3 ', ' C   ASN A   3 ', ' CA  GLN A   4 ', 166.21, '0', ' C   ASN A   3 ', ' N   ASN A   3 ', ' CB  ASN A   3 ', -122.56, '0'] ,
  'chiral': [' CA  GLN A   5 ', ' N   GLN A   5 ', ' CB  GLN A   5 ', 2.51, 0.12, ' CA  ASN A   3 ', ' N   ASN A   3 ', ' CB  ASN A   3 ', 2.51, '-0.00'] ,
  'plane': [' CB  TYR A   7 ', ' CG  TYR A   7 ', ' CD2 TYR A   7 ', ' CE2 TYR A   7 ', ' OH  TYR A   7 ', ' CA  ASN A   3 ', ' C   ASN A   3 ', ['0.000', -0.002, 0.001], [2500.0, 2500.0, 2500.0], [0.00997, 0.00997, 0.00997]] ,
  'nonbonded': [' OH  TYR A   7 ', ' O   HOH A  11 ', 3.04, 0, ' N   ASN A   6 ', ' CG  TYR A   7 ', 3.34] ,
  }
  grm = model.restraints_manager.geometry
  
  buffer = StringIO()
  site_labels = model.get_xray_structure().scatterers().extract_labels()
  grm.write_geo_file(model.get_sites_cart(),site_labels=site_labels,file_descriptor=buffer)
  geo_str = buffer.getvalue()
  geo_container = GeoParseContainer.from_geo_str(geo_str,model=None)
  
  
  results = extract_results(geo_container)
  
  # Check values
  assert expected==results
  
  # Check numbers
  records = geo_container.records_list
  entries = geo_container.entries_list
  assert len(records) == len(entries)
  assert len(entries) ==   1369
  assert not geo_container.has_proxies

def tst_05():
  # Test reading complicated geo file
  # YES labels and NO model, cannot build proxies
  
  expected= {
  'bond': [' C   KBEDW   1 ', ' N   DPPDW   2 ', 1.498, 0.014, 145.0, ' NG  DPPHW   2 ', ' C   5OHHW   6 ', 1.502, 0.01, 52.3] ,
  'angle': [' N   GLNC5 122 ', ' CA  GLNC5 122 ', 108.78, -12.47, 1.49, ' C   PHE E  13 ', ' O   PHE E  13 ', 155.0, -14.75, 0.01] ,
  'dihedral': [' CA  TRPBV 218 ', ' C   TRPBV 218 ', ' CA  HISBV 219 ', -133.25, '0', ' N   UALFW   5 ', ' C   UALFW   5 ', ' CB  UALFW   5 ', 179.97, '0'] ,
  'chiral': [' CB  VALE5 108 ', ' CA  VALE5 108 ', ' CG2 VALE5 108 ', -2.63, -2.92, ' P     AGA1866 ', ' OP1   AGA1866 ', " O5'   AGA1866 ", 2.41, 2.17] ,
  'plane': [' N   UALFW   5 ', ' CA  UALFW   5 ', ' CB  UALFW   5 ', [-0.055, -0.042, 0.115, -0.171, 0.153], [2500.0, 2500.0, 2500.0, 2500.0, 2500.0], " C1'   AGA1095 ", ' N9    AGA1095 ', ' N7    AGA1095 ', ' C6    AGA1095 ', ' N1    AGA1095 '] ,
  'parallelity': [" C1'  DG B  11 ", ' N9   DG B  11 ', ' N7   DG B  11 ', ' C6   DG B  11 ', ' N1   DG B  11 ', " C1'  DC A  10 ", ' N1   DC A  10 ', ' O2   DC A  10 ', ' C4   DC A  10 ', ' C5   DC A  10 '] ,
  'nonbonded': ['MG    MGAA3098 ', ' O   HOHAA3655 ', 2.17, 'MG    MGDA1642 ', ' O   HOHDA1879 ', 2.17] ,
  }
  
  geo_container = GeoParseContainer.from_geo_str(tst_1_geo)
  
  
  results = extract_results(geo_container)
  
  
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

def tst_06():
  # Test reading complicated geo file
  # Use 'dummy' i_seqs and 1yjp to simulate a small model with complex a .geo file
  # NO labels (so i_seqs) and NO model, will build proxies
  
  expected= {
  'bond': [0, 1, '1', 1.498, 0.014, 12, 13, '13', 1.502, 0.01] ,
  'angle': [14, 15, '14', '16', 121.25, 26, 27, '26', '28', 169.75] ,
  'dihedral': [29, 30, 1, '30', '1', 14, 15, 17, '15', '17'] ,
  'chiral': [18, 19, 21, '19', '21', 22, 23, 25, '23', '25'] ,
  'plane': [26, 27, 29, '26', '28', 5, 6, 8, 10, 12] ,
  'parallelity': [16, 18, 22, 26, 30, 17, 19, 23, 27, 0] ,
  'nonbonded': [7, 8, '8', 2.17, 13, 14, '14', 2.17] ,
  }
  
  tst_1_geo_iseqs = replace_idstr_with_int(tst_1_geo,max_int=30)
  geo_container = GeoParseContainer.from_geo_str(tst_1_geo_iseqs)
  
  results = extract_results(geo_container)
  
  
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



def main():
  tst_01()
  tst_02()
  tst_03()
  tst_04()
  tst_05()
  tst_06()
  print('OK')

if __name__ == '__main__':
  main()