import json
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

# Start tests

dm = DataManager()
dm.process_model_str("1yjp",pdb_str_1yjp)
model= dm.get_model()

def tst_01():
  # test for no crashes on id_str geo
  container = GeoParseContainer()
  container.parse_str(tst_1_geo)
  entries = container.entries_list
  records = container.records_list
  assert len(entries) == len(records)
  assert len(entries) ==30


def tst_02():
  # test for no crashes on i_seq geo
  container = GeoParseContainer()
  # replace id_strs with ints as a substitute for i_seqs
  tst_1_geo_int = container._replace_idstr_with_int(tst_1_geo,max_int=50)
  container.parse_str(tst_1_geo_int)
  records = container.records
  assert len(container.entries_list) == 30
  entries = container.entries_list
  

  # test add model
  container.model = model
  records = container.records_list

  # test build proxies
  proxies = container.proxies_list

  assert len(entries) ==30
  assert len(entries) == len(records)

  # nonbonded entries are not converted to proxies. Subtract them
  assert len(records)-len(container.entries["nonbonded"]) == len(proxies), (
    f"{len(records)}, {len(proxies)}"
  )

def main():
  tst_01()
  tst_02()
  print('OK')

if __name__ == '__main__':
  main()