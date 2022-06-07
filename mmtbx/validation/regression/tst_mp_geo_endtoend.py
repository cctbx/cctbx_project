from __future__ import absolute_import, division, print_function
from mmtbx.validation.molprobity import mp_geo
import time
import libtbx.load_env
import os

def test_size_of_mp_geo_result():
  """This is an end-to-end test of mp_geo.
The sample file contains alternate conformations.
"""
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/2ERL_noH.pdb",
    test=os.path.isfile)
  args = ['pdb='+regression_pdb,
          'out_file=chiral_volume_validation_endtoend.out',
          'outliers_only=False',
          'bonds_and_angles=True']
  mp_geo.run(args)
  with open('chiral_volume_validation_endtoend.out', 'r') as f:
    lines = [i.rstrip('\n\r') for i in f.readlines()]
  count_bonds = 0
  count_angles = 0
  count_chirals = 0
  for line in lines:
    geom_identifier = line.split(':')[6]
    if "--" in geom_identifier:
      #2ERL_noH.pdb:  :   1: : :ASP:C--O:1.230:0.071:PROTEIN
      count_bonds += 1
    elif "-" in geom_identifier:
      #2ERL_noH.pdb:  :   1: : :ASP:C-CA-CB:111.672:0.827:PROTEIN
      count_angles += 1
    elif len(geom_identifier) == 2:
      #2ERL_noH.pdb:  :  21: :A:GLU:CA:2.514:0.020:PROTEIN
      count_chirals += 1
  assert count_bonds == 327
  assert count_angles == 442
  assert count_chirals == 49
  assert count_bonds + count_angles + count_chirals == len(lines),"there are lines in mp_geo output that are not bonds, angles, or chirals"

def test_size_of_mp_geo_result_outliers_only():
  """This is an end-to-end test of mp_geo using the outlier_only flag.
The sample file contains alternate conformations.
"""
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/2ERL_noH.pdb",
    test=os.path.isfile)
  args = ['pdb='+regression_pdb,
          'out_file=chiral_volume_validation_endtoend_outliers_only.out',
          'outliers_only=True',
          'bonds_and_angles=True']
  mp_geo.run(args)
  with open('chiral_volume_validation_endtoend_outliers_only.out', 'r') as f:
    lines = [i.rstrip('\n\r') for i in f.readlines()]
  count_bonds = 0
  count_angles = 0
  count_chirals = 0
  for line in lines:
    geom_identifier = line.split(':')[6]
    if "--" in geom_identifier:
      #2ERL_noH.pdb:  :   1: : :ASP:C--O:1.230:0.071:PROTEIN
      count_bonds += 1
    elif "-" in geom_identifier:
      #2ERL_noH.pdb:  :   1: : :ASP:C-CA-CB:111.672:0.827:PROTEIN
      count_angles += 1
    elif len(geom_identifier) == 2:
      #2ERL_noH.pdb:  :  21: :A:GLU:CA:2.514:0.020:PROTEIN
      count_chirals += 1
  assert count_bonds == 7
  assert count_angles == 5
  assert count_chirals == 0
  assert count_bonds + count_angles + count_chirals == len(lines),"there are lines in mp_geo outliers_only output that are not bonds, angles, or chirals"

def exercise():
  test_size_of_mp_geo_result()
  test_size_of_mp_geo_result_outliers_only()

if (__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("OK. Time: %8.3f"%(time.time()-t0))
