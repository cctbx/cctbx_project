
from __future__ import division
from mmtbx.validation.molprobity import mp_geo
from iotbx import pdb
import time

pdb_str_1 = """
ATOM     60  N   GLU A   9       7.757   9.623  27.829  1.00 12.74           N
ATOM     61  CA  GLU A   9       7.695  10.769  26.927  1.00 12.72           C
ATOM     62  C   GLU A   9       7.351  12.056  27.684  1.00 12.59           C
ATOM     63  O   GLU A   9       6.321  12.134  28.348  1.00 14.32           O
ATOM     64  CB  GLU A   9       6.633  10.459  25.859  1.00 13.52           C
ATOM     65  CG  GLU A   9       6.323  11.587  24.899  1.00 15.93           C
ATOM     66  CD  GLU A   9       7.533  12.006  24.105  1.00 16.48           C
ATOM     67  OE1 GLU A   9       8.059  11.171  23.337  1.00 19.96           O
ATOM     68  OE2 GLU A   9       7.960  13.164  24.257  1.00 18.90           O
ATOM     69  N   ASP A  10       8.216  13.062  27.589  1.00 14.68           N
ATOM     70  CA  ASP A  10       7.992  14.328  28.283  1.00 15.47           C
ATOM     71  C   ASP A  10       6.987  15.255  27.598  1.00 15.82           C
ATOM     72  O   ASP A  10       6.315  16.049  28.258  1.00 15.43           O
ATOM     73  CB  ASP A  10       9.322  15.075  28.453  1.00 16.97           C
ATOM     74  CG AASP A  10       9.139  16.482  28.993  0.50 16.48           C
ATOM     75  CG BASP A  10       9.941  15.472  27.128  0.50 17.05           C
ATOM     76  OD1AASP A  10       9.189  17.442  28.197  0.50 17.34           O
ATOM     77  OD1BASP A  10      11.349  14.573  26.359  0.50 18.03           O
ATOM     78  OD2AASP A  10       8.934  16.630  30.215  0.50 17.88           O
ATOM     79  OD2BASP A  10      10.016  16.689  26.852  0.50 20.20           O
"""

def exercise_mp_geo():
  open('mp_geo.pdb', 'w').write(pdb_str_1)
  args = ['pdb=mp_geo.pdb',
          'out_file=mp_geo.out',
          'outliers_only=False',
          'bonds_and_angles=True']
  mp_geo.run(args)
  f = file('mp_geo.out', 'rb')
  lines = f.readlines()
  assert 'mp_geo.pdb: A:  10: :B:ASP:CG--OD1:1.839:31.054:PROTEIN\n' in lines
  assert \
    'mp_geo.pdb: A:  10: :B:ASP:OD1-CG-OD2:109.733:5.486:PROTEIN\n' in lines
  f.close()

if (__name__ == "__main__"):
  t0 = time.time()
  exercise_mp_geo()
  print "OK. Time: %8.3f"%(time.time()-t0)
