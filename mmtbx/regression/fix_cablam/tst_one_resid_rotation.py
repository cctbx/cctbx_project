from __future__ import absolute_import, division, print_function

import iotbx.pdb
import mmtbx.model
from mmtbx.building.cablam_idealization import cablam_idealization, master_phil
import sys
import libtbx.load_env

pdb_str = """\
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
ATOM   2327  N   GLY A 318     169.195 115.930  63.690  1.00216.32           N
ATOM   2328  CA  GLY A 318     169.975 114.907  64.348  1.00193.16           C
ATOM   2329  C   GLY A 318     169.246 113.598  64.539  1.00197.19           C
ATOM   2330  O   GLY A 318     168.148 113.399  64.016  1.00193.16           O
ATOM   2331  N   GLN A 319     169.849 112.700  65.308  1.00184.03           N
ATOM   2332  CA  GLN A 319     169.232 111.415  65.589  1.00195.95           C
ATOM   2333  C   GLN A 319     169.246 111.137  67.080  1.00193.64           C
ATOM   2334  O   GLN A 319     168.185 111.047  67.708  1.00229.34           O
ATOM   2335  CB  GLN A 319     169.941 110.308  64.822  1.00201.09           C
ATOM   2336  CG  GLN A 319     169.719 110.407  63.336  1.00236.37           C
ATOM   2337  CD  GLN A 319     168.255 110.312  62.966  1.00254.36           C
ATOM   2338  OE1 GLN A 319     167.506 109.520  63.536  1.00280.71           O
ATOM   2339  NE2 GLN A 319     167.836 111.126  62.007  1.00220.80           N
ATOM   2340  N   ALA A 320     170.446 111.006  67.646  1.00140.99           N
ATOM   2341  CA  ALA A 320     170.595 110.942  69.090  1.00197.51           C
ATOM   2342  C   ALA A 320     169.906 109.734  69.704  1.00203.65           C
ATOM   2343  O   ALA A 320     168.789 109.863  70.203  1.00242.54           O
ATOM   2344  CB  ALA A 320     170.069 112.226  69.727  1.00240.45           C
ATOM   2345  N   LYS A 321     170.554 108.566  69.662  1.00164.18           N
ATOM   2346  CA  LYS A 321     169.963 107.306  70.104  1.00134.95           C
ATOM   2347  C   LYS A 321     169.103 107.477  71.344  1.00134.95           C
ATOM   2348  O   LYS A 321     167.904 107.194  71.302  1.00134.95           O
ATOM   2349  CB  LYS A 321     171.040 106.265  70.421  1.00145.47           C
ATOM   2350  CG  LYS A 321     171.950 105.868  69.279  1.00164.10           C
ATOM   2351  CD  LYS A 321     171.197 105.229  68.138  1.00145.47           C
ATOM   2352  CE  LYS A 321     172.173 104.777  67.070  1.00145.80           C
ATOM   2353  NZ  LYS A 321     171.487 104.149  65.918  1.00145.93           N
ATOM   2354  N   ARG A 322     169.682 107.900  72.454  1.00185.62           N
ATOM   2355  CA  ARG A 322     168.888 108.089  73.652  1.00142.57           C
ATOM   2356  C   ARG A 322     169.546 109.124  74.551  1.00128.56           C
ATOM   2357  O   ARG A 322     170.758 109.341  74.474  1.00128.56           O
ATOM   2358  CB  ARG A 322     168.719 106.769  74.369  1.00115.16           C
ATOM   2359  CG  ARG A 322     167.669 106.817  75.430  1.00133.21           C
ATOM   2360  CD  ARG A 322     167.578 105.605  76.270  1.00149.50           C
ATOM   2361  NE  ARG A 322     168.665 105.482  77.219  1.00115.16           N
ATOM   2362  CZ  ARG A 322     168.912 104.370  77.883  1.00115.16           C
ATOM   2363  NH1 ARG A 322     168.133 103.302  77.715  1.00116.90           N
ATOM   2364  NH2 ARG A 322     169.915 104.340  78.745  1.00115.16           N
ATOM   2365  N   VAL A 323     168.740 109.783  75.382  1.00121.29           N
ATOM   2366  CA  VAL A 323     169.198 110.884  76.220  1.00121.29           C
ATOM   2367  C   VAL A 323     168.668 110.712  77.632  1.00139.86           C
ATOM   2368  O   VAL A 323     167.480 110.430  77.835  1.00121.29           O
ATOM   2369  CB  VAL A 323     168.795 112.246  75.654  1.00117.56           C
ATOM   2370  CG1 VAL A 323     168.912 113.310  76.710  1.00151.93           C
ATOM   2371  CG2 VAL A 323     169.721 112.603  74.534  1.00121.07           C
"""

def exercise_no_sidechains(prefix="tst_one_resid_rotation_no_sidechains"):
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split('\n'), source_info=None)
  model = mmtbx.model.manager(
      model_input = pdb_inp)
  with open("%s_start.pdb" % prefix, 'w') as f:
    f.write(model.model_as_pdb())
  s = model.selection("name N or name CA or name C or name O")
  model = model.select(s)
  ci = cablam_idealization(model = model, params=master_phil.extract().cablam_idealization, log=sys.stdout)
  pdb_txt = model.model_as_pdb()

def exercise_yes_sidechains(prefix="tst_one_resid_rotation_yes_sidechains"):
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split('\n'), source_info=None)
  model = mmtbx.model.manager(
      model_input = pdb_inp)
  with open("%s_start.pdb" % prefix, 'w') as f:
    f.write(model.model_as_pdb())
  ci = cablam_idealization(model = model, params=master_phil.extract().cablam_idealization, log=sys.stdout)
  pdb_txt = model.model_as_pdb()

if __name__ == '__main__':
  if (not libtbx.env.has_module(name="probe")):
    print("Skipping: probe not configured")
  else:
    exercise_no_sidechains()
    exercise_yes_sidechains()
