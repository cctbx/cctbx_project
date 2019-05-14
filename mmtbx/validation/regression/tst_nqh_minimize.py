
from __future__ import absolute_import, division, print_function
from mmtbx.validation.molprobity import nqh_minimize
from iotbx import pdb
from libtbx.utils import null_out
import time

pdb_str_1 = \
"""USER  MOD Single : A  49 GLN     :FLIP  amide:sc=   0.595  F(o=-3.2!,f=0.59)
USER  MOD Single : A  60 ASN    B:FLIP  amide:sc= -0.0622  F(o=-0.92,f=-0.062)
ATOM    321  N   ARG A  42      32.478  30.917  22.269  1.00  5.73           N
ATOM    322  CA  ARG A  42      31.200  30.329  22.780  1.00  6.97           C
ATOM    323  C   ARG A  42      30.210  30.509  21.650  1.00  7.15           C
ATOM    324  O   ARG A  42      29.978  31.726  21.269  1.00  7.33           O
ATOM    325  CB  ARG A  42      30.847  30.931  24.118  1.00 13.23           C
ATOM    326  CG  ARG A  42      29.412  30.796  24.598  1.00 21.27           C
ATOM    327  CD  ARG A  42      29.271  31.314  26.016  1.00 26.14           C
ATOM    328  NE  ARG A  42      27.875  31.317  26.443  1.00 32.26           N
ATOM    329  CZ  ARG A  42      27.132  32.423  26.574  1.00 34.32           C
ATOM    330  NH1 ARG A  42      27.630  33.656  26.461  1.00 35.30           N
ATOM    331  NH2 ARG A  42      25.810  32.299  26.732  1.00 36.39           N
ATOM      0  HA  ARG A  42      31.236  29.383  22.989  1.00  6.97           H   new
ATOM      0  HB2 ARG A  42      31.424  30.530  24.787  1.00 13.23           H   new
ATOM      0  HB3 ARG A  42      31.063  31.876  24.086  1.00 13.23           H   new
ATOM      0  HG2 ARG A  42      28.821  31.289  24.008  1.00 21.27           H   new
ATOM      0  HG3 ARG A  42      29.139  29.866  24.560  1.00 21.27           H   new
ATOM      0  HD2 ARG A  42      29.794  30.762  26.618  1.00 26.14           H   new
ATOM      0  HD3 ARG A  42      29.630  32.213  26.071  1.00 26.14           H   new
ATOM      0  HE  ARG A  42      27.507  30.561  26.622  1.00 32.26           H   new
ATOM      0 HH11 ARG A  42      28.467  33.768  26.296  1.00 35.30           H   new
ATOM      0 HH12 ARG A  42      27.114  34.338  26.553  1.00 35.30           H   new
ATOM      0 HH21 ARG A  42      25.450  31.518  26.748  1.00 36.39           H   new
ATOM      0 HH22 ARG A  42      25.320  33.000  26.817  1.00 36.39           H   new
ATOM    377  N   GLN A  49      23.880  26.727  23.851  1.00  8.89           N
ATOM    378  CA  GLN A  49      25.349  26.872  23.643  1.00  7.18           C
ATOM    379  C   GLN A  49      25.743  25.586  22.922  1.00  8.23           C
ATOM    380  O   GLN A  49      25.325  24.489  23.378  1.00  9.70           O
ATOM    381  CB  GLN A  49      26.070  27.025  24.960  1.00 11.67           C
ATOM    382  CG  GLN A  49      27.553  27.356  24.695  1.00 15.82           C
ATOM    383  CD  GLN A  49      28.262  27.576  26.020  1.00 20.21           C
ATOM    384  OE1 GLN A  49      27.777  28.585  26.739  1.00 23.23           O   flip
ATOM    385  NE2 GLN A  49      29.189  26.840  26.335  1.00 20.67           N   flip
ATOM      0  HA  GLN A  49      25.584  27.663  23.134  1.00  7.18           H   new
ATOM      0  HB2 GLN A  49      25.660  27.730  25.485  1.00 11.67           H   new
ATOM      0  HB3 GLN A  49      25.998  26.207  25.477  1.00 11.67           H   new
ATOM      0  HG2 GLN A  49      27.975  26.632  24.207  1.00 15.82           H   new
ATOM      0  HG3 GLN A  49      27.624  28.150  24.142  1.00 15.82           H   new
ATOM      0 HE21 GLN A  49      29.437  26.210  25.805  1.00 20.67           H   new
ATOM      0 HE22 GLN A  49      29.592  26.954  27.086  1.00 20.67           H   new
ATOM    467  N  BASN A  60      19.993  20.884  14.049  1.00 12.38           N
ATOM    468  CA BASN A  60      19.065  21.352  12.999  1.00 13.94           C
ATOM    469  C  BASN A  60      19.442  22.745  12.510  1.00 14.16           C
ATOM    470  O  BASN A  60      18.571  23.610  12.289  1.00 14.26           O
ATOM    471  CB BASN A  60      17.586  21.282  13.461  1.00 19.23           C
ATOM    472  CG BASN A  60      16.576  21.258  12.315  1.00 22.65           C
ATOM    473  OD1BASN A  60      16.924  20.586  11.216  1.00 25.45           O   flip
ATOM    474  ND2BASN A  60      15.440  21.819  12.378  1.00 24.09           N   flip
ATOM      0  H  BASN A  60      20.401  20.147  13.878  1.00 12.38           H   new
ATOM      0  HA BASN A  60      19.150  20.747  12.246  1.00 13.94           H   new
ATOM      0  HB2BASN A  60      17.464  20.487  14.004  1.00 19.23           H   new
ATOM      0  HB3BASN A  60      17.398  22.045  14.029  1.00 19.23           H   new
ATOM      0 HD21BASN A  60      15.213  22.253  13.085  1.00 24.09           H   new
ATOM      0 HD22BASN A  60      14.897  21.766  11.714  1.00 24.09           H   new
END
"""

def run_validation(pdb_file, ignore_hd=True):
  from mmtbx.validation import restraints
  import mmtbx.command_line
  cmdline = mmtbx.command_line.load_model_and_data(
    args=[pdb_file],
    master_phil=mmtbx.command_line.generic_simple_input_phil(),
    process_pdb_file=True,
    require_data=False,
    out=null_out())
  validation = restraints.combined(
    pdb_hierarchy=cmdline.pdb_hierarchy,
    xray_structure=cmdline.xray_structure,
    geometry_restraints_manager=cmdline.geometry,
    ignore_hd=ignore_hd)
  return validation

def exercise_nqh_minimize():
  open('modelFH.pdb', 'w').write(pdb_str_1)
  args = ['modelFH.pdb', 'modelFH_reg.pdb', '.']
  nqh_minimize.run(args)
  v1 = run_validation('modelFH.pdb', ignore_hd=True)
  v2 = run_validation('modelFH_reg.pdb', ignore_hd=True)
  assert v1.bonds.n_outliers == 3
  assert v2.bonds.n_outliers == 0

if (__name__ == "__main__"):
  t0 = time.time()
  exercise_nqh_minimize()
  print("OK. Time: %8.3f"%(time.time()-t0))
