
from __future__ import division
from iotbx.command_line import merging_statistics
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.utils import Sorry
import libtbx.load_env
from cStringIO import StringIO
import os
import sys

def exercise (debug=False) :
  if (not libtbx.env.has_module("phenix_regression")) :
    print "phenix_regression not configured, skipping."
    return
  hkl_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/wizards/p9_se_w2.sca",
    test=os.path.isfile)
  args = [
    hkl_file,
    "space_group=I4",
    "unit_cell=113.949,113.949,32.474,90,90,90",
    "loggraph=True",
  ]
  if (debug) :
    args.append("debug=True")
    print " ".join(args)
  out = StringIO()
  result = merging_statistics.run(args, out=out)
  if (debug) :
    print out.getvalue()
  assert ("R-merge: 0.073" in out.getvalue())
  assert ("R-meas:  0.079" in out.getvalue())
  assert ("""  1.81   1.74  12528   2073    6.04  97.05    1449.2     5.5  0.252  0.275  0.110  0.968""" in out.getvalue()), out.getvalue()
  cif_block = result.as_cif_block()
  assert "_reflns_shell" in cif_block
  assert approx_equal(float(cif_block["_reflns.pdbx_Rpim_I_obs"]), result.overall.r_pim)
  assert approx_equal(float(cif_block["_reflns.phenix_cc_star"]), result.overall.cc_star)
  assert approx_equal(float(cif_block["_reflns.phenix_cc_1/2"]), result.overall.cc_one_half)
  assert approx_equal(
    flex.int(cif_block["_reflns_shell.number_measured_obs"]),
    [15737, 15728, 15668, 15371, 14996, 14771, 13899, 13549, 13206, 12528])
  assert "_reflns_shell.phenix_cc_star" in cif_block
  assert "_reflns_shell.phenix_cc_1/2" in cif_block
  remark_200 = result.as_remark_200(wavelength=0.9792).splitlines()
  assert ("REMARK 200  <I/SIGMA(I)> FOR SHELL         : 5.4942" in remark_200)
  assert ("REMARK 200  WAVELENGTH OR RANGE        (A) : 0.9792" in remark_200)
  # test resolution cutoffs
  args2 = list(args[:-1]) + ["high_resolution=2.5", "low_resolution=15"]
  out = StringIO()
  result = merging_statistics.run(args2, out=out)
  if (debug) :
    print out.getvalue()
  assert ("Resolution: 14.96 - 2.50" in out.getvalue())
  # extend binning
  args2 = list(args[:-1]) + ["high_resolution=1.5", "low_resolution=100",
    "--extend_d_max_min"]
  out = StringIO()
  result = merging_statistics.run(args2, out=out)
  if (debug) :
    print out.getvalue()
  assert ("Resolution: 100.00 - 1.50" in out.getvalue())
  assert ("  1.55   1.50      0      0    0.00   0.00       0.0     0.0   None   None   None  0.000  0.000""" in out.getvalue())
  # these should crash
  args2 = list(args[:-1]) + ["high_resolution=15", "low_resolution=2.5"]
  try :
    result = merging_statistics.run(args2, out=out)
  except Sorry, s :
    pass
  else :
    raise Exception_expected
  args2 = list(args[:-1]) + ["high_resolution=1.5", "low_resolution=1.6"]
  try :
    result = merging_statistics.run(args2, out=out)
  except Sorry, s :
    pass
  else :
    raise Exception_expected
  # change space group
  args = [
    hkl_file,
    "space_group=I422",
    "unit_cell=113.949,113.949,32.474,90,90,90",
    "loggraph=True",
  ]
  if (debug) :
    args.append("debug=True")
    print " ".join(args)
  out = StringIO()
  result = merging_statistics.run(args, out=out)
  if (debug) :
    print out.getvalue()
  assert ("28.49   3.76  15737   1224   12.86  99.84   47967.0    27.7  0.482  0.500  0.135  0.9" in out.getvalue()), out.getvalue()
  # exercise 2: estimate resolution cutoffs (and symmetry_file argument)
  hkl_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/harvesting/unmerged.sca",
    test=os.path.isfile)
  symm_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/harvesting/refine.mtz",
    test=os.path.isfile)
  args = [
    hkl_file,
    "symmetry_file=\"%s\"" % symm_file,
    "--estimate_cutoffs",
  ]
  out = StringIO()
  result = merging_statistics.run(args, out=out)
  if (debug) :
    print out.getvalue()
  for line in """\
  resolution of all data          :   2.000
  based on CC(1/2) >= 0.33        :   2.000
  based on mean(I/sigma) >= 2.0   :   2.155
  based on R-merge < 0.5          :   2.372
  based on R-meas < 0.5           :   2.520
  based on completeness >= 90%    :   2.000
  based on completeness >= 50%    :   2.000""".splitlines() :
    assert line in out.getvalue(), out.getvalue()

if (__name__ == "__main__") :
  exercise(debug=("--debug" in sys.argv))
  print "OK"
