from __future__ import division
import libtbx.load_env
from libtbx.test_utils import approx_equal
from iotbx.command_line import merging_statistics
from cctbx.array_family import flex
import os
import sys
from cStringIO import StringIO

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
  assert ("""  1.81   1.74  12528   2073    6.04  97.10    1449.2     5.5  0.252  0.275  0.110  0.968""" in out.getvalue())
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
Crude resolution cutoff estimates:
  resolution of all data          :   2.000
  based on CC(1/2) > 0.5          : (use all data)
  based on mean(I/sigma) > 2      :   2.155
  based on R-merge < 0.5          :   2.372
  based on R-meas < 0.5           :   2.521
  based on completeness > 0.9     : (use all data)
  based on completeness > 0.5     : (use all data)""".splitlines() :
    assert line in out.getvalue()

if (__name__ == "__main__") :
  exercise(debug=("--debug" in sys.argv))
  print "OK"
