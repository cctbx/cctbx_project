
from __future__ import division, print_function
from iotbx.command_line import merging_statistics
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.utils import Sorry
import libtbx.load_env
from cStringIO import StringIO
import os
import sys

def exercise(debug=False):
  if (not libtbx.env.has_module("phenix_regression")):
    print("phenix_regression not configured, skipping.")
    return
  hkl_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/wizards/help_tests/test_help/p9_se_w2.sca",
    test=os.path.isfile)
  assert hkl_file is not None
  args = [
    hkl_file,
    "space_group=I4",
    "unit_cell=113.949,113.949,32.474,90,90,90",
    "loggraph=True",
  ]
  if (debug):
    args.append("debug=True")
    print(" ".join(args))
  out = StringIO()
  result = merging_statistics.run(args, out=out)
  if (debug):
    print(out.getvalue())
  assert ("R-merge: 0.073" in out.getvalue())
  assert ("R-meas:  0.079" in out.getvalue())
  assert ("""  1.81   1.74  12528   2073    6.04  97.05    1449.2     5.2    0.252    0.275    0.110   0.967   0.281""" in out.getvalue()), out.getvalue()
  cif_block = result.as_cif_block()
  assert "_reflns_shell" in cif_block
  assert approx_equal(float(cif_block["_reflns.pdbx_Rpim_I_obs"]), result.overall.r_pim)
  assert approx_equal(float(cif_block["_reflns.pdbx_CC_half"]), result.overall.cc_one_half)
  assert approx_equal(
    flex.int(cif_block["_reflns_shell.number_measured_obs"]),
    [15737, 15728, 15668, 15371, 14996, 14771, 13899, 13549, 13206, 12528])
  assert "_reflns_shell.pdbx_CC_half" in cif_block
  remark_200 = result.as_remark_200(wavelength=0.9792).splitlines()
  assert ("REMARK 200  <I/SIGMA(I)> FOR SHELL         : 5.1536" in remark_200),"\n".join(remark_200)
  assert ("REMARK 200  WAVELENGTH OR RANGE        (A) : 0.9792" in remark_200)
  # test resolution cutoffs
  args2 = list(args[:-1]) + ["high_resolution=2.5", "low_resolution=15"]
  out = StringIO()
  result = merging_statistics.run(args2, out=out)
  if (debug):
    print(out.getvalue())
  assert ("Resolution: 14.96 - 2.50" in out.getvalue())
  # extend binning
  args2 = list(args[:-1]) + ["high_resolution=1.5", "low_resolution=100",
    "--extend_d_max_min"]
  out = StringIO()
  result = merging_statistics.run(args2, out=out)
  if (debug):
    print(out.getvalue())
  assert ("Resolution: 100.00 - 1.50" in out.getvalue())
  assert ("  1.55   1.50      0      0    0.00   0.00       0.0     0.0     None     None     None   0.000   0.000""" in out.getvalue())
  args2 = args + ["json.file_name=merging_stats.json", "json.indent=2",
                  "mmcif.file_name=merging_stats.mmcif", "mmcif.data_name=test"]
  out = StringIO()
  result = merging_statistics.run(args2, out=out)
  assert os.path.exists("merging_stats.json")
  with open("merging_stats.json", "rb") as f:
    import json
    d = json.load(f)
    d.keys()
    expected_keys = [
      'n_obs', 'd_star_sq_max', 'i_over_sigma_mean', 'completeness',
      'cc_one_half', 'r_meas', 'd_star_sq_min', 'cc_anom', 'r_pim',
      'i_mean', 'cc_one_half_critical_value', 'r_merge', 'multiplicity',
      'cc_one_half_significance', 'n_uniq']
    for k in expected_keys:
      assert k in d
    assert approx_equal(
      d['i_over_sigma_mean'],
      [20.66560600258035, 20.345307828448572, 18.349521604862186,
       16.032340432560876, 15.111031030727927, 12.289599369298855,
       10.059003053756124, 8.042714696129208, 6.297114205813105,
       5.153613598686054]
    )
    assert approx_equal(
      d['n_obs'],
      [15737, 15728, 15668, 15371, 14996, 14771, 13899, 13549, 13206, 12528]
    )
    assert 'overall' in d
    assert approx_equal(d['overall']['i_mean'], 18672.03367141032)
    assert approx_equal(d['overall']['multiplicity'], 6.747367444449599)
  import iotbx.cif
  assert os.path.exists("merging_stats.mmcif")
  cif = iotbx.cif.reader("merging_stats.mmcif").model()
  assert approx_equal(float(cif["test"]["_reflns.pdbx_CC_half"]), 0.997601585888)
  assert list(cif["test"]["_reflns_shell.number_measured_obs"]) == ['15737', '15728', '15668', '15371', '14996', '14771', '13899', '13549', '13206', '12528']
  # these should crash
  args2 = list(args[:-1]) + ["high_resolution=15", "low_resolution=2.5"]
  try :
    result = merging_statistics.run(args2, out=out)
  except Sorry as s :
    pass
  else :
    raise Exception_expected
  args2 = list(args[:-1]) + ["high_resolution=1.5", "low_resolution=1.6"]
  try :
    result = merging_statistics.run(args2, out=out)
  except Sorry as s :
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
  if (debug):
    args.append("debug=True")
    print(" ".join(args))
  out = StringIO()
  result = merging_statistics.run(args, out=out)
  if (debug):
    print(out.getvalue())
  assert (" 28.49   3.76  15737   1224   12.86  99.84   47967.0    11.6    0.482    0.500    0.135   0.973  -0.513" in out.getvalue()), out.getvalue()
  # exercise 2: estimate resolution cutoffs (and symmetry_file argument)
  hkl_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/harvesting/unmerged.sca",
    test=os.path.isfile)
  assert hkl_file is not None
  symm_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/harvesting/refine.mtz",
    test=os.path.isfile)
  assert symm_file is not None
  args = [
    hkl_file,
    "symmetry_file=\"%s\"" % symm_file,
    "--estimate_cutoffs",
  ]
  out = StringIO()
  result = merging_statistics.run(args, out=out)
  if (debug):
    print(out.getvalue())
  for line in """\
  resolution of all data          :   2.000
  based on CC(1/2) >= 0.33        :   2.000
  based on mean(I/sigma) >= 2.0   :   2.253
  based on R-merge < 0.5          :   2.372
  based on R-meas < 0.5           :   2.520
  based on completeness >= 90%    :   2.000
  based on completeness >= 50%    :   2.000""".splitlines():
    assert line in out.getvalue(), out.getvalue()
  # check suitable error emitted given merged input mtz (not containing M_ISYM column)
  out = StringIO()
  hkl_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/i_anomalous.mtz",
    test=os.path.isfile)
  assert hkl_file is not None
  args = [hkl_file]
  try:
    result = merging_statistics.run(args, out=out)
  except Sorry as e:
    assert str(e) == 'The data in i_anomalous(+),SIGi_anomalous(+),i_anomalous(-),SIGi_anomalous(-) are already merged.  Only unmerged (but scaled) data may be used in this program.'
  else: raise Exception_expected
  # test use_internal_variance option
  out = StringIO()
  hkl_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/wizards/help_tests/test_help/unmerged.mtz",
    test=os.path.isfile)
  assert hkl_file is not None
  args = [hkl_file, "use_internal_variance=False"]
  result = merging_statistics.run(args, out=out)
  assert approx_equal(result.overall.i_over_sigma_mean, 4.4867598237199)
  assert approx_equal(result.overall.unmerged_i_over_sigma_mean, 2.559049577429115)
  args = [hkl_file,
          #"use_internal_variance=True" # this is the default behaviour
          ]
  result = merging_statistics.run(args, out=out)
  assert approx_equal(result.overall.i_over_sigma_mean, 4.063574245292925)
  assert approx_equal(result.overall.unmerged_i_over_sigma_mean, 2.559049577429115)
  # test eliminate_sys_absent option
  out = StringIO()
  hkl_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/AUTOMATIC_DEFAULT_scaled_unmerged_WAVE1.mtz",
    test=os.path.isfile)
  assert hkl_file is not None
  args = [hkl_file,
          #"eliminate_sys_absent=True" # default behaviour
          ]
  result = merging_statistics.run(args, out=out)
  assert approx_equal(result.overall.d_max, 43.069972142418365)
  assert result.overall.n_obs == 118981, result.overall.n_obs
  args = [hkl_file, "binning_method=counting_sorted",
          "eliminate_sys_absent=False", "n_bins=20",
          "cc_one_half_significance_level=0.01"]
  result = merging_statistics.run(args, out=out)
  assert approx_equal(result.overall.d_max, 52.445602416992195)
  assert result.overall.n_obs == 119045, result.overall.n_obs
  assert approx_equal(result.bins[0].cc_anom, 0.8445946103668791)
  assert result.bins[0].cc_anom_significance is True
  assert approx_equal(result.bins[0].cc_anom_critical_value, 0.08001889224417998)
  assert approx_equal(result.cc_one_half_overall, 0.9277155216469152)
  assert approx_equal(result.cc_one_half_sigma_tau_overall, 0.9334218927846825)
  assert approx_equal(result.bins[0].cc_one_half, 0.9965897219959194)
  assert approx_equal(result.bins[0].cc_one_half_sigma_tau, 0.9971174078562669)
  assert result.bins[0].cc_one_half_significance is True
  assert result.bins[0].cc_one_half_sigma_tau_significance is True
  assert approx_equal(result.bins[-1].cc_one_half, 0.716558466728842)
  assert approx_equal(result.bins[-1].cc_one_half_sigma_tau, 0.772170987527236)
  assert result.bins[-1].cc_one_half_significance is True
  assert result.bins[-1].cc_one_half_sigma_tau_significance is True

if (__name__ == "__main__"):
  exercise(debug=("--debug" in sys.argv))
  print("OK")
