
from __future__ import absolute_import, division, print_function
from iotbx.command_line import merging_statistics
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal, Exception_expected, show_diff
from libtbx.utils import Sorry
import libtbx.load_env
from six.moves import cStringIO as StringIO
import os
import sys


def exercise(debug=False):
  if not libtbx.env.has_module("phenix_regression"):
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
  if debug:
    args.append("debug=True")
    print(" ".join(args))
  out = StringIO()
  result = merging_statistics.run(args, out=out)
  if debug:
    print(out.getvalue())
  assert ("R-merge: 0.073" in out.getvalue())
  assert ("R-meas:  0.079" in out.getvalue())
  assert ("""  1.81   1.74  12528   2073    6.04  97.05    1449.2     5.2    0.252    0.275    0.110    0.294   0.967   0.281""" in out.getvalue()), out.getvalue()
  cif_block = result.as_cif_block()
  assert "_reflns_shell" in cif_block
  assert approx_equal(float(cif_block["_reflns.pdbx_Rpim_I_all"]), result.overall.r_pim)
  assert approx_equal(float(cif_block["_reflns.pdbx_CC_half"]), result.overall.cc_one_half)
  assert approx_equal(float(cif_block["_reflns.percent_possible_obs"]), result.overall.completeness * 100.0)
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
  merging_statistics.run(args2, out=out)
  if debug:
    print(out.getvalue())
  assert ("Resolution: 14.96 - 2.50" in out.getvalue())
  # extend binning
  args2 = list(args[:-1]) + ["high_resolution=1.5", "low_resolution=100",
    "--extend_d_max_min"]
  out = StringIO()
  merging_statistics.run(args2, out=out)
  if debug:
    print(out.getvalue())
  assert ("Resolution: 100.00 - 1.50" in out.getvalue())
  assert ("  1.55   1.50      0      0    0.00   0.00       0.0     0.0     None     None     None     0.000   0.000""" in out.getvalue())
  args2 = args + ["json.file_name=merging_stats.json", "json.indent=2",
                  "mmcif.file_name=merging_stats.mmcif", "mmcif.data_name=test"]
  out = StringIO()
  merging_statistics.run(args2, out=out)
  assert os.path.exists("merging_stats.json")
  with open("merging_stats.json", "rb") as f:
    import json
    d = json.load(f)
    list(d.keys())  #FIXME why am I here?
    expected_keys = [
      'n_obs', 'd_star_sq_max', 'i_over_sigma_mean', 'completeness',
      'cc_one_half', 'r_meas', 'd_star_sq_min', 'cc_anom', 'r_pim',
      'i_mean', 'cc_one_half_critical_value', 'r_merge', 'multiplicity',
      'cc_one_half_significance', 'n_uniq', 'anom_completeness', 'anom_signal',
      'delta_i_mean_over_sig_delta_i_mean',
    ]
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
  #
  args3 = args + ["anomalous=True"]
  result = merging_statistics.run(args3, out=out)
  assert approx_equal(result.overall.anom_signal, 0.0835294960933567)
  assert approx_equal(result.overall.anom_completeness, 0.9836098441180893)
  assert approx_equal(result.overall.delta_i_mean_over_sig_delta_i_mean, 1.1764137800824204)
  # these should crash
  args2 = list(args[:-1]) + ["high_resolution=15", "low_resolution=2.5"]
  try :
    merging_statistics.run(args2, out=out)
  except Sorry:
    pass
  else :
    raise Exception_expected
  args2 = list(args[:-1]) + ["high_resolution=1.5", "low_resolution=1.6"]
  try :
    merging_statistics.run(args2, out=out)
  except Sorry:
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
  if debug:
    args.append("debug=True")
    print(" ".join(args))
  out = StringIO()
  merging_statistics.run(args, out=out)
  if debug:
    print(out.getvalue())
  assert (" 28.49   3.76  15737   1224   12.86  99.84   47967.0    11.6    0.482    0.500    0.135    0.136   0.973  -0.513" in out.getvalue()), out.getvalue()
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
  merging_statistics.run(args, out=out)
  if debug:
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

  # Test that merged anomalous data are already merged.
  try:
    merging_statistics.run([hkl_file, "anomalous=True"], out=out)
  except Sorry as e:
    assert str(e) == 'The data in i_anomalous(+),SIGi_anomalous(+),i_anomalous(-),SIGi_anomalous(-) are already merged.  Only unmerged (but scaled) data may be used in this program.'
  else: raise Exception_expected

  # Test that merged anomalous data can still be merged as non-anomalous data.
  merging_statistics.run([hkl_file], out=out)

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
  assert len(result.bins) == 20
  assert approx_equal(result.overall.d_max, 52.445602416992195)
  assert result.overall.n_obs == 119045, result.overall.n_obs
  assert approx_equal(result.bins[0].cc_anom, 0.879550520045)
  assert approx_equal(result.bins[0].r_anom, 0.12112581576975405)
  assert result.bins[0].cc_anom_significance is True
  assert approx_equal(result.bins[0].cc_anom_critical_value, 0.0873548986308)
  assert approx_equal(result.cc_one_half_overall, 0.931122967496)
  assert approx_equal(result.cc_one_half_sigma_tau_overall, 0.9280192675969664)
  assert approx_equal(result.bins[0].cc_one_half, 0.9969293192434535)
  assert approx_equal(result.bins[0].cc_one_half_sigma_tau, 0.9968045160775104)
  assert result.bins[0].cc_one_half_significance is True
  assert result.bins[0].cc_one_half_sigma_tau_significance is True
  assert approx_equal(result.bins[-1].cc_one_half, 0.675340867481686)
  assert approx_equal(result.bins[-1].cc_one_half_sigma_tau, 0.6711734115834956)
  assert result.bins[-1].cc_one_half_significance is True
  assert result.bins[-1].cc_one_half_sigma_tau_significance is True
  #
  args2 = args + ["anomalous=True"]
  result = merging_statistics.run(args2, out=out)
  app = result.overall.anom_probability_plot_all_data
  assert approx_equal(app.slope, 1.6747371892218779)
  assert approx_equal(app.intercept, 0.03443233016003127)
  assert app.n_pairs == 24317
  assert app.expected_delta is None
  app = result.overall.anom_probability_plot_expected_delta
  assert approx_equal(app.slope, 1.2190614967160294)
  assert approx_equal(app.intercept, 0.031151653693628906)
  assert app.n_pairs == 15365
  assert app.expected_delta == 0.9
  d = result.overall.as_dict()
  for k in ("anom_probability_plot_all_data", "anom_probability_plot_expected_delta"):
    assert list(d[k].keys()) == ["slope", "intercept", "n_pairs", "expected_delta"], list(d[k].keys())
  out = StringIO()
  result.overall.show_anomalous_probability_plot(out)
  assert not show_diff(out.getvalue(),
                       """\
Anomalous probability plot (all data):
  slope:     1.675
  intercept: 0.034
  n_pairs:   24317
Anomalous probability plot (expected delta = 0.9):
  slope:     1.219
  intercept: 0.031
  n_pairs:   15365
""")

  args = [hkl_file, "binning_method=counting_sorted",
          "eliminate_sys_absent=False", "n_bins=20",
          "reflections_per_bin=2000",
          "cc_one_half_significance_level=0.01"]
  result = merging_statistics.run(args, out=out)
  assert len(result.bins) == 12


if __name__ == "__main__":
  exercise(debug=("--debug" in sys.argv))
  print("OK")
