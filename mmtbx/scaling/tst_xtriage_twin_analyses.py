
# XXX This overlaps heavily with tst_scaling.py, which is much more
# comprehensive.  Maybe remove this one?

from __future__ import absolute_import, division, print_function
from mmtbx.scaling import twin_analyses
import iotbx.pdb
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal, show_diff
from libtbx.utils import null_out
from six.moves import cStringIO as StringIO
import sys

def exercise_twin_detection(verbose=False):
  # simple model with translational pseudosymmetry
  # XXX one big disadvantage: no centric reflections!
  pdb_str = """\
CRYST1   12.000    8.000   12.000  90.02  89.96  90.05 P 1           1
ATOM     39  N   ASN A   6       5.514   2.664   4.856  1.00 11.99           N
ATOM     40  CA  ASN A   6       6.831   2.310   4.318  1.00 12.30           C
ATOM     41  C   ASN A   6       7.854   2.761   5.324  1.00 13.40           C
ATOM     42  O   ASN A   6       8.219   3.943   5.374  1.00 13.92           O
ATOM     43  CB  ASN A   6       7.065   3.016   2.993  1.00 12.13           C
ATOM     44  CG  ASN A   6       5.961   2.735   2.003  1.00 12.77           C
ATOM     45  OD1 ASN A   6       5.798   1.604   1.551  1.00 14.27           O
ATOM     46  ND2 ASN A   6       5.195   3.747   1.679  1.00 10.07           N
ATOM     47  N   TYR A   7       8.292   1.817   6.147  1.00 14.70           N
ATOM     48  CA  TYR A   7       9.159   2.144   7.299  1.00 15.18           C
ATOM     49  C   TYR A   7      10.603   2.331   6.885  1.00 15.91           C
ATOM     50  O   TYR A   7      11.041   1.811   5.855  1.00 15.76           O
ATOM     51  CB  TYR A   7       9.061   1.065   8.369  1.00 15.35           C
ATOM     52  CG  TYR A   7       7.665   0.929   8.902  1.00 14.45           C
ATOM     53  CD1 TYR A   7       6.771   0.021   8.327  1.00 15.68           C
ATOM     54  CD2 TYR A   7       7.210   1.756   9.920  1.00 14.80           C
ATOM     55  CE1 TYR A   7       5.480  -0.094   8.796  1.00 13.46           C
ATOM     56  CE2 TYR A   7       5.904   1.649  10.416  1.00 14.33           C
ATOM     57  CZ  TYR A   7       5.047   0.729   9.831  1.00 15.09           C
ATOM     58  OH  TYR A   7       3.766   0.589  10.291  1.00 14.39           O
ATOM     59  OXT TYR A   7      11.358   2.999   7.612  1.00 17.49           O
TER
ATOM      1  N   ASN B   6       1.414   5.113   6.019  1.00 12.99           N
ATOM      2  CA  ASN B   6       2.720   4.776   5.445  1.00 13.30           C
ATOM      3  C   ASN B   6       3.763   5.209   6.438  1.00 14.40           C
ATOM      4  O   ASN B   6       4.125   6.391   6.507  1.00 14.92           O
ATOM      5  CB  ASN B   6       2.922   5.513   4.131  1.00 13.13           C
ATOM      6  CG  ASN B   6       1.798   5.250   3.160  1.00 13.77           C
ATOM      7  OD1 ASN B   6       1.629   4.129   2.686  1.00 15.27           O
ATOM      8  ND2 ASN B   6       1.022   6.266   2.875  1.00 11.07           N
ATOM      9  N   TYR B   7       4.222   4.248   7.230  1.00 15.70           N
ATOM     10  CA  TYR B   7       5.113   4.552   8.370  1.00 16.18           C
ATOM     11  C   TYR B   7       6.547   4.754   7.929  1.00 16.91           C
ATOM     12  O   TYR B   7       6.964   4.259   6.878  1.00 16.76           O
ATOM     13  CB  TYR B   7       5.042   3.449   9.417  1.00 16.35           C
ATOM     14  CG  TYR B   7       3.659   3.296   9.977  1.00 15.45           C
ATOM     15  CD1 TYR B   7       2.756   2.398   9.402  1.00 16.68           C
ATOM     16  CD2 TYR B   7       3.224   4.098  11.023  1.00 15.80           C
ATOM     17  CE1 TYR B   7       1.476   2.267   9.896  1.00 14.46           C
ATOM     18  CE2 TYR B   7       1.929   3.975  11.545  1.00 15.33           C
ATOM     19  CZ  TYR B   7       1.063   3.065  10.959  1.00 16.09           C
ATOM     20  OH  TYR B   7      -0.207   2.910  11.443  1.00 15.39           O
ATOM     21  OXT TYR B   7       7.316   5.408   8.654  1.00 18.49           O
END
"""
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str)
  xrs = pdb_in.xray_structure_simple()
  fc = abs(xrs.structure_factors(d_min=2.5).f_calc())
  fc = fc.set_observation_type_xray_amplitude()
  sigf = flex.double(fc.size(), 0.1) + (fc.data() * 0.03)
  fc = fc.customized_copy(sigmas=sigf)
  # and now add twinning
  fc_twin = fc.twin_data(twin_law='-l,-k,-h', alpha=0.4)
  fc_twin.as_mtz_dataset(column_root_label="F").mtz_object().write(
    "tmp_xtriage_twinned.mtz")
  # twin_laws
  laws = twin_analyses.twin_laws(miller_array=fc_twin)
  assert (len(laws.operators) == 7)
  delta_santoro = [ tl.delta_santoro for tl in laws.operators ]
  assert approx_equal(delta_santoro,
    [0.104655, 0.104655, 0.066599, 0.104655, 0.076113, 0.066599, 0.028542])
  delta_le_page = [ tl.delta_le_page for tl in laws.operators ]
  assert approx_equal(delta_le_page,
    [0.053839, 0.053839, 0.053839, 0.064020, 0.044706, 0.049480, 0.021221])
  delta_lebedev = [ tl.delta_lebedev for tl in laws.operators ]
  assert approx_equal(delta_lebedev,
    [0.000787, 0.000787, 0.000767, 0.000912, 0.000637, 0.000705, 0.000302])
  # TNCS
  tps = twin_analyses.detect_pseudo_translations(
    miller_array=fc_twin, out=null_out())
  assert ([tps.mod_h, tps.mod_k, tps.mod_l] == [3,2,3])
  assert approx_equal(tps.high_peak, 47.152352)
  assert approx_equal(tps.high_p_value, 0.000103)
  # L test
  fc_norm_acentric = twin_analyses.wilson_normalised_intensities(
    miller_array=fc_twin, normalise=True, out=null_out()).acentric
  ltest = twin_analyses.l_test(miller_array=fc_norm_acentric,
    parity_h=3, parity_l=3)
  assert approx_equal(ltest.mean_l, 0.498974)
  assert approx_equal(ltest.mean_l2, 0.371674)
  # we need to go through the wrapper class for a lot of these...
  out = StringIO()
  tw = twin_analyses.twin_analyses(miller_array=fc_twin, out=out)
  tw2 = twin_analyses.twin_analyses(miller_array=fc_twin, out=out,
    d_hkl_for_l_test=[3,2,3])
  # Wilson moments
  wm = tw.wilson_moments
  assert ([wm.centric_i_ratio, wm.centric_f_ratio, wm.centric_e_sq_minus_one]
          == [None, None, None])
  assert approx_equal(wm.acentric_i_ratio, 2.878383)
  assert approx_equal(wm.acentric_f_ratio, 0.717738)
  assert approx_equal(wm.acentric_abs_e_sq_minus_one, 0.808899)
  assert (not wm.centric_present)
  # Overall analyses
  out = StringIO()
  tw.show(out=out)
  assert ("significantly different than is expected" in out.getvalue())
  summary = tw.twin_summary
  assert approx_equal(summary.l_mean, 0.4989738)
  assert approx_equal(summary.patterson_height, 47.152352)
  for k, twin_r_factor in enumerate(summary.r_obs):
    if (k == 6) : assert twin_r_factor < 0.11
    else : assert twin_r_factor > 0.3
  out2 = StringIO()
  tw.l_test.show(out=out2)
  assert ("""\
  Mean |L|   :0.499  (untwinned: 0.500; perfect twin: 0.375)
  Mean  L^2  :0.372  (untwinned: 0.333; perfect twin: 0.200)
""" in out2.getvalue()), out2.getvalue()
  out3 = StringIO()
  tw2.l_test.show(out=out3)
  assert not show_diff(out2.getvalue(), out3.getvalue())
  if (verbose):
    print(out.getvalue())
  # twin_results_interpretation object via cctbx.miller.array API extension
  # XXX I get slightly different numbers here versus running through the
  # twin_analyses call above - this seems to be caused by the resolution
  # limits passed here.  Need to confirm these with PHZ.
  result = fc_twin.analyze_intensity_statistics()
  assert approx_equal(result.maha_l, 27.580674)
  assert approx_equal(result.r_obs,
    [0.514901, 0.514901, 0.397741, 0.580353, 0.579294, 0.420914, 0.114999])
  assert approx_equal(result.h_alpha,
    [0.019980, 0.019980, 0.211788, 0.073926, 0.079920, 0.150849, 0.389610])
  assert result.has_twinning()

if (__name__ == "__main__"):
  exercise_twin_detection(verbose=("--verbose" in sys.argv))
  print("OK")
