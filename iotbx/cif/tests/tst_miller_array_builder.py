from __future__ import absolute_import, division, print_function
# Regression test for iotbx.cif.builders.miller_array_builder:
#  - the anomalous flag is computed once per distinct row selection per
#    reflection loop, not once per column: auto_anomalous() runs a full
#    Bijvoet-mate matching of the selected indices, and calling it for
#    every column dominated the build time of large structure-factor files;
#  - crystal_id / wavelength_id / scale_group_code columns are converted
#    once, not twice.
# Neither change may alter the built arrays: every array's anomalous flag
# is checked against miller.set(...).auto_anomalous() on its own indices.
import random
from cctbx import crystal, miller
from cctbx.array_family import flex
import iotbx.cif
from iotbx.cif import builders
from libtbx.test_utils import approx_equal

CRYSTAL_SYMMETRY = crystal.symmetry(
  unit_cell=(10, 20, 30, 90, 90, 90), space_group_symbol="P 21 21 21")

def make_cif(anomalous):
  """A one-block reflection cif with three distinct null patterns:
  none (intensities, ids, status), A (F and sigma F) and B (FWT and PHWT).
  In the anomalous case pattern A blanks the minus mate of every Bijvoet
  pair, so the F array alone holds no pairs while the loop does."""
  ms = miller.build_set(CRYSTAL_SYMMETRY, anomalous_flag=anomalous, d_min=4.0)
  n = ms.size()
  rnd = random.Random(0)
  f = [rnd.uniform(1, 100) for i in range(n)]
  if anomalous:
    asu, matches = ms.match_bijvoet_mates()
    missing_a = set(matches.pairs_hemisphere_selection("-"))
    assert len(missing_a) > 0
  else:
    missing_a = set(range(0, n, 5))
  missing_b = set(range(0, n, 7))
  lines = [
    "data_test",
    "_cell.length_a 10", "_cell.length_b 20", "_cell.length_c 30",
    "_cell.angle_alpha 90", "_cell.angle_beta 90", "_cell.angle_gamma 90",
    "_symmetry.space_group_name_H-M 'P 21 21 21'",
    "loop_",
    "_refln.crystal_id", "_refln.wavelength_id", "_refln.scale_group_code",
    "_refln.index_h", "_refln.index_k", "_refln.index_l",
    "_refln.status",
    "_refln.F_meas_au", "_refln.F_meas_sigma_au",
    "_refln.intensity_meas", "_refln.intensity_sigma",
    "_refln.pdbx_FWT", "_refln.pdbx_PHWT",
  ]
  for i, hkl in enumerate(ms.indices()):
    fa = ("%.10g %.10g" % (f[i], 0.1 * f[i])) if i not in missing_a else "? ?"
    fb = ("%.10g %.1f" % (f[i], 10.0 * (i % 36))) if i not in missing_b else "? ?"
    lines.append("1 1 1 %d %d %d %s %s %.10g %.10g %s" % (
      hkl[0], hkl[1], hkl[2], "f" if i % 10 == 0 else "o", fa, f[i] ** 2, f[i], fb))
  return "\n".join(lines) + "\n", f, missing_a, missing_b

def build_with_counters(cif_text):
  """Build the arrays while counting auto_anomalous() calls and column
  conversions (flex_std_string_as_miller_array)."""
  auto_anomalous_calls = []
  orig_auto = miller.set.auto_anomalous
  def counting_auto_anomalous(self, *args, **kw):
    auto_anomalous_calls.append(self.indices().size())
    return orig_auto(self, *args, **kw)
  n_conversions = [0]
  orig_conv = builders.miller_array_builder.flex_std_string_as_miller_array
  def counting_conversion(self, *args, **kw):
    n_conversions[0] += 1
    return orig_conv(self, *args, **kw)
  miller.set.auto_anomalous = counting_auto_anomalous
  builders.miller_array_builder.flex_std_string_as_miller_array = counting_conversion
  try:
    arrays = iotbx.cif.reader(input_string=cif_text).build_miller_arrays()["test"]
  finally:
    miller.set.auto_anomalous = orig_auto
    builders.miller_array_builder.flex_std_string_as_miller_array = orig_conv
  return arrays, auto_anomalous_calls, n_conversions[0]

def find_array(arrays, tag):
  hits = [a for a in arrays.values() if tag in a.info().label_string()]
  assert len(hits) == 1, (tag, [a.info().label_string() for a in arrays.values()])
  return hits[0]

def exercise(anomalous):
  cif_text, f, missing_a, missing_b = make_cif(anomalous)
  arrays, auto_anomalous_calls, n_conversions = build_with_counters(cif_text)
  n = len(f)
  labels = [a.info().label_string() for a in arrays.values()]
  assert len(arrays) == 7, labels  # 3 ids, status, F/sigF, I/sigI, FWT/PHWT

  # Parity: every array carries the flag auto_anomalous() gives on its own indices.
  for a in arrays.values():
    expected = miller.set(CRYSTAL_SYMMETRY, a.indices()).auto_anomalous().anomalous_flag()
    assert a.anomalous_flag() == expected, (a.info().label_string(), a.anomalous_flag(), expected)
  f_array = find_array(arrays, "_refln.F_meas_au")
  i_array = find_array(arrays, "_refln.intensity_meas")
  fwt_array = find_array(arrays, "_refln.pdbx_FWT")
  assert i_array.anomalous_flag() == anomalous
  assert f_array.anomalous_flag() == False  # pattern A removed every minus mate
  assert fwt_array.anomalous_flag() == anomalous
  assert i_array.size() == n
  assert f_array.size() == n - len(missing_a), (f_array.size(), n, len(missing_a))
  assert fwt_array.size() == n - len(missing_b), (fwt_array.size(), n, len(missing_b))
  assert approx_equal(list(f_array.data()), [f[i] for i in range(n) if i not in missing_a])
  assert approx_equal(list(f_array.sigmas()), [0.1 * f[i] for i in range(n) if i not in missing_a])
  assert approx_equal(list(i_array.data()), [f[i] ** 2 for i in range(n)])
  assert approx_equal(list(flex.abs(fwt_array.data())), [f[i] for i in range(n) if i not in missing_b])
  for tag in ("crystal_id", "wavelength_id", "scale_group_code", "status"):
    assert find_array(arrays, "_refln." + tag).size() == n, tag

  # One matching per distinct selection: none, A, B.
  assert len(auto_anomalous_calls) == 3, (len(auto_anomalous_calls), auto_anomalous_calls)
  assert sorted(auto_anomalous_calls) == sorted([n, n - len(missing_a), n - len(missing_b)]), \
    (auto_anomalous_calls, n, len(missing_a), len(missing_b))
  # One conversion per column (FWT, PHWT, F, sigF, I, sigI, 3 ids, status).
  assert n_conversions == 10, n_conversions

def run():
  exercise(anomalous=True)
  exercise(anomalous=False)
  print("OK")

if __name__ == "__main__":
  run()
