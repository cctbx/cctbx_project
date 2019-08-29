from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
import mmtbx.hydrogens
from libtbx.utils import null_out


def run():
  correct_H_position_with_cdl()


def compare_models(pdb_str,
                   contains     = None,
                   not_contains = None,
                   number_h     = None):
  #
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model_initial = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  xrs = model_initial.get_xray_structure()
  hd_sel_initial = model_initial.get_hd_selection()
  number_h_expected = hd_sel_initial.count(True)

  model_without_h = model_initial.select(~hd_sel_initial)
  hd_sel_without_h = model_without_h.get_hd_selection()
  assert (hd_sel_without_h is not None)
  assert (hd_sel_without_h.count(True) == 0)

  model_h_added = mmtbx.hydrogens.add(model = model_without_h)
  hd_sel_h_added = model_h_added.get_hd_selection()
  number_h_added = hd_sel_h_added.count(True)
  ph_initial = model_initial.get_hierarchy()
  ph_h_added = model_h_added.get_hierarchy()
  assert ph_initial.is_similar_hierarchy(other=ph_h_added)

  if number_h:
    assert(number_h == number_h_added)

  if not_contains:
    h_atoms_added = model_h_added.get_hierarchy().select(hd_sel_h_added).atoms()
    h_names_added = list(h_atoms_added.extract_name())
    assert (not_contains not in h_names_added)

  if contains:
    h_atoms_added = model_h_added.get_hierarchy().select(hd_sel_h_added).atoms()
    h_names_added = list(h_atoms_added.extract_name())
    assert (contains in h_names_added)


def correct_H_position_with_cdl():
  compare_models(pdb_str = pdb_str_1)



pdb_str_1 = """
CRYST1   72.240   72.010   86.990  90.00  90.00  90.00 P 21 21 21
SCALE1      0.013843  0.000000  0.000000        0.00000
SCALE2      0.000000  0.013887  0.000000        0.00000
SCALE3      0.000000  0.000000  0.011496        0.00000
ATOM      1  N   PRO H  14      52.628 -74.147  33.427  1.00 20.43           N
ATOM      2  CA  PRO H  14      53.440 -73.630  34.533  1.00 20.01           C
ATOM      3  C   PRO H  14      54.482 -72.584  34.124  1.00 20.76           C
ATOM      4  O   PRO H  14      55.025 -72.627  33.021  1.00 16.34           O
ATOM      5  CB  PRO H  14      54.055 -74.895  35.134  1.00 22.06           C
ATOM      6  CG  PRO H  14      54.084 -75.862  33.972  1.00 25.16           C
ATOM      7  CD  PRO H  14      52.770 -75.608  33.294  1.00 17.36           C
ATOM      8  HD2 PRO H  14      52.801 -75.872  32.361  1.00 17.36           H
ATOM      9  HD1 PRO H  14      52.048 -76.072  33.746  1.00 17.36           H
ATOM     10  HG2 PRO H  14      54.830 -75.664  33.385  1.00 25.16           H
ATOM     11  HG1 PRO H  14      54.147 -76.775  34.292  1.00 25.16           H
ATOM     12  HA  PRO H  14      52.877 -73.212  35.203  1.00 20.01           H
ATOM     13  HB1 PRO H  14      54.949 -74.711  35.461  1.00 22.06           H
ATOM     14  HB2 PRO H  14      53.499 -75.228  35.856  1.00 22.06           H
ATOM     15  N   SER H  15      54.727 -71.646  35.038  1.00 21.70           N
ATOM     16  CA  SER H  15      55.670 -70.537  34.874  1.00 25.33           C
ATOM     17  C   SER H  15      55.049 -69.401  34.057  1.00 24.78           C
ATOM     18  O   SER H  15      55.581 -68.291  34.023  1.00 27.51           O
ATOM     19  CB  SER H  15      56.982 -71.005  34.219  1.00 25.20           C
ATOM     20  OG  SER H  15      56.914 -70.938  32.802  1.00 28.91           O
ATOM     21  H   SER H  15      54.335 -71.634  35.803  1.00 21.70           H
ATOM     22  HA  SER H  15      55.899 -70.163  35.739  1.00 25.33           H
ATOM     23  HB1 SER H  15      57.705 -70.434  34.524  1.00 25.20           H
ATOM     24  HG  SER H  15      56.225 -70.518  32.567  1.00 28.91           H
ATOM     25  HB2 SER H  15      57.151 -71.924  34.481  1.00 25.20           H
ATOM     26  N   GLN H  16      53.918 -69.678  33.412  1.00 24.55           N
ATOM     27  CA  GLN H  16      53.224 -68.673  32.611  1.00 29.39           C
ATOM     28  C   GLN H  16      52.340 -67.778  33.475  1.00 28.13           C
ATOM     29  O   GLN H  16      52.234 -67.987  34.681  1.00 26.35           O
ATOM     30  CB  GLN H  16      52.371 -69.346  31.533  1.00 31.67           C
ATOM     31  CG  GLN H  16      53.196 -70.112  30.524  1.00 44.80           C
ATOM     32  CD  GLN H  16      54.379 -69.303  30.030  1.00 48.55           C
ATOM     33  OE1 GLN H  16      54.213 -68.269  29.386  1.00 52.45           O
ATOM     34  NE2 GLN H  16      55.584 -69.766  30.342  1.00 55.07           N
ATOM     35  H   GLN H  16      53.530 -70.445  33.423  1.00 24.55           H
ATOM     36  HG2 GLN H  16      53.533 -70.922  30.937  1.00 44.80           H
ATOM     37  HG1 GLN H  16      52.641 -70.335  29.761  1.00 44.80           H
ATOM     38 HE22 GLN H  16      55.661 -70.489  30.801  1.00 55.07           H
ATOM     39 HE21 GLN H  16      56.287 -69.343  30.085  1.00 55.07           H
ATOM     40  HA  GLN H  16      53.888 -68.112  32.179  1.00 29.39           H
ATOM     41  HB1 GLN H  16      51.871 -68.665  31.056  1.00 31.67           H
ATOM     42  HB2 GLN H  16      51.761 -69.970  31.957  1.00 31.67           H
TER
"""

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
