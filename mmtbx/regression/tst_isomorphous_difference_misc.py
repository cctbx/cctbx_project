
from __future__ import absolute_import, division, print_function
from iotbx.pdb import hierarchy
from libtbx import easy_run
import iotbx.pdb
from iotbx.file_reader import any_file
from mmtbx.regression.make_fake_anomalous_data import generate_calcium_inputs
import time
import os

def exercise_calcium_substitution():
  anom_mtz_file, pdb_file = generate_calcium_inputs(
    file_base = "tst_isomorphous_difference_misc_anom", anonymize = False)
  hoh_file = generate_calcium_inputs(
    file_base = "tst_isomorphous_difference_misc_hoh", anonymize = True)[1]
  mtz_file = "tst_isomorphous_difference_misc.mtz"
  args = [
    "phenix.fmodel",
    hoh_file,
    "high_resolution=1.5",
    "r_free_flags_fraction=0.1",
    "type=real",
    "output.file_name=" + mtz_file
    ]
  easy_run.fully_buffered(args).raise_if_errors()
  time.sleep(2)
  assert os.path.isfile(mtz_file)
  minus_mtz_file = "tst_isomorphous_difference_misc_anon_minus.mtz"
  args = [
    "phenix.fobs_minus_fobs_map",
    "f_obs_1_file=" + anom_mtz_file,
    "f_obs_2_file=" + mtz_file,
    pdb_file,
    "omit_selection=\"element CA\"",
    "multiscale=True",
    "output_file=" + minus_mtz_file,
  ]
  result = easy_run.fully_buffered(args).raise_if_errors()
  assert os.path.isfile(minus_mtz_file)
  assert ("1 atoms selected for removal" in result.stdout_lines)
  f = any_file(minus_mtz_file)
  coeffs = f.file_server.miller_arrays[0]
  fft_map = coeffs.fft_map(resolution_factor=0.25)
  minus_map = fft_map.apply_sigma_scaling().real_map_unpadded()
  pdb_in = iotbx.pdb.input(pdb_file)
  hierarchy = pdb_in.construct_hierarchy()
  xrs = pdb_in.xray_structure_simple()
  for i_seq, atom in enumerate(hierarchy.atoms()):
    if (atom.element == "CA"):
      site_frac = xrs.sites_frac()[i_seq]
      val = minus_map.eight_point_interpolation(site_frac)
      assert (val > 40)

def exercise_anomalous_isomorphous_difference_map():
  pdb_1 = """\
CRYST1   32.247   37.398   33.613  90.00  90.00  90.00 P 1
HETATM 1983  PG  AGS A 340       5.029   2.780   3.852  1.00 60.68           P
HETATM 1984  S1G AGS A 340       3.762   2.077   2.900  1.00 64.20           S
HETATM 1985  O2G AGS A 340       6.396   1.935   3.849  1.00 61.32           O
HETATM 1986  O3G AGS A 340       4.528   2.885   5.348  1.00 57.78           O
HETATM 1987  PB  AGS A 340       4.068   5.371   3.258  1.00 53.67           P
HETATM 1988  O1B AGS A 340       3.533   5.492   1.878  1.00 52.45           O
HETATM 1989  O2B AGS A 340       2.939   5.014   4.243  1.00 53.83           O
HETATM 1990  O3B AGS A 340       5.255   4.278   3.263  1.00 56.71           O
HETATM 1991  PA  AGS A 340       5.581   6.888   5.132  1.00 42.66           P
HETATM 1992  O1A AGS A 340       4.625   6.671   6.259  1.00 43.69           O
HETATM 1993  O2A AGS A 340       6.834   5.959   5.235  1.00 45.93           O
HETATM 1994  O3A AGS A 340       4.800   6.740   3.720  1.00 48.66           O
HETATM 1995  O5' AGS A 340       6.098   8.409   5.107  1.00 36.10           O
HETATM 1996  C5' AGS A 340       7.059   8.753   4.138  1.00 31.02           C
HETATM 1997  C4' AGS A 340       7.848   9.969   4.591  1.00 29.79           C
HETATM 1998  O4' AGS A 340       6.947  11.004   4.917  1.00 30.10           O
HETATM 1999  C3' AGS A 340       8.682   9.711   5.831  1.00 30.12           C
HETATM 2000  O3' AGS A 340       9.998   9.322   5.503  1.00 29.10           O
HETATM 2001  C2' AGS A 340       8.619  11.034   6.578  1.00 30.77           C
HETATM 2002  O2' AGS A 340       9.677  11.886   6.190  1.00 29.37           O
HETATM 2003  C1' AGS A 340       7.319  11.667   6.111  1.00 30.26           C
HETATM 2004  N9  AGS A 340       6.261  11.498   7.129  1.00 27.43           N
HETATM 2005  C8  AGS A 340       5.454  10.400   7.302  1.00 27.85           C
HETATM 2006  N7  AGS A 340       4.609  10.631   8.330  1.00 29.03           N
HETATM 2007  C5  AGS A 340       4.862  11.863   8.817  1.00 26.18           C
HETATM 2008  C6  AGS A 340       4.295  12.588   9.859  1.00 26.06           C
HETATM 2009  N6  AGS A 340       3.275  12.079  10.559  1.00 28.35           N
HETATM 2010  N1  AGS A 340       4.767  13.859  10.138  1.00 29.55           N
HETATM 2011  C2  AGS A 340       5.794  14.401   9.385  1.00 29.19           C
HETATM 2012  N3  AGS A 340       6.356  13.672   8.351  1.00 28.57           N
HETATM 2013  C4  AGS A 340       5.896  12.422   8.071  1.00 25.89           C
HETATM 2014 MN    MN A 341       2.683   3.335   5.745  1.00 47.43          MN
"""
  pdb_2 = """\
CRYST1   32.247   37.398   33.613  90.00  90.00  90.00 P 1
HETATM 1983  PG  AGS A 340       5.029   2.780   3.852  1.00 60.68           P
HETATM 1984  S1G AGS A 340       3.762   2.077   2.900  1.00 64.20           S
HETATM 1985  O2G AGS A 340       6.396   1.935   3.849  1.00 61.32           O
HETATM 1986  O3G AGS A 340       4.528   2.885   5.348  1.00 57.78           O
HETATM 1987  PB  AGS A 340       4.068   5.371   3.258  1.00 53.67           P
HETATM 1988  O1B AGS A 340       3.533   5.492   1.878  1.00 52.45           O
HETATM 1989  O2B AGS A 340       2.939   5.014   4.243  1.00 53.83           O
HETATM 1990  O3B AGS A 340       5.255   4.278   3.263  1.00 56.71           O
HETATM 1991  PA  AGS A 340       5.581   6.888   5.132  1.00 42.66           P
HETATM 1992  O1A AGS A 340       4.625   6.671   6.259  1.00 43.69           O
HETATM 1993  O2A AGS A 340       6.834   5.959   5.235  1.00 45.93           O
HETATM 1994  O3A AGS A 340       4.800   6.740   3.720  1.00 48.66           O
HETATM 1995  O5' AGS A 340       6.098   8.409   5.107  1.00 36.10           O
HETATM 1996  C5' AGS A 340       7.059   8.753   4.138  1.00 31.02           C
HETATM 1997  C4' AGS A 340       7.848   9.969   4.591  1.00 29.79           C
HETATM 1998  O4' AGS A 340       6.947  11.004   4.917  1.00 30.10           O
HETATM 1999  C3' AGS A 340       8.682   9.711   5.831  1.00 30.12           C
HETATM 2000  O3' AGS A 340       9.998   9.322   5.503  1.00 29.10           O
HETATM 2001  C2' AGS A 340       8.619  11.034   6.578  1.00 30.77           C
HETATM 2002  O2' AGS A 340       9.677  11.886   6.190  1.00 29.37           O
HETATM 2003  C1' AGS A 340       7.319  11.667   6.111  1.00 30.26           C
HETATM 2004  N9  AGS A 340       6.261  11.498   7.129  1.00 27.43           N
HETATM 2005  C8  AGS A 340       5.454  10.400   7.302  1.00 27.85           C
HETATM 2006  N7  AGS A 340       4.609  10.631   8.330  1.00 29.03           N
HETATM 2007  C5  AGS A 340       4.862  11.863   8.817  1.00 26.18           C
HETATM 2008  C6  AGS A 340       4.295  12.588   9.859  1.00 26.06           C
HETATM 2009  N6  AGS A 340       3.275  12.079  10.559  1.00 28.35           N
HETATM 2010  N1  AGS A 340       4.767  13.859  10.138  1.00 29.55           N
HETATM 2011  C2  AGS A 340       5.794  14.401   9.385  1.00 29.19           C
HETATM 2012  N3  AGS A 340       6.356  13.672   8.351  1.00 28.57           N
HETATM 2013  C4  AGS A 340       5.896  12.422   8.071  1.00 25.89           C
HETATM 2014 MN    MN A 341       1.319   3.281   7.076  0.50 47.43          MN
HETATM   60  C   ACT     1       2.759   8.395   7.862  1.00 18.56           C
HETATM   61  O   ACT     1       3.812   8.507   7.193  1.00 19.26           O
HETATM   62  OXT ACT     1       2.218   9.388   8.383  1.00 18.36           O
HETATM   63  CH3 ACT     1       2.127   7.078   8.030  1.00 18.08           C
"""
  with open("tst_anom_iso_diff_1.pdb", "w") as f:
    f.write(pdb_1)
  with open("tst_anom_iso_diff_2.pdb", "w") as f:
    f.write(pdb_2)
  base_args = [
    "phenix.fmodel",
    "high_resolution=1.5",
    "type=real",
    "label=F",
    "r_free_flags_fraction=0.1",
    "wavelength=1.77",
    "add_random_error_to_amplitudes_percent=5",
  ]
  args_1 = base_args + [ "tst_anom_iso_diff_1.pdb",
    "output.file_name=tst_anom_iso_diff_1.mtz" ]
  args_2 = base_args + [ "tst_anom_iso_diff_2.pdb",
    "output.file_name=tst_anom_iso_diff_2.mtz" ]
  print(" ".join(args_1))
  rc = easy_run.fully_buffered(" ".join(args_1)).raise_if_errors().return_code
  assert (rc == 0)
  print(" ".join(args_2))
  rc = easy_run.fully_buffered(" ".join(args_2)).raise_if_errors().return_code
  assert (rc == 0)
  base_args = [
    "phenix.fobs_minus_fobs_map",
    "f_obs_1_file=tst_anom_iso_diff_2.mtz",
    "f_obs_2_file=tst_anom_iso_diff_1.mtz",
    "tst_anom_iso_diff_2.pdb",
    "omit_selection=\"element MN\"",
  ]
  args_1 = base_args + ["output_file=tst_anom_iso_diff_map_coeffs.mtz",
    "anomalous=True",]
  args_2 = base_args + ["output_file=tst_anom_iso_diff_map_coeffs_control.mtz"]
  print(" " .join(args_1))
  rc = easy_run.fully_buffered(" " .join(args_1)).raise_if_errors().return_code
  assert (rc == 0)
  print(" " .join(args_2))
  rc = easy_run.fully_buffered(" " .join(args_2)).raise_if_errors().return_code
  assert (rc == 0)
  # The input structures differ in the position of the MN ion and in the
  # presence of an acetate ion in the second model.  I am calculating two
  # isomorphous difference maps: one using the anomalous differences of each
  # dataset, and a control using the merged amplitudes. The MN should be
  # prominent in both, and by far the strongest feature in the anomalous map,
  # but the ACT should only have density in the control map.
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_2)
  xrs = pdb_in.xray_structure_simple()
  mtz_1 = any_file("tst_anom_iso_diff_map_coeffs.mtz")
  map_1 = mtz_1.file_server.miller_arrays[0].fft_map(
    resolution_factor=0.25).apply_sigma_scaling().real_map_unpadded()
  mtz_2 = any_file("tst_anom_iso_diff_map_coeffs_control.mtz")
  map_2 = mtz_2.file_server.miller_arrays[0].fft_map(
    resolution_factor=0.25).apply_sigma_scaling().real_map_unpadded()
  sites_frac = xrs.sites_frac()
  anom_max = mn_anom = 0
  hierarchy = pdb_in.construct_hierarchy()
  for i_seq, atom in enumerate(hierarchy.atoms()):
    site_frac = sites_frac[i_seq]
    anom_diff = map_1.eight_point_interpolation(site_frac)
    fobs_diff = map_2.eight_point_interpolation(site_frac)
    labels = atom.fetch_labels()
    if (labels.resname.strip() == "MN"):
      assert (anom_diff > 10) and (fobs_diff > 5)
      mn_anom = anom_diff
    elif (labels.resname.strip() == "ACT"):
      assert (anom_diff < 3) and (fobs_diff > 5)
    if (anom_diff > anom_max):
      anom_max = anom_diff
  assert (anom_max == mn_anom)

if (__name__ == "__main__"):
  exercise_anomalous_isomorphous_difference_map()
  exercise_calcium_substitution()
  print("OK")
