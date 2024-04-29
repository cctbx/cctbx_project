
from __future__ import absolute_import, division, print_function
import os
from mmtbx.regression import make_fake_anomalous_data
from mmtbx.programs import fmodel
from iotbx.cli_parser import run_program
from iotbx import file_reader
from cctbx import miller
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.utils import null_out, Sorry
import iotbx.pdb


def exercise():
  if (os.path.isfile("tst_fmodel_anomalous.mtz")):
    os.remove("tst_fmodel_anomalous.mtz")
  pdb_file = make_fake_anomalous_data.write_pdb_input_cd_cl(
    file_base="tst_fmodel_anomalous")
  # phenix.fmodel (with wavelength)
  args = [
    pdb_file,
    "high_resolution=1.0",
    "wavelength=1.116",
    "label=F",
    "type=real",
    "output.file_name=tst_fmodel_anomalous.mtz",
    "r_free_flags_fraction=0.1",
  ]
  run_program(program_class=fmodel.Program, args=args, logger=null_out())
  assert os.path.isfile("tst_fmodel_anomalous.mtz")
  mtz_in = file_reader.any_file("tst_fmodel_anomalous.mtz")
  array = mtz_in.file_server.miller_arrays[0]
  assert (array.anomalous_flag())
  anom_diffs = array.anomalous_differences()
  assert approx_equal(flex.max(anom_diffs.data()), 5.72, eps=0.01)
  print("OK")
  os.remove('tst_fmodel_anomalous.mtz')
  os.remove('tst_fmodel_anomalous.pdb')

def exercise_intensity_output():
  if (os.path.isfile("tst_fmodel_anomalous.mtz")):
    os.remove("tst_fmodel_anomalous.mtz")
  pdb_file = make_fake_anomalous_data.write_pdb_input_cd_cl(
    file_base="tst_fmodel_anomalous")
  # phenix.fmodel (with wavelength)
  args = [
    pdb_file,
    "high_resolution=1.0",
    "wavelength=1.116",
    "obs_type=intensities",
    "type=real",
    "output.file_name=tst_fmodel_intensity.mtz",
    "r_free_flags_fraction=0.1",
  ]
  args2 = args + ["label=Imodel"]
  run_program(program_class=fmodel.Program, args=args2)
  assert os.path.isfile("tst_fmodel_intensity.mtz")
  mtz_in = file_reader.any_file("tst_fmodel_intensity.mtz")
  assert mtz_in.file_server.miller_arrays[0].is_xray_intensity_array()
  try :
    run_program(program_class=fmodel.Program, args=args)
  except Sorry :
    pass
  else :
    raise Exception_expected
  try :
    run_program(program_class=fmodel.Program, args=args+["format=cns"],
      logger=null_out())
  except Sorry :
    pass
  else :
    raise Exception_expected

  os.remove('tst_fmodel_intensity.mtz')

def exercise_selection_consistency():
  """
  Test that the atom selections for anomalous scatterers actually correspond
  to the intended xray scatterers.  This will be dependent on any re-ordering
  done by pdb_inp.construct_hierarchy(), which will move the SE atom in the
  structure below.
  """
  pdb_str = """\
ATOM    920  N   LYS A 123      -3.350   9.199  20.988  1.00 21.23           N
ATOM    921  CA  LYS A 123      -2.285   8.581  20.226  1.00 21.43           C
ATOM    922  C   LYS A 123      -1.116   9.544  20.127  1.00 21.24           C
ATOM    923  O   LYS A 123      -1.316  10.753  19.966  1.00 22.19           O
ATOM    924  CB  LYS A 123      -2.800   8.221  18.845  1.00 21.81           C
ATOM    925  CG  LYS A 123      -3.908   7.201  18.856  1.00 21.27           C
ATOM    926  CD  LYS A 123      -4.383   6.912  17.424  1.00 23.64           C
ATOM    927  CE  LYS A 123      -5.534   5.898  17.333  1.00 26.54           C
ATOM    928  NZ  LYS A 123      -5.023   4.519  17.519  1.00 30.98           N
HETATM  929  N   MSE A 124       0.105   8.997  20.163  1.00 22.33           N
HETATM  930  CA AMSE A 124       1.288   9.830  20.208  0.50 19.27           C
HETATM  932  C   MSE A 124       2.339   9.438  19.174  1.00 21.61           C
HETATM  933  O   MSE A 124       2.500   8.250  18.848  1.00 22.32           O
HETATM  934  CB AMSE A 124       1.906   9.779  21.634  0.50 22.80           C
HETATM  936  CG AMSE A 124       0.966  10.314  22.733  0.50 18.31           C
HETATM  938 SE   MSE A 124       1.972  10.293  24.432  0.37 22.53          SE
HETATM  940  CE AMSE A 124       1.818   8.557  24.764  0.50 27.54           C
ATOM    942  N   LEU A 125       3.059  10.448  18.722  1.00 21.43           N
ATOM    943  CA  LEU A 125       4.237  10.292  17.862  1.00 21.15           C
ATOM    944  C   LEU A 125       5.364  10.848  18.730  1.00 22.73           C
ATOM    945  O   LEU A 125       5.445  12.049  18.960  1.00 23.20           O
ATOM    946  CB  LEU A 125       4.108  11.070  16.559  1.00 22.37           C
ATOM    947  CG  LEU A 125       5.197  11.005  15.504  1.00 22.27           C
ATOM    948  CD1 LEU A 125       4.812  11.903  14.306  1.00 22.88           C
ATOM    949  CD2 LEU A 125       6.533  11.408  15.996  1.00 26.02           C
"""
  with open("tst_fmodel_misc.pdb", "w") as f:
    f.write(pdb_str)
  with open("tst_fmodel_misc.eff", "w") as f:
    f.write("""\
high_resolution = 1.0
output.file_name = tst_fmodel_misc.mtz
generate_fake_p1_symmetry = True
anomalous_scatterers {
  group {
    selection = element SE
    f_prime = -8.0
    f_double_prime = 4.5
  }
}
""")
  run_program(program_class=fmodel.Program, args=["tst_fmodel_misc.pdb","tst_fmodel_misc.eff"], logger=null_out())
  mtz_in = file_reader.any_file("tst_fmodel_misc.mtz")
  f_model = mtz_in.file_server.miller_arrays[0]
  dano = abs(f_model).anomalous_differences()
  f_model = f_model.average_bijvoet_mates()
  dano, f_model = dano.common_sets(other=f_model)
  map_coeffs = dano.phase_transfer(phase_source=f_model)
  map_coeffs = miller.array(
    miller_set=map_coeffs,
    data=map_coeffs.data()/(2j))
  map_coeffs.as_mtz_dataset(column_root_label="ANOM").mtz_object().write("anom.mtz")
  fft_map = map_coeffs.fft_map(resolution_factor=0.25).apply_sigma_scaling()
  real_map = fft_map.real_map_unpadded()
  hierarchy = iotbx.pdb.input("tst_fmodel_misc.pdb").construct_hierarchy()
  for atom in hierarchy.atoms():
    if (atom.element == "SE"):
      site = f_model.unit_cell().fractionalize(site_cart=atom.xyz)
      map_val = real_map.eight_point_interpolation(site)
      assert (map_val > 100)
    elif (atom.element.strip() == "O"):
      site = f_model.unit_cell().fractionalize(site_cart=atom.xyz)
      map_val = real_map.eight_point_interpolation(site)
      assert (map_val < 5)

  os.remove('tst_fmodel_misc.mtz')
  os.remove('tst_fmodel_misc.pdb')
  os.remove('tst_fmodel_misc.eff')
  os.remove('anom.mtz')

if (__name__ == "__main__"):
  exercise_intensity_output()
  exercise_selection_consistency()
  exercise()
