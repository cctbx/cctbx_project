from __future__ import division, print_function
import iotbx.pdb
import time, os, shutil
import iotbx.phil
from iotbx.data_manager import DataManager
from libtbx import easy_run
from libtbx.utils import Sorry, null_out
from libtbx.test_utils import Exception_expected
from iotbx.cli_parser import run_program
from mmtbx.programs import polder

pdb_str= """
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH A   1       1.000   2.000   3.000  1.00 10.00           O
HETATM    2  O   HOH A   2       4.000   5.000   6.000  1.00 20.00           O
END
"""

def write_model_and_data(prefix):
  pdb_inp = iotbx.pdb.input(lines = pdb_str, source_info=None)
  pdb_inp.write_pdb_file(file_name="%s.pdb"%prefix)
  xrs = pdb_inp.xray_structure_simple()
  fo = abs(xrs.structure_factors(d_min=1).f_calc())
  mtz_dataset = fo.as_mtz_dataset(column_root_label = "FOBS")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "%s.mtz"%prefix)

def exercise_00(prefix="tst_polder_4"):
  """
  Make sure EXACTLY THIS comand below runs (= does not crash):

  phenix.polder model.pdb data.mtz scattering_table=electron selectioon="resseq 2"
  """
  write_model_and_data(prefix)
  cmd = " ".join([
    "phenix.polder",
    "%s.mtz"%prefix,
    "%s.pdb"%prefix,
    "selection='resseq 1' ",
    "scattering_table=electron",
    "> %s.zlog"%prefix
    ])
  #print(cmd)
  assert not easy_run.call(cmd)
  assert os.path.isfile("tst_polder_4_polder_map_coeffs.mtz")

  os.remove('tst_polder_4.mtz')
  os.remove('tst_polder_4.pdb')
  os.remove('tst_polder_4.zlog')
  os.remove('tst_polder_4_polder_map_coeffs.mtz')

def exercise_01(prefix="tst_polder_4_twins"):
  """
  Make sure Sorry happens in case of twinning.
  """
  write_model_and_data(prefix)
  args = [
    "%s.mtz"%prefix,
    "%s.pdb"%prefix,
    "selection=resseq 1",
    "data_manager.fmodel.xray_data.twin_law=h,-h-k,-l"]
  try:
    run_program(program_class=polder.Program, args=args, logger=null_out())
  except Sorry as e:
    assert str(e) == "Twinning is not supported in Polder map calculation.", \
      str(e)
  else:
    raise Exception_expected

  os.remove('%s.mtz'%prefix)
  os.remove('%s.pdb'%prefix)

def exercise_02():
  """
  Program.result_message: every branch and the 0.7/0.8 cc13 boundaries.
  """
  params = iotbx.phil.parse(polder.Program.master_phil_str,
    process_includes=True).extract()
  program = polder.Program(data_manager=DataManager(), params=params,
    logger=null_out())
  shows   = "likely to show the omitted atoms."
  solvent = "very likely to show bulk-solvent or noise."
  maybe_atoms   = "more likely to show the omitted atoms than bulk solvent."
  maybe_solvent = "more likely to show bulk-solvent or noise instead of"
  # cc12, cc13, cc23, expected
  cases = [
    (0.40, 0.998, 0.38, shows),
    (0.40, 0.80,  0.38, shows),
    (0.90, 0.85,  0.38, shows),         # cc13 neither smallest nor cc23 top
    (0.30, 0.38,  0.95, solvent),
    (0.30, 0.69,  0.20, solvent),       # cc13 < 0.7
    (0.40, 0.85,  0.90, solvent),       # cc23 largest
    (0.90, 0.85,  0.88, solvent),       # cc13 smallest
    (0.40, 0.70,  0.48, maybe_atoms),   # cc23 < 0.7*cc13
    (0.40, 0.75,  0.52, maybe_atoms),
    (0.40, 0.75,  0.53, maybe_solvent), # cc23 >= 0.7*cc13
    (0.40, 0.79,  0.60, maybe_solvent)]
  for cc12, cc13, cc23, expected in cases:
    msg = program.result_message(cc12=cc12, cc13=cc13, cc23=cc23)
    assert expected in " ".join(msg.split()), (cc12, cc13, cc23, msg)

def exercise_03(prefix="tst_polder_4_errors"):
  """
  Make sure Sorry happens for invalid input, with the expected message.
  """
  write_model_and_data(prefix)
  pdb_fn, mtz_fn = "%s.pdb" % prefix, "%s.mtz" % prefix
  extra_files = ["%s_2.mtz" % prefix, "%s_3.mtz" % prefix,
                 "%s_cell.pdb" % prefix]
  shutil.copy(mtz_fn, extra_files[0])
  shutil.copy(mtz_fn, extra_files[1])
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str.replace(
    "10.000   10.000   10.000", "12.000   12.000   12.000"))
  pdb_inp.write_pdb_file(file_name=extra_files[2])
  base = [pdb_fn, mtz_fn, "selection=resseq 1"]
  # arguments, expected start of the Sorry message
  cases = [
    ([pdb_fn, mtz_fn],
     "Selection for atoms to be omitted is required."),
    (base + ["sphere_radius=2"],
     "Solvent exclusion radius must be at least 3 A."),
    (base + ["sphere_radius=None"],
     "Solvent exclusion radius must be at least 3 A."),
    (base + ["compute_box=True", "box_buffer=6"],
     "Box buffer out of range"),
    (base + ["compute_box=True", "box_buffer=-1"],
     "Box buffer out of range"),
    (base + ["resolution_factor=None"],
     "Resolution gridding factor is required."),
    (base + ["resolution_factor=-0.1"],
     "Use a positive value for the resolution gridding factor."),
    (base + ["resolution_factor=0.6"],
     "Use a value < 0.5 for the resolution gridding factor."),
    ([pdb_fn, mtz_fn, "selection=resseq 5"],
     "No atoms where selected."),
    ([pdb_fn, mtz_fn, "selection=all"],
     "More than half of total number of atoms selected."),
    (base + extra_files[:2],
     "Dont input more than 2 reflection files."),
    ([extra_files[2], mtz_fn, "selection=resseq 1"],
     "Crystal symmetry mismatch between different files."),
    # 2*5 A buffer spans the 10 A cell
    (base + ["compute_box=True", "box_buffer=5"],
     "The box is too big.")]
  for args, expected in cases:
    try:
      run_program(program_class=polder.Program, args=args, logger=null_out())
    except Sorry as e:
      assert str(e).startswith(expected), (args, str(e))
    else:
      raise Exception_expected
  # box mode does not use sphere_radius
  r = run_program(program_class=polder.Program,
    args=base + ["compute_box=True", "sphere_radius=None"], logger=null_out())
  os.remove(r.output_file)

  for fn in [pdb_fn, mtz_fn] + extra_files:
    os.remove(fn)

if (__name__ == "__main__"):
  t0 = time.time()
  exercise_01()
  exercise_00()
  exercise_02()
  exercise_03()
  print("Time: %6.4f" % (time.time()-t0))
