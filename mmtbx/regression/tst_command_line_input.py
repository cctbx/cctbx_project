
# TODO make this much more comprehensive

from __future__ import division
import mmtbx.command_line
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.utils import null_out, Sorry
from libtbx import easy_run
import libtbx.load_env
from cStringIO import StringIO
import os.path

def exercise_simple () :
  import iotbx.pdb
  pdb_str="""
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
END
"""
  pdb_in = iotbx.pdb.input(source_info=None,lines=pdb_str)
  hierarchy = pdb_in.construct_hierarchy()
  xrs = pdb_in.xray_structure_simple()
  xrs.scattering_type_registry(
    d_min=1.5,
    table="n_gaussian")
  file_base = "tmp_mmtbx_cmdline"
  open(file_base+".pdb", "w").write(
    hierarchy.as_pdb_string(crystal_symmetry=xrs))
  fc = abs(xrs.structure_factors(d_min=1.5).f_calc())
  flags = fc.generate_r_free_flags()
  mtz = fc.as_mtz_dataset(column_root_label="F")
  mtz.add_miller_array(flags, column_root_label="FreeR_flag")
  mtz.mtz_object().write(file_base+".mtz")
  open(file_base+".fa", "w").write(">Tyr\nY\n")
  cmdline = mmtbx.command_line.load_model_and_data(
    update_f_part1_for=None,
    args=[ file_base + ext for ext in [".pdb",".mtz",".fa"] ],
    master_phil=mmtbx.command_line.generate_master_phil_with_inputs(""),
    out=StringIO(),
    create_log_buffer=True)
  assert (cmdline.params.input.xray_data.file_name is not None)
  assert (cmdline.sequence is not None)
  r_factor = cmdline.fmodel.r_work()
  assert (r_factor < 0.001)
  log = cmdline.start_log_file("tst_mmtbx_cmdline.log")

def exercise_combine_symmetry () :
  """
  Test the extraction of symmetry from both a PDB file and an MTZ file.
  """
  from mmtbx.regression import model_1yjp
  import mmtbx.command_line
  import iotbx.pdb.hierarchy
  from cctbx import sgtbx
  from cctbx import uctbx
  # 1yjp, as usual
  pdb_in = iotbx.pdb.hierarchy.input(pdb_string=model_1yjp)
  xrs = pdb_in.input.xray_structure_simple()
  f = open("tst_combine_symmetry.pdb", "w")
  f.write(pdb_in.hierarchy.as_pdb_string(crystal_symmetry=xrs))
  f.close()
  f_calc = abs(xrs.structure_factors(d_min=1.5).f_calc())
  # Make up slightly more exact unit cell, but set SG to P2
  f_calc = f_calc.customized_copy(
    crystal_symmetry=f_calc.crystal_symmetry().customized_copy(
      space_group_info=sgtbx.space_group_info("P2"),
      unit_cell=uctbx.unit_cell((21.9371, 4.8659, 23.4774, 90.0, 107.0832,
        90.00))))
  flags = f_calc.generate_r_free_flags()
  mtz = f_calc.as_mtz_dataset(column_root_label="F")
  mtz.add_miller_array(flags, column_root_label="FreeR_flag")
  mtz.mtz_object().write("tst_combine_symmetry.mtz")
  cmdline = mmtbx.command_line.load_model_and_data(
    args=["tst_combine_symmetry.pdb", "tst_combine_symmetry.mtz"],
    master_phil=mmtbx.command_line.generate_master_phil_with_inputs(""),
    process_pdb_file=False,
    create_fmodel=True,
    out=null_out())
  symm = cmdline.xray_structure.crystal_symmetry()
  assert (approx_equal(symm.unit_cell().parameters(),
          (21.9371, 4.8659, 23.4774, 90.0, 107.0832, 90.0)))
  assert (str(symm.space_group_info()) == "P 1 21 1")
  # Part 2: incompatible space groups
  f_calc_2 = f_calc.customized_copy(
    crystal_symmetry=f_calc.crystal_symmetry().customized_copy(
      space_group_info=sgtbx.space_group_info("P1")))
  flags_2 = f_calc_2.generate_r_free_flags()
  mtz = f_calc_2.as_mtz_dataset(column_root_label="F")
  mtz.add_miller_array(flags_2, column_root_label="FreeR_flag")
  mtz.mtz_object().write("tst_combine_symmetry_2.mtz")
  try :
    cmdline = mmtbx.command_line.load_model_and_data(
      args=["tst_combine_symmetry.pdb", "tst_combine_symmetry_2.mtz"],
      master_phil=mmtbx.command_line.generic_simple_input_phil(),
      process_pdb_file=False,
      create_fmodel=True,
      out=null_out())
  except Sorry, s :
    assert ("Incompatible space groups" in str(s))
  else :
    raise Exception_expected
  # Part 3: unit cell mismatch
  f_calc_3 = f_calc.customized_copy(
    crystal_symmetry=f_calc.crystal_symmetry().customized_copy(
      unit_cell=uctbx.unit_cell((21.9371, 4.8659, 23.4774, 90.0, 104.0123,
        90.00))))
  flags_3 = f_calc_3.generate_r_free_flags()
  mtz = f_calc_3.as_mtz_dataset(column_root_label="F")
  mtz.add_miller_array(flags_3, column_root_label="FreeR_flag")
  mtz.mtz_object().write("tst_combine_symmetry_3.mtz")
  try :
    cmdline = mmtbx.command_line.load_model_and_data(
      args=["tst_combine_symmetry.pdb", "tst_combine_symmetry_3.mtz"],
      master_phil=mmtbx.command_line.generic_simple_input_phil(),
      process_pdb_file=False,
      create_fmodel=True,
      out=null_out())
  except Sorry, s :
    assert ("Unit cell mismatch" in str(s))
  else :
    raise Exception_expected

def exercise_example () :
  from mmtbx.regression import make_fake_anomalous_data
  pdb_file, mtz_file = make_fake_anomalous_data.generate_cd_cl_inputs(
    file_base="tst_cmdline_example")
  script = os.path.join(libtbx.env.dist_path("mmtbx"), "examples",
    "simple_command_line_cc.py")
  assert os.path.isfile(script)
  args = [
    "mmtbx.python",
    "\"%s\"" % script,
    pdb_file,
    mtz_file,
    "skip_twin_detection=True",
  ]
  result = easy_run.fully_buffered(" ".join(args)).raise_if_errors()
  assert ("""CC(obs-calc): 0.999""" in result.stdout_lines)

def exercise_cns_input () :
  from mmtbx.regression import make_fake_anomalous_data
  pdb_file, mtz_file = make_fake_anomalous_data.generate_cd_cl_inputs(
    file_base="tst_cmdline_cns")
  from iotbx.file_reader import any_file
  mtz_in = any_file("tst_cmdline_cns.mtz")
  f_obs = mtz_in.file_server.miller_arrays[0].average_bijvoet_mates()
  flags = mtz_in.file_server.miller_arrays[1].average_bijvoet_mates()
  f = open("tst_cmdline_cns.hkl", "w")
  out = StringIO()
  f_obs.export_as_cns_hkl(
    file_object=out,
    r_free_flags=flags)
  # get rid of embedded symmetry
  for line in out.getvalue().splitlines() :
    if (not "{" in line) :
      f.write("%s\n" % line)
  f.close()
  cmdline = mmtbx.command_line.load_model_and_data(
    args=["tst_cmdline_cns.pdb", "tst_cmdline_cns.hkl"],
    master_phil=mmtbx.command_line.generic_simple_input_phil(),
    process_pdb_file=False,
    create_fmodel=True,
    out=null_out())
  cmdline.crystal_symmetry.show_summary()

if (__name__ == "__main__") :
  exercise_simple()
  exercise_combine_symmetry()
  exercise_example()
  exercise_cns_input()
  print "OK"
