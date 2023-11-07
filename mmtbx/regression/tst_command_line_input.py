
# TODO make this much more comprehensive

from __future__ import absolute_import, division, print_function
from mmtbx.regression import model_1yjp
import mmtbx.command_line
from iotbx.scalepack import no_merge_original_index
import iotbx.pdb
from cctbx import sgtbx
from cctbx import uctbx
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.utils import null_out, Sorry
from libtbx import easy_run
import libtbx.load_env
from six.moves import cStringIO as StringIO
import random
import os.path

def exercise_simple():
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
TER
"""
  pdb_in = iotbx.pdb.input(source_info=None,lines=pdb_str)
  hierarchy = pdb_in.construct_hierarchy()
  xrs = pdb_in.xray_structure_simple()
  xrs.scattering_type_registry(
    d_min=1.5,
    table="n_gaussian")
  xrs.set_inelastic_form_factors(
    photon=1.54,
    table="sasaki")
  file_base = "tmp_mmtbx_cmdline"
  with open(file_base+".pdb", "w") as f:
    f.write(hierarchy.as_pdb_string(crystal_symmetry=xrs))
  fc = abs(xrs.structure_factors(d_min=1.5).f_calc())
  flags = fc.generate_r_free_flags()
  mtz = fc.as_mtz_dataset(column_root_label="F")
  mtz.add_miller_array(flags, column_root_label="FreeR_flag")
  mtz.mtz_object().write(file_base+".mtz")
  with open(file_base+".fa", "w") as f:
    f.write(">Tyr\nY\n")
  base_args = [ file_base + ext for ext in [".pdb",".mtz",".fa"] ]
  cmdline = mmtbx.command_line.load_model_and_data(
    args=base_args+["wavelength=1.54"],
    master_phil=mmtbx.command_line.generate_master_phil_with_inputs(""),
    out=StringIO(),
    create_log_buffer=True)
  assert (cmdline.params.input.xray_data.file_name is not None)
  assert (cmdline.sequence is not None)
  r_factor = cmdline.fmodel.r_work()
  assert (r_factor < 0.002)
  cmdline.save_data_mtz("tmp_mmtbx_cmdline_data.mtz")
  assert os.path.isfile("tmp_mmtbx_cmdline_data.mtz")
  model = cmdline.create_model_manager()
  # energy input
  cmdline = mmtbx.command_line.load_model_and_data(
    args=base_args+["energy=8050"],
    master_phil=mmtbx.command_line.generate_master_phil_with_inputs(""),
    out=StringIO(),
    create_log_buffer=True)
  assert approx_equal(cmdline.params.input.wavelength, 1.54018, eps=0.0001)
  # UNMERGED DATA INPUT
  log = cmdline.start_log_file("tst_mmtbx_cmdline.log")
  log.close()
  fc2 = xrs.structure_factors(d_min=1.3).f_calc().generate_bijvoet_mates()
  fc2 = fc2.randomize_amplitude_and_phase(amplitude_error=0.01,
    phase_error_deg=5, random_seed=12345).customized_copy(
      sigmas=flex.random_double(fc2.size(), 10))
  i_obs = abs(fc2).f_as_f_sq()
  i_obs = i_obs.expand_to_p1().customized_copy(
    crystal_symmetry=fc2).set_observation_type_xray_intensity()
  with open(file_base + ".sca", "w") as f:
    no_merge_original_index.writer(i_obs, file_object=f)
  master_phil = mmtbx.command_line.generate_master_phil_with_inputs(
    phil_string="",
    enable_unmerged_data=True)
  cmdline = mmtbx.command_line.load_model_and_data(
    args=[ file_base + ext for ext in [".pdb",".mtz",".fa",] ] +
      ["unmerged_data.file_name=%s.sca" % file_base ],
    master_phil=master_phil,
    out=StringIO(),
    create_log_buffer=True)
  assert (cmdline.unmerged_i_obs is not None)
  # test with unknown scatterers
  pdb_in = iotbx.pdb.input(source_info=None,lines=pdb_str+"""\
ATOM     59  UNK UNL A   7       0.000   0.000   0.000  1.00 20.00           X
""")
  hierarchy = pdb_in.construct_hierarchy()
  file_base = "tmp_mmtbx_cmdline"
  with open(file_base+".pdb", "w") as f:
    f.write(hierarchy.as_pdb_string(crystal_symmetry=xrs))
  try :
    cmdline = mmtbx.command_line.load_model_and_data(
      args=[ file_base + ext for ext in [".pdb",".mtz",".fa",] ],
      master_phil=master_phil,
      out=StringIO(),
      process_pdb_file=False,
      create_log_buffer=True)
  except Sorry :
    pass
  else :
    raise Exception_expected
  cmdline = mmtbx.command_line.load_model_and_data(
    args=[ file_base + ext for ext in [".pdb",".mtz",".fa",] ],
    master_phil=master_phil,
    out=StringIO(),
    process_pdb_file=False,
    create_log_buffer=True,
    remove_unknown_scatterers=True)

def exercise_combine_symmetry():
  """
  Test the extraction of symmetry from both a PDB file and an MTZ file.
  """
  from mmtbx.regression import model_1yjp
  import mmtbx.command_line
  import iotbx.pdb
  from cctbx import sgtbx
  from cctbx import uctbx
  # 1yjp, as usual
  pdb_in = iotbx.pdb.input(source_info=None, lines=model_1yjp)
  hierarchy = pdb_in.construct_hierarchy()
  xrs = pdb_in.xray_structure_simple()
  f = open("tst_combine_symmetry.pdb", "w")
  f.write(hierarchy.as_pdb_string(crystal_symmetry=xrs))
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
          (21.9371, 4.8659, 23.4774, 90.0, 107.083, 90.0), eps=0.001))
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
  except Sorry as s :
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
  except Sorry as s :
    assert ("Unit cell mismatch" in str(s))
  else :
    raise Exception_expected

def exercise_example():
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

def exercise_cns_input():
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
  for line in out.getvalue().splitlines():
    if (not "{" in line):
      f.write("%s\n" % line)
  f.close()
  cmdline = mmtbx.command_line.load_model_and_data(
    args=["tst_cmdline_cns.pdb", "tst_cmdline_cns.hkl"],
    master_phil=mmtbx.command_line.generic_simple_input_phil(),
    process_pdb_file=False,
    create_fmodel=True,
    out=null_out())
  out = StringIO()
  cmdline.crystal_symmetry.show_summary(f=out)
  assert (out.getvalue() == """\
Unit cell: (21.362, 23.436, 23.594, 90, 90, 90)
Space group: P 1 (No. 1)
"""), out.getvalue()

# XXX this is essentially identical to tst_cc_star_space_group.py
def exercise_load_unmerged():
  flex.set_random_seed(123456)
  random.seed(123456)
  base = "tst_load_unmerged"
  pdb_in = iotbx.pdb.input(source_info=None, lines=model_1yjp)
  xrs = pdb_in.xray_structure_simple()
  xrs.set_inelastic_form_factors(
    photon=1.54,
    table="sasaki")
  fc = abs(xrs.structure_factors(d_min=1.5).f_calc()).average_bijvoet_mates()
  fc.set_observation_type_xray_amplitude()
  flags = fc.generate_r_free_flags()
  mtz = fc.as_mtz_dataset(column_root_label="F")
  mtz.add_miller_array(flags, column_root_label="FreeR_flag")
  mtz.mtz_object().write(base + ".mtz")
  xrs_p1 = xrs.expand_to_p1()
  xrs_p1.shake_sites_in_place(rms_difference=0.1)
  fc_p1 = xrs_p1.structure_factors(d_min=1.4).f_calc()
  fc_p1_extra = fc_p1.randomize_amplitude_and_phase(amplitude_error=1.0,
    phase_error_deg=0,
    random_seed=123456)
  fc_p1 = abs(fc_p1.concatenate(other=fc_p1_extra)).sort(
    by_value="packed_indices")
  fc_p1.set_observation_type_xray_amplitude()
  sg_p2 = sgtbx.space_group_info("P2")
  ic = fc_p1.f_as_f_sq().customized_copy(
    space_group_info=sg_p2,
    sigmas=flex.double(fc_p1.size(), 10.0))
  ic.export_as_scalepack_unmerged(file_name=base + ".sca")
  with open(base + ".pdb", "w") as f:
    f.write(model_1yjp)
  args = [
    base + ".mtz",
    base + ".pdb",
    "unmerged_data=%s.sca" % base,
  ]
  master_phil = mmtbx.command_line.generate_master_phil_with_inputs(
    phil_string="",
    enable_unmerged_data=True)
  cmdline = mmtbx.command_line.load_model_and_data(
    args=[ base + ext for ext in [".pdb",".mtz",] ] +
      ["unmerged_data.file_name=%s.sca" % base ],
    master_phil=master_phil,
    out=StringIO(),
    create_fmodel=False,
    process_pdb_file=False,
    create_log_buffer=True)
  # now with .sca in P1 (raises Sorry)
  ic2 = fc_p1.f_as_f_sq().customized_copy(
    sigmas=flex.double(fc_p1.size(), 10.0))
  ic2.export_as_scalepack_unmerged(file_name=base + "_p1.sca")
  args = [
    base + ".mtz",
    base + ".pdb",
    "unmerged_data=%s_p1.sca" % base,
  ]
  try :
    cmdline = mmtbx.command_line.load_model_and_data(
      args=[ base + ext for ext in [".pdb",".mtz",] ] +
      ["unmerged_data.file_name=%s_p1.sca" % base ],
      master_phil=master_phil,
      out=StringIO(),
      create_fmodel=False,
      process_pdb_file=False,
      create_log_buffer=True)
  except Sorry as s :
    assert (str(s) == "Incompatible space groups in merged and unmerged data:P 1 21 1 versus P 1"), s
  else :
    raise Exception_expected
  # XXX
  f = open(base + ".cif", "w")
  ic.as_cif_simple(array_type="meas",
    out=f)
  f.close()
  args = [
    base + ".mtz",
    base + ".pdb",
    "unmerged_data=%s.cif" % base,
  ]
  cmdline = mmtbx.command_line.load_model_and_data(
    args=[ base + ext for ext in [".pdb",".mtz",] ] +
      ["unmerged_data.file_name=%s.cif" % base ],
    master_phil=master_phil,
    out=StringIO(),
    create_fmodel=False,
    process_pdb_file=False,
    create_log_buffer=True)
  # bad unit cell
  uc2 = uctbx.unit_cell((23,6.5,23.5,90,108,90))
  ic3 = ic.customized_copy(unit_cell=uc2)
  f = open(base + "_new_uc.cif", "w")
  ic3.as_cif_simple(array_type="meas",
    out=f)
  f.close()
  args = [
    base + ".mtz",
    base + ".pdb",
    "unmerged_data=%s_new_uc.cif" % base,
  ]
  try :
    cmdline = mmtbx.command_line.load_model_and_data(
      args=[ base + ext for ext in [".pdb",".mtz",] ] +
        ["unmerged_data.file_name=%s_new_uc.cif" % base ],
      master_phil=master_phil,
      out=StringIO(),
      create_fmodel=False,
      process_pdb_file=False,
      create_log_buffer=True)
  except Sorry as s :
    assert ("Incompatible symmetry definitions" in str(s)), s
  else :
    raise Exception_expected

if (__name__ == "__main__"):
  exercise_simple()
  exercise_cns_input()
  exercise_load_unmerged()
  exercise_combine_symmetry()
  exercise_example()
  print("OK")
