from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from mmtbx.bulk_solvent import mosaic
from mmtbx.bulk_solvent import tst_mosaic as tm
from libtbx.test_utils import approx_equal
from six.moves import cStringIO as StringIO
import mmtbx.f_model
import iotbx.pdb
import random, re, time

"""
Programs that use the mosaic bulk-solvent model. Each one is run on a model
with three cavities (one empty, two partially filled) and has to
- ask for the mosaic model (mosaic.number_of_runs), and
- report the R-factor of the mosaic model, not of the flat one.

R-work of this model with default scaling (resolution bins) is about 0.056
with the flat model and 0.01-0.03 with the mosaic one, depending on details
such as the choice of free reflections.
"""

PDB, MTZ = "tst_mosaic_programs.pdb", "tst_mosaic_programs.mtz"

def make_inputs():
  o = tm.synthetic(q = (0., 0.5, 0.8))
  with open(PDB, "w") as f: f.write(o.xray_structure.as_pdb_file())
  mtz = o.f_obs.as_mtz_dataset(column_root_label = "FOBS")
  mtz.add_miller_array(miller_array = o.r_free_flags,
    column_root_label = "R-free-flags")
  mtz.mtz_object().write(file_name = MTZ)
  # the model as programs get it (coordinates as written to the file)
  xrs = iotbx.pdb.input(file_name = PDB).xray_structure_simple()
  def r_work(mosaic):
    fm = mmtbx.f_model.manager(f_obs = o.f_obs, r_free_flags = o.r_free_flags,
      xray_structure = xrs)
    fm.update_all_scales(mosaic = mosaic)
    return fm.r_work()
  r_flat, r_mosaic = r_work(False), r_work(True)
  assert r_flat > 0.05 and r_mosaic < r_flat - 0.02, [r_flat, r_mosaic]
  return o, r_flat

class counted(object):
  """
  with counted(r_flat, name) as c: ...; c.r = R-work the program reports
  Checks that the mosaic model was computed and that R is not the one of
  the flat model.
  """
  def __init__(self, r_flat, name):
    self.r_flat, self.name, self.r = r_flat, name, None
  def __enter__(self):
    self.n0 = mosaic.number_of_runs
    return self
  def __exit__(self, exc_type, exc_value, tb):
    if(exc_type is not None): return False
    assert mosaic.number_of_runs > self.n0, self.name
    assert self.r is not None, self.name
    assert self.r < self.r_flat - 0.02, [self.name, self.r, self.r_flat]
    return False

def r_work_in(text, pattern = r"r_work[=:]\s*([\d.]+)"):
  found = re.findall(pattern, text)
  assert len(found) > 0, text[-500:]
  return float(found[-1])

def exercise_maps(r_flat):
  from mmtbx.command_line import maps
  with counted(r_flat, "phenix.maps") as c:
    log = StringIO()
    maps.run(args = [PDB, MTZ], log = log)
    c.r = r_work_in(log.getvalue())
  # no bulk solvent, no mosaic
  n0 = mosaic.number_of_runs
  log = StringIO()
  maps.run(args = [PDB, MTZ, "bulk_solvent_correction=False"], log = log)
  assert mosaic.number_of_runs == n0

def exercise_model_vs_data(r_flat):
  from mmtbx import model_vs_data
  with counted(r_flat, "phenix.model_vs_data") as c:
    log = StringIO()
    model_vs_data.run(args = [PDB, MTZ], log = log)
    c.r = r_work_in(log.getvalue())

def exercise_reciprocal_space_arrays(r_flat):
  from mmtbx.command_line import reciprocal_space_arrays
  from iotbx import reflection_file_reader
  out = "tst_mosaic_programs_arrays.mtz"
  with counted(r_flat, "phenix.reciprocal_space_arrays") as c:
    reciprocal_space_arrays.run(args = [PDB, MTZ,
      "output_file_name=%s"%out], log = StringIO())
    arrays = {}
    for ma in reflection_file_reader.any_reflection_file(
        file_name = out).as_miller_arrays():
      arrays[ma.info().labels[0]] = ma
    f_obs, f_model = arrays["FOBS"].common_sets(arrays["FMODEL"])
    flags = arrays["R_FREE_FLAGS"].common_set(f_obs).data()
    if(not isinstance(flags, flex.bool)): flags = flags != 0
    assert 0 < flags.count(True) < flags.size()//5
    work = ~flags
    fo, fc = f_obs.data().select(work), flex.abs(f_model.data()).select(work)
    c.r = flex.sum(flex.abs(fo-fc))/flex.sum(fo)

def exercise_molprobity(r_flat):
  from mmtbx.command_line import molprobity
  with counted(r_flat, "phenix.molprobity") as c:
    result = molprobity.run(args = [PDB, MTZ], out = StringIO())
    c.r = result.validation.r_work()

def exercise_fem(o, r_flat):
  # The fmodel of the feature-enhanced map (the whole procedure is not for
  # a model of carbon shells)
  from mmtbx.maps import fem
  with counted(r_flat, "phenix.fem") as c:
    obj = fem.run.__new__(fem.run)
    obj.f_obs, obj.r_free_flags = o.f_obs, o.r_free_flags
    obj.xray_structure = iotbx.pdb.input(
      file_name = PDB).xray_structure_simple()
    obj.log = StringIO()
    c.r = obj.create_fmodel().r_work()

def exercise_loader(r_flat):
  """
  mmtbx.command_line.load_model_and_data: mosaic model only if asked for.
  """
  import mmtbx.command_line
  for mosaic_flag in [False, True]:
    n0 = mosaic.number_of_runs
    kw = {}
    if(mosaic_flag): kw["mosaic"] = True
    cmdline = mmtbx.command_line.load_model_and_data(
      args = [PDB, MTZ], master_phil = mmtbx.command_line.generic_simple_input_phil(),
      process_pdb_file = False, out = StringIO(), **kw)
    assert (mosaic.number_of_runs > n0) == mosaic_flag
    r = cmdline.fmodel.r_work()
    if(mosaic_flag): assert r < r_flat - 0.02, [r, r_flat]
    else:            assert approx_equal(r, r_flat, 0.01), [r, r_flat]

def run():
  random.seed(0)
  flex.set_random_seed(0)
  o, r_flat = make_inputs()
  exercise_maps(r_flat)
  exercise_model_vs_data(r_flat)
  exercise_reciprocal_space_arrays(r_flat)
  exercise_loader(r_flat)
  exercise_molprobity(r_flat)
  exercise_fem(o, r_flat)

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %.2f"%(time.time()-t0))
  print("OK")
