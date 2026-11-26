from cctbx.array_family import flex
from mmtbx.bulk_solvent import tst_squeeze as ts
from six.moves import cStringIO as StringIO
import mmtbx.f_model
import iotbx.pdb
import random, re, time

"""
Programs that can use squeeze (mmtbx.bulk_solvent.squeeze): off by default,
and with squeeze=True they report the R-factor of the model with it.
Simulated data: blobs of density in the solvent that the model does not have.
"""

PDB, MTZ = "tst_squeeze_programs.pdb", "tst_squeeze_programs.mtz"

def make_inputs():
  o = ts.synthetic()
  ts.write_files(o = o, pdb_file_name = PDB, mtz_file_name = MTZ)
  # the model as programs get it
  pdb_inp = iotbx.pdb.input(file_name = PDB)
  h = pdb_inp.construct_hierarchy()
  water = h.atom_selection_cache().selection("water")
  assert water.count(True) == o.water_selection.count(True) > 0
  def r_work(squeeze):
    fm = mmtbx.f_model.manager(f_obs = o.f_obs, r_free_flags = o.r_free_flags,
      xray_structure = pdb_inp.xray_structure_simple())
    fm.update_all_scales(squeeze = squeeze, water_selection = water)
    return fm.r_work()
  r_off, r_on = r_work(False), r_work(True)
  assert r_off > 0.08 and r_on < 0.7*r_off, [r_off, r_on]
  return r_off

def r_work_in(text, pattern = r"r_work[=:]\s*([\d.]+)"):
  found = re.findall(pattern, text)
  assert len(found) > 0, text[-500:]
  return float(found[-1])

def check(name, r_off, r_default, r_squeeze):
  assert abs(r_default - r_off) < 0.01, [name, r_default, r_off]
  assert r_squeeze < 0.7*r_off, [name, r_squeeze, r_off]

def exercise_maps(r_off):
  from mmtbx.command_line import maps
  r = []
  for args in [[], ["squeeze=True"]]:
    log = StringIO()
    maps.run(args = [PDB, MTZ] + args, log = log)
    r.append(r_work_in(log.getvalue()))
    has = "Squeeze: unmodelled density" in log.getvalue()
    assert has == (len(args) == 1)
  check("phenix.maps", r_off, r[0], r[1])
  # not applicable: said so, maps computed
  log = StringIO()
  maps.run(args = [PDB, MTZ, "squeeze=True", "high_resolution=1.5"],
    log = log)
  assert "squeeze: not used (data resolution 1.50 A" in log.getvalue()

def exercise_model_vs_data(r_off):
  from mmtbx import model_vs_data
  r = []
  for args in [[], ["squeeze=True"]]:
    log = StringIO()
    model_vs_data.run(args = [PDB, MTZ] + args, log = log)
    r.append(r_work_in(log.getvalue()))
  check("phenix.model_vs_data", r_off, r[0], r[1])

def exercise_loader(r_off):
  import mmtbx.command_line
  r = []
  for kw in [{}, {"squeeze":True}]:
    cmdline = mmtbx.command_line.load_model_and_data(
      args = [PDB, MTZ],
      master_phil = mmtbx.command_line.generic_simple_input_phil(),
      process_pdb_file = False, out = StringIO(), **kw)
    fm = cmdline.fmodel
    assert (fm.squeeze_result() is not None) == ("squeeze" in kw)
    if("squeeze" in kw): assert fm.squeeze_result().water_excluded
    r.append(fm.r_work())
  check("load_model_and_data", r_off, r[0], r[1])

def run():
  random.seed(0)
  flex.set_random_seed(0)
  r_off = make_inputs()
  exercise_maps(r_off)
  exercise_model_vs_data(r_off)
  exercise_loader(r_off)

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.2f"%(time.time()-t0))
  print("OK")
