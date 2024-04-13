
from __future__ import absolute_import, division, print_function
import os, time
from six.moves import cStringIO as StringIO
from iotbx.data_manager import DataManager
from mmtbx.regression.make_fake_anomalous_data import generate_magnessium_inputs
import iotbx.phil
import mmtbx.ions.identify

def exercise():
  base = "tst_validate_mg"
  mtz_file, pdb_file = generate_magnessium_inputs(file_base=base, anonymize=False)
  time.sleep(2)
  dm = DataManager()
  m = dm.get_model(pdb_file)
  m.process(make_restraints=True)
  grm = m.get_restraints_manager()
  ma = dm.get_miller_arrays(filename = mtz_file)
  fmo = dm.get_fmodel(scattering_table="n_gaussian")
  from mmtbx.ions.identify import ion_identification_phil_str
  params = iotbx.phil.parse(input_string = ion_identification_phil_str).extract()
  out = StringIO()
  manager = mmtbx.ions.identify.create_manager(
    pdb_hierarchy = m.get_hierarchy(),
    fmodel=fmo,
    geometry_restraints_manager=grm.geometry,
    wavelength=None,
    params=params,
    nproc = 1,
    log= out)
  result = manager.validate_ions(out = out)
  n_mg, n_bad = 0, 0
  for line in out.getvalue().splitlines():
    if "| MG" in line:
      n_mg += 1
    if "!!!" in line:
      n_bad += 1
  assert n_mg == 2 and n_bad == 0
  for ext in [".pdb", ".mtz", "_fmodel.eff"]:
    os.remove(base + ext)
  print("OK")

if (__name__ == "__main__"):
  exercise()
