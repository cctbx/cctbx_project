
from __future__ import absolute_import, division, print_function
import os, time
from six.moves import cStringIO as StringIO
from iotbx.data_manager import DataManager
from mmtbx.regression.make_fake_anomalous_data import generate_calcium_inputs
import iotbx.phil
import mmtbx.ions.identify

def exercise():
  base = "tst_validate_ca"
  mtz_file, pdb_file = generate_calcium_inputs(file_base=base, anonymize=False)
  time.sleep(2)
  dm = DataManager()
  m = dm.get_model(pdb_file)
  m.process(make_restraints=True)
  grm = m.get_restraints_manager()
  dm.process_miller_array_file(mtz_file)
  fmo = dm.get_fmodel(scattering_table="n_gaussian")
  from mmtbx.ions.identify import ion_identification_phil_str
  params = iotbx.phil.parse(input_string = ion_identification_phil_str).extract()
  out = StringIO()
  manager = mmtbx.ions.identify.create_manager(
    pdb_hierarchy = m.get_hierarchy(),
    fmodel=fmo,
    geometry_restraints_manager=grm.geometry,
    wavelength=1.12,
    params=params,
    nproc = 1,
    log= out)
  result = manager.validate_ions(out = out)
  n_ca, n_bad = 0, 0
  for line in out.getvalue().splitlines():
    if "| CA" in line:
      n_ca += 1
    if "!!!" in line:
      n_bad += 1
  assert n_ca == 1 and n_bad == 0
  for ext in [".pdb", ".mtz", "_fmodel.eff"]:
    os.remove(base + ext)
  print("OK")

if (__name__ == "__main__"):
  #print("WARNING: TEST TOO SLOW. MAKE IT RUN UNDER 300s AND ENABLE BACK.")
  t0 = time.time()
  exercise()
  print("Time: %6.2f"%(time.time()-t0))
  print("OK")
