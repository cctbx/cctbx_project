
from __future__ import absolute_import, division, print_function
import os,time
from six.moves import cStringIO as StringIO
from mmtbx.ions.svm.dump_sites import master_phil
from mmtbx.ions.tst_symmetry_axis import get_analyze_waters_result
from iotbx.data_manager import DataManager

def exercise():
  from mmtbx.regression.make_fake_anomalous_data import generate_calcium_inputs
  base = "tst_pick_ca"
  mtz_file, pdb_file = generate_calcium_inputs(file_base=base, anonymize=True)
  time.sleep(2)

  wavelength = 1.1
  dm = DataManager()
  m = dm.get_model(pdb_file)
  dm.process_miller_array_file(mtz_file)
  fmo = dm.get_fmodel(scattering_table="n_gaussian")
  params = master_phil().extract()
  params.input.wavelength = 1.1
  params.phaser.fpp_ratio_max=1.2
  params.use_phaser=False
  params.nproc=1
  out = StringIO()
  result = get_analyze_waters_result(m,fmo,params,out)

  n_ca = 0
  for line in out.getvalue().splitlines():
    if "Probable cation: CA+2" in line:
      n_ca += 1
  if (n_ca != 1):
    print(out.getvalue().splitlines())
    raise RuntimeError("Expected 1 Ca2+, found %d" % n_ca)
  os.remove(pdb_file)
  os.remove(mtz_file)
  os.remove(os.path.splitext(pdb_file)[0][:-4] + ".pdb")
  os.remove(os.path.splitext(pdb_file)[0][:-4] + "_fmodel.eff")
  print("OK")

if (__name__ == "__main__"):
  exercise()
