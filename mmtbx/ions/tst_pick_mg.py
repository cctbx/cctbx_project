
from __future__ import absolute_import, division, print_function
import os, time
from six.moves import cStringIO as StringIO
from mmtbx.ions.svm.dump_sites import master_phil
from mmtbx.ions.tst_symmetry_axis import get_analyze_waters_result
from iotbx.data_manager import DataManager


def exercise():
  from mmtbx.regression.make_fake_anomalous_data import generate_magnessium_inputs
  base = "tst_pick_mg"
  mtz_file, pdb_file = generate_magnessium_inputs(file_base=base, anonymize=True)
  time.sleep(2)

  dm = DataManager()
  m = dm.get_model(pdb_file)
  dm.process_miller_array_file(mtz_file)
  fmo = dm.get_fmodel(scattering_table="n_gaussian")
  params = master_phil().extract()
  params.use_phaser=False
  params.elements='MG'
  out = StringIO()
  result = get_analyze_waters_result(m,fmo,params,out)

  n_mg = 0
  for line in out.getvalue().splitlines():
    if ("Probable cation: MG+2" in line):
      n_mg += 1
  if n_mg != 2:
    print("\n".join(result.stdout_lines))
    raise RuntimeError("Expected 2 MG+2, found %d" % n_mg)
  os.remove(pdb_file)
  os.remove(mtz_file)
  # "zn_frag_hoh.pdb" => "zn_frag_fmodel.eff"
  os.remove(os.path.splitext(pdb_file)[0][:-4] + ".pdb")
  os.remove(os.path.splitext(pdb_file)[0][:-4] + "_fmodel.eff")
  print("OK")

if __name__ == "__main__":
  exercise()
