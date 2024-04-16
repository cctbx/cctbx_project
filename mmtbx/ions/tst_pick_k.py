
from __future__ import absolute_import, division, print_function
import os, time
from libtbx import easy_run
from six.moves import cStringIO as StringIO
from mmtbx.ions.svm.dump_sites import master_phil
from mmtbx.ions.tst_symmetry_axis import get_analyze_waters_result
from iotbx.data_manager import DataManager

def exercise():
  # FIXME
  print("Temporarily disabled, skipping")
  return
  # the function below does not exist????
  from mmtbx.regression.make_fake_anomalous_data import generate_potassium_inputs
  base = "tst_pick_k"
  mtz_file, pdb_file = generate_potassium_inputs(file_base=base, anonymize=True)
  time.sleep(2)

#  dm = DataManager()
#  m = dm.get_model(pdb_file)
#  dm.process_miller_array_file(mtz_file)
#  fmo = dm.get_fmodel(scattering_table="n_gaussian")
#  params = master_phil().extract()
#  params.elements='K,MG'
#  params.use_phaser=False
#  params.nproc=1
#  out = StringIO()
#  result = get_analyze_waters_result(m,fmo,params,out)

  args = [pdb_file, mtz_file, "nproc=1", "elements=K,MG", "use_phaser=False"]
  result = easy_run.fully_buffered("mmtbx.water_screen %s" % " ".join(args)
    ).raise_if_errors()
  n_k = 0
  for line in result.stdout_lines :
    if ("Probable cation: K+1" in line):
      n_k += 1
  assert n_k == 3
  os.remove(pdb_file)
  os.remove(mtz_file)
  # "zn_frag_hoh.pdb" => "zn_frag_fmodel.eff"
  os.remove(os.path.splitext(pdb_file)[0][:-4] + "_fmodel.eff")
  print("OK")

if (__name__ == "__main__"):
  exercise()
