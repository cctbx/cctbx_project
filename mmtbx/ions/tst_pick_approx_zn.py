
from __future__ import absolute_import, division, print_function
import os, time
from six.moves import cStringIO as StringIO
from mmtbx.ions.svm.dump_sites import master_phil
from mmtbx.ions.tst_symmetry_axis import get_analyze_waters_result
from iotbx.data_manager import DataManager

def exercise():
  from mmtbx.regression.make_fake_anomalous_data import generate_zinc_inputs
  base = "tst_pick_approx_zn"
  mtz_file, pdb_file = generate_zinc_inputs(file_base=base, anonymize = True)
  time.sleep(2)

  dm = DataManager()
  m = dm.get_model(pdb_file)
  dm.process_miller_array_file(mtz_file)
  fmo = dm.get_fmodel(scattering_table="n_gaussian")
  xrs = m.get_xray_structure()
  xrs.set_inelastic_form_factors(photon = 1.54, table = "sasaki")
  fmo.update_xray_structure(xrs, update_f_calc = True)
  params = master_phil().extract()
  params.input.wavelength = 1.54
  params.use_phaser=False
  params.elements='CA,ZN'
  params.nproc=1
  out = StringIO()
  result = get_analyze_waters_result(m,fmo,params,out)

  n_zn = 0
  for line in out.getvalue().splitlines():
    if "Probable cation: ZN+2" in line:
      n_zn += 1
  if n_zn != 1:
    print(out.getvalue().splitlines())
    raise RuntimeError("Expected 1 ZN+2, found %d" % n_zn)
  os.remove(pdb_file)
  os.remove(mtz_file)
  # "zn_frag_hoh.pdb" => "zn_frag_fmodel.eff"
  os.remove(os.path.splitext(pdb_file)[0][:-4] + ".pdb")
  os.remove(os.path.splitext(pdb_file)[0][:-4] + "_fmodel.eff")
  print("OK")

if (__name__ == "__main__"):
  exercise()
