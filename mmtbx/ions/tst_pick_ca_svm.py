
from __future__ import absolute_import, division, print_function
import os
from six.moves import cStringIO as StringIO
from mmtbx.ions.svm.dump_sites import master_phil
from mmtbx.ions.tst_symmetry_axis import get_analyze_waters_result
from iotbx.data_manager import DataManager
import mmtbx.ions.svm

def exercise():
  try :
    from libsvm import svm      # import dependency
    from libsvm import svmutil  # import dependency
  except ImportError :
    print("libsvm not available, skipping this test")
    return
  from mmtbx.regression.make_fake_anomalous_data import generate_calcium_inputs
  base = "tst_pick_ca_svm"
  mtz_file, pdb_file = generate_calcium_inputs(file_base=base, anonymize=True)

  wavelength = 1.1
  dm = DataManager()
  m = dm.get_model(pdb_file)
  dm.process_miller_array_file(mtz_file)
  fmo = dm.get_fmodel(scattering_table="n_gaussian")
  xrs = m.get_xray_structure()
  xrs.set_inelastic_form_factors(photon = 1.54, table = "sasaki")
  fmo.update_xray_structure(xrs, update_f_calc = True)
  params = master_phil().extract()
  params.input.wavelength = wavelength
  params.use_svm=True
  params.elements='CA,ZN'
  params.use_phaser=False
  params.nproc=1
  manager_class = mmtbx.ions.svm.manager

  out = StringIO()
  result = get_analyze_waters_result(m,fmo,params,out,
    manager_class = manager_class)

  os.remove(pdb_file)
  os.remove(mtz_file)
  # "zn_frag_hoh.pdb" => "zn_frag_fmodel.eff"
  os.remove(os.path.splitext(pdb_file)[0][:-4] + ".pdb")
  os.remove(os.path.splitext(pdb_file)[0][:-4] + "_fmodel.eff")

  n_ca = 0
  for line in out.getvalue().splitlines():
    if ("Final choice: CA" in line):
      n_ca += 1
  if (n_ca != 1):
    print(out.getvalue().splitlines())
    raise RuntimeError("Expected 1 Ca2+, found %d" % n_ca)
  print("OK")

if (__name__ == "__main__"):
  exercise()
