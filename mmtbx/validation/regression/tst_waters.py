
from __future__ import absolute_import, division, print_function
from six.moves import cStringIO as StringIO
from libtbx.test_utils import show_diff
from libtbx.utils import null_out
from libtbx import easy_pickle
from libtbx import group_args
from iotbx.data_manager import DataManager
from mmtbx.regression import make_fake_anomalous_data
import mmtbx.ions.utils
import iotbx.pdb
import mmtbx.model
from mmtbx.validation import waters

def exercise_heavy():

  file_base = "tst_validate_waters_1"
  pdb_file = make_fake_anomalous_data.write_pdb_input_cd_cl(file_base=file_base)
  mtz_file = make_fake_anomalous_data.generate_mtz_file(
    file_base="tst_validate_waters_1",
    d_min=1.5,
    anomalous_scatterers=[
      group_args(selection="element CD", fp=-0.29, fdp=2.676),
      group_args(selection="element CL", fp=0.256, fdp=0.5),
    ])

  pdb_in = iotbx.pdb.input(pdb_file)
  m1 = mmtbx.model.manager(model_input = pdb_in, log = null_out())
  hierarchy, n = mmtbx.ions.utils.anonymize_ions(m1.get_hierarchy(),
                                                 log=null_out())
  fn_anonymized = file_base + '_start.pdb'
  hierarchy.write_pdb_file(fn_anonymized,crystal_symmetry=m1.crystal_symmetry())
  dm = DataManager()
  m = dm.get_model(fn_anonymized)
  ma = dm.get_miller_arrays(filename = mtz_file)
  fmo = dm.get_fmodel(scattering_table="n_gaussian")
  results = waters.waters(
    pdb_hierarchy=m.get_hierarchy(),
    xray_structure=m.get_xray_structure(),
    fmodel=fmo,
    collect_all=True)
  out = StringIO()
  results.show(out=out)
  s = easy_pickle.dumps(results)
  r2 = easy_pickle.loads(s)
  out2 = StringIO()
  r2.show(out=out2)
  assert not show_diff(out.getvalue(), out2.getvalue())
  assert (results.n_bad >= 1) and (results.n_heavy == 2)
  # XXX statistics are approximate (probably a precision issue), so we can't
  # directly compare outputs

if (__name__ == "__main__"):
  exercise_heavy()
  print("OK")
