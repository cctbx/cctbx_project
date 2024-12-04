from __future__ import division, print_function
from scitbx.array_family import flex
import iotbx.pdb
from iotbx.data_manager import DataManager
import mmtbx.model
from libtbx.utils import null_out
import mmtbx.f_model

def get_r(x,y):
  x = abs(x).data()
  y = abs(y).data()
  scale = flex.sum(x*y)/flex.sum(y*y)
  return flex.sum(flex.abs(x-scale*y))/flex.sum(flex.abs(x+scale*y))*2*100.

pdb_str = """
CRYST1    4.870   10.060   30.660  94.85  90.26  99.99 P 1           1
ATOM      1  N   GLN A   1       1.153   1.465  -2.445  1.00  9.04           N
ATOM      2  CA  GLN A   1       1.681   2.444  -1.508  1.00  4.57           C
ATOM      3  C   GLN A   1       1.010   2.351  -0.141  1.00  3.21           C
ATOM      4  O   GLN A   1      -0.214   2.379  -0.045  1.00  5.48           O
ATOM      5  CB  GLN A   1       1.514   3.851  -2.090  1.00  6.45           C
ATOM      6  CG  GLN A   1       2.245   4.927  -1.325  1.00  8.01           C
ATOM      7  CD  GLN A   1       1.858   6.297  -1.796  1.00  5.89           C
ATOM      8  OE1 GLN A   1       0.717   6.720  -1.624  1.00  7.78           O
ATOM      9  NE2 GLN A   1       2.793   6.997  -2.405  1.00  6.95           N
"""

phil_str = '''
data_manager {
    model {
      file = %s
      type = x_ray neutron *electron reference
    }
    miller_array {
      file = %s
      labels {
        name = %s
        type = x_ray neutron *electron
      }
      labels {
        name = %s
        type = x_ray neutron *electron
      }
  }
}'''

def run():
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model.setup_scattering_dictionaries(scattering_table="electron")
  fc_1 = model.get_xray_structure().structure_factors(d_min=1,
    algorithm="direct").f_calc()
  model.setup_scattering_dictionaries(scattering_table="n_gaussian")
  fc_2 = fc_1.structure_factors_from_scatterers(
    xray_structure = model.get_xray_structure()).f_calc()
  r = get_r(fc_1, fc_2)
  assert r>10
  model.setup_scattering_dictionaries(scattering_table="electron")
  fc_2 = fc_1.structure_factors_from_scatterers(
    xray_structure = model.get_xray_structure(), algorithm="direct").f_calc()
  r = get_r(fc_1, fc_2)
  assert r<1.e-6
  #
  pdb_fn = 'tst_fmodel_and_dm.pdb'
  mtz_fn = 'tst_fmodel_and_dm.mtz'
  label_fobs = 'FOBS'
  label_flag = 'R-free-flags'
  with open(pdb_fn,"w") as fo:
    fo.write(model.model_as_pdb())
  mtz_dataset = abs(fc_2).as_mtz_dataset(column_root_label = 'FOBS')
  mtz_dataset.add_miller_array(
    miller_array      = fc_2.generate_r_free_flags(),
    column_root_label = label_flag)
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = mtz_fn)
  dm_phil = iotbx.phil.parse(phil_str%(pdb_fn, mtz_fn, label_fobs, label_flag))
  dm = DataManager(['model', 'miller_array'])
  dm.load_phil_scope(dm_phil)
  fmodel = dm.get_fmodel(scattering_table="electron")
  params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  params.algorithm="direct"
  fmodel.update(sf_and_grads_accuracy_params=params)
  fmodel.update_xray_structure(update_f_calc=True)
  fc_3 = fmodel.f_calc()
  r = get_r(fc_2, fc_3)
  assert r<1.e-4, r # we lost a couple of digits due to file reading?

if(__name__ == "__main__"):
  run()
