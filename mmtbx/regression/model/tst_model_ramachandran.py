from __future__ import absolute_import, division, print_function
import mmtbx.model
import iotbx.pdb

from mmtbx.geometry_restraints.ramachandran import master_phil

from tst_model_remove_alternative_conformations import tst_pdb_str

def exercise_1(prefix="tst_model_ramachandran_1"):
  """ Make sure that set_ramachandran_plot_restraints() works. """
  inp = iotbx.pdb.input(lines=tst_pdb_str, source_info=None)
  model = mmtbx.model.manager(
      model_input = inp)
  model.process(make_restraints=True)
  rama_params = master_phil.fetch().extract().ramachandran_plot_restraints
  rama_params.favored = 'oldfield'
  rama_params.allowed = 'oldfield'
  rama_params.outlier = 'oldfield'
  model.set_ramachandran_plot_restraints(rama_params=rama_params)
  model.process(make_restraints=True)
  n = model.get_restraints_manager().geometry.ramachandran_manager.get_n_oldfield_proxies()
  assert n == 2

  rama_params.favored = 'emsley'
  rama_params.allowed = 'emsley'
  rama_params.outlier = 'emsley'
  model.set_ramachandran_plot_restraints(rama_params=rama_params)
  model.process(make_restraints=True)
  n = model.get_restraints_manager().geometry.ramachandran_manager.get_n_emsley_proxies()
  assert n == 2

if __name__ == '__main__':
  exercise_1()
