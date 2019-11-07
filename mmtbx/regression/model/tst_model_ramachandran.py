from __future__ import absolute_import, division, print_function
import mmtbx.model
import iotbx.pdb

from tst_model_remove_alternative_conformations import tst_pdb_str

def exercise_1(prefix="tst_model_ramachandran_1"):
  """ Make sure that set_ramachandran_plot_restraints() works. """
  inp = iotbx.pdb.input(lines=tst_pdb_str, source_info=None)
  model = mmtbx.model.manager(
      model_input = inp)
  model.set_ramachandran_plot_restraints(rama_potential='oldfield')
  n = model.get_restraints_manager().geometry.ramachandran_manager.get_n_oldfield_proxies()
  assert n == 2

if __name__ == '__main__':
  exercise_1()
