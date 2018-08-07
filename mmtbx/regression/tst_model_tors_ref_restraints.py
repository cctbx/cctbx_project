from __future__ import division
import mmtbx.model
import libtbx.load_env
from libtbx.utils import format_cpu_times
from libtbx.test_utils import approx_equal
import iotbx.pdb
from mmtbx.regression import tst_model_cart_ref_restraints


def exercise_adopting_ref_tors_restraints_h():
  inp_1 = iotbx.pdb.input(lines=tst_model_cart_ref_restraints.pdb_str_h, source_info=None)
  h_model = mmtbx.model.manager(model_input = inp_1)

  inp_2 = iotbx.pdb.input(lines=tst_model_cart_ref_restraints.pdb_str, source_info=None)
  model = mmtbx.model.manager(model_input = inp_2)

  # case 1: big model adopting small model.
  h_model.set_reference_torsion_restraints(
      ref_model = model)


def run():
  if (not libtbx.env.has_module("reduce")):
    print "Reduce not installed"
    return
  exercise_adopting_ref_tors_restraints_h()

  print format_cpu_times()

if (__name__ == "__main__"):
  run()
