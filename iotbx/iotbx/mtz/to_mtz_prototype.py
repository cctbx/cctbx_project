from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx import miller
from cctbx.array_family import flex
import math
import random
import sys
from iotbx.mtz import MtzWriter

def to_mtz(miller_array, label_data, label_sigmas=None):
  miller_array.show_summary()
  W = MtzWriter()
  W.setTitle("mtz writer test")
  W.setSpaceGroup(miller_array.space_group())
  W.oneCrystal("test_crystal","test_project",miller_array.unit_cell())
  W.oneDataset("test_dataset",1.00)#set the wavelength

  anomalous_flag = miller_array.anomalous_flag()
  miller_indices = miller_array.indices()
  data = miller_array.data()
  sigmas = miller_array.sigmas() # none if no sigmas
  assert [sigmas, label_sigmas].count(None) in (0, 2)
  if (anomalous_flag):
    asu, matches = miller_array.match_bijvoet_mates()
    for i,j in matches.pairs():
      for h,k in zip(asu.indices()[i], asu.indices()[j]): assert h == -k
      print asu.indices()[i], asu.indices()[j],
      print asu.data()[i], asu.data()[j]
    reciprocal_space_asu = asu.space_group_info().reciprocal_space_asu()
    for i in matches.singles():
      h = asu.indices()[i]
      f = asu.data()[i]
      s = reciprocal_space_asu.which(h)
      assert s != 0
      if (reciprocal_space_asu.which(h) > 0):
        print "+", h, f
      else:
        print "-", h, f
  else:
    for i,h in miller_array.indices().items():
      print h, miller_array.data()[i]
  return W

def exercise(space_group_info, n_scatterers=8, d_min=5,
             anomalous_flag=False, verbose=0):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["const"]*n_scatterers,
    volume_per_atom=500,
    min_distance=2.,
    general_positions_only=0001)
  if (0 or verbose):
    structure.show_summary().show_scatterers()
  f_calc = structure.structure_factors(
    d_min=d_min, anomalous_flag=anomalous_flag).f_calc()
  if (f_calc.anomalous_flag()):
    selection = flex.bool(f_calc.indices().size(), 0001)
    for i in xrange(f_calc.indices().size()/10):
      j = random.randrange(f_calc.indices().size())
      selection[j] = 00000
    f_calc = f_calc.apply_selection(selection)
  to_mtz(f_calc, "f_calc").write('f_calc.mtz')
  to_mtz(abs(f_calc), "f_obs").write('f_obs.mtz')
  to_mtz(miller.array(
    miller_set=f_calc,
    data=flex.abs(f_calc.data()),
    sigmas=flex.abs(f_calc.data())/10), "f_obs", "sigma").write('sigma.mtz')

def run_call_back(flags, space_group_info):
  for anomalous_flag in (False, True):
    exercise(
      space_group_info,
      anomalous_flag=anomalous_flag,
      verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
