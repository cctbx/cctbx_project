from iotbx import mtz
from iotbx.mtz.dump import dump
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex
from scitbx.test_utils import approx_equal
import random
import sys

def to_mtz(miller_array, mtz_label):
  w = mtz.MtzWriter()
  w.setTitle("mtz writer test")
  w.setSpaceGroup(miller_array.space_group_info())
  w.oneCrystal("test_crystal","test_project",miller_array.unit_cell())
  wavelength = 1.0
  w.oneDataset("test_dataset",wavelength)
  w.add_miller_array(miller_array, mtz_label)
  return w

def recycle(miller_array, mtz_label, verbose=0):
  w = to_mtz(miller_array, mtz_label)
  w.write("tmp.mtz")
  label_sigmas = w.label_sigmas(mtz_label)
  label_phases = w.label_phases(mtz_label)
  if (0 or verbose):
    p = dump("tmp.mtz")
  else:
    p = mtz.Mtz("tmp.mtz")
  assert p.title() == "mtz writer test"
  assert p.nsym() == miller_array.space_group().order_z()
  assert p.ncrystals() == 1
  assert p.getSgtbxSpaceGroup() == miller_array.space_group()
  assert p.history().size() == 0
  cryst = p.getCrystal(0)
  assert cryst.crystal_name() == "test_crystal"
  assert cryst.project_name() == "test_project"
  assert cryst.UnitCell().is_similar_to(miller_array.unit_cell())
  crystal_symmetry = crystal.symmetry(
    unit_cell=cryst.UnitCell(),
    space_group_info=p.get_space_group_info())
  assert cryst.ndatasets() == 1
  dataset = cryst.getDataset(0)
  assert dataset.dataset_name() == "test_dataset"
  assert dataset.wavelength() - 1.0 < 1.e-5
  if (not miller_array.anomalous_flag()):
    if (miller_array.sigmas() == None):
      if (miller_array.is_complex()):
        assert dataset.ncolumns() == 3+2
        r = miller.array(
          miller_set=miller.set(
            crystal_symmetry=crystal_symmetry,
            indices=p.valid_indices(mtz_label),
            anomalous_flag=False),
          data=p.valid_complex(mtz_label, label_phases))
      else:
        assert dataset.ncolumns() == 3+1
        r = miller.array(
          miller_set=miller.set(
            crystal_symmetry=crystal_symmetry,
            indices=p.valid_indices(mtz_label),
            anomalous_flag=False),
          data=p.valid_values(mtz_label))
    else:
      assert dataset.ncolumns() == 3+2
      r = miller.array(
        miller_set=miller.set(
          crystal_symmetry=crystal_symmetry,
          indices=p.valid_indices(mtz_label),
          anomalous_flag=False),
        data=p.valid_values(mtz_label),
        sigmas=p.valid_values(label_sigmas))
  else:
    if (miller_array.sigmas() == None):
      if (miller_array.is_complex()):
        assert dataset.ncolumns() == 3+4
        r = miller.array(
          miller_set=miller.set(
            crystal_symmetry=crystal_symmetry,
            indices=p.valid_indices(
              w.label_plus(mtz_label), w.label_minus(mtz_label)),
            anomalous_flag=True),
          data=p.valid_complex(
            w.label_plus(mtz_label), w.label_plus(label_phases),
            w.label_minus(mtz_label), w.label_minus(label_phases)))
      else:
        assert dataset.ncolumns() == 3+2
        r = miller.array(
          miller_set=miller.set(
            crystal_symmetry=crystal_symmetry,
            indices=p.valid_indices(
              w.label_plus(mtz_label), w.label_minus(mtz_label)),
            anomalous_flag=True),
          data=p.valid_values(
            w.label_plus(mtz_label), w.label_minus(mtz_label)))
    else:
      assert dataset.ncolumns() == 3+4
      r = miller.array(
        miller_set=miller.set(
          crystal_symmetry=crystal_symmetry,
          indices=p.valid_indices(
            w.label_plus(mtz_label), w.label_minus(mtz_label)),
          anomalous_flag=True),
        data=p.valid_values(
          w.label_plus(mtz_label), w.label_minus(mtz_label)),
        sigmas=p.valid_values(
          w.label_plus(label_sigmas), w.label_minus(label_sigmas)))
  verify_miller_arrays(miller_array, r)
  r = p.as_miller_arrays()
  assert len(r) == 1
  verify_miller_arrays(miller_array, r[0])

def verify_miller_arrays(a1, a2):
  v = a2.adopt_set(a1)
  assert flex.max(flex.abs(a1.data() - v.data())) < 1.e-5
  if (v.sigmas() != None):
    assert flex.max(flex.abs(a1.sigmas() - v.sigmas())) < 1.e-5

def exercise(space_group_info, n_scatterers=8, d_min=5,
             anomalous_flag=False, verbose=0):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["const"]*n_scatterers,
    volume_per_atom=500,
    min_distance=2.,
    general_positions_only=True)
  if (0 or verbose):
    structure.show_summary().show_scatterers()
  f_calc = structure.structure_factors(
    d_min=d_min, anomalous_flag=anomalous_flag).f_calc()
  f_calc = miller.array(
    miller_set=f_calc,
    data=f_calc.data()/flex.mean(flex.abs(f_calc.data())))
  if (f_calc.anomalous_flag()):
    selection = flex.bool(f_calc.indices().size(), True)
    for i in xrange(f_calc.indices().size()/10):
      j = random.randrange(f_calc.indices().size())
      selection[j] = False
    f_calc = f_calc.apply_selection(selection)
  recycle(f_calc, "f_calc", verbose=verbose)
  recycle(abs(f_calc), "f_obs", verbose=verbose)
  recycle(miller.array(
    miller_set=f_calc,
    data=flex.abs(f_calc.data()),
    sigmas=flex.abs(f_calc.data())/10), "f_obs", verbose=verbose)

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
