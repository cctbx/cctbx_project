from iotbx import mtz
from iotbx.mtz.dump import dump
from iotbx.regression.utils import random_f_calc
from cctbx.development import debug_utils
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex
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
  assert p.ncrystals() == 2
  assert p.getSgtbxSpaceGroup() == miller_array.space_group()
  assert p.history().size() == 0
  cryst = p.getCrystal(0)
  assert cryst.crystal_name() == "HKL_base"
  assert cryst.project_name() == "HKL_base"
  cryst = p.getCrystal(1)
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
    if (miller_array.sigmas() is None):
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
    if (miller_array.sigmas() is None):
      if (miller_array.is_complex()):
        assert dataset.ncolumns() == 3+4
        r = miller.array(
          miller_set=miller.set(
            crystal_symmetry=crystal_symmetry,
            indices=p.valid_indices_anomalous(
              w.label_plus(mtz_label), w.label_minus(mtz_label)),
            anomalous_flag=True),
          data=p.valid_complex_anomalous(
            w.label_plus(mtz_label), w.label_plus(label_phases),
            w.label_minus(mtz_label), w.label_minus(label_phases)))
      else:
        assert dataset.ncolumns() == 3+2
        r = miller.array(
          miller_set=miller.set(
            crystal_symmetry=crystal_symmetry,
            indices=p.valid_indices_anomalous(
              w.label_plus(mtz_label), w.label_minus(mtz_label)),
            anomalous_flag=True),
          data=p.valid_values_anomalous(
            w.label_plus(mtz_label), w.label_minus(mtz_label)))
    else:
      assert dataset.ncolumns() == 3+4
      r = miller.array(
        miller_set=miller.set(
          crystal_symmetry=crystal_symmetry,
          indices=p.valid_indices_anomalous(
            w.label_plus(mtz_label), w.label_minus(mtz_label)),
          anomalous_flag=True),
        data=p.valid_values_anomalous(
          w.label_plus(mtz_label), w.label_minus(mtz_label)),
        sigmas=p.valid_values_anomalous(
          w.label_plus(label_sigmas), w.label_minus(label_sigmas)))
  verify_miller_arrays(miller_array, r)
  r = p.as_miller_arrays()
  assert len(r) == 1
  verify_miller_arrays(miller_array, r[0])
  miller_array.export_as_mtz("tmp.mtz", column_label="data")

def verify_miller_arrays(a1, a2):
  v = a2.adopt_set(a1)
  assert flex.max(flex.abs(a1.data() - v.data())) < 1.e-5
  if (v.sigmas() is not None):
    assert flex.max(flex.abs(a1.sigmas() - v.sigmas())) < 1.e-5

def exercise(space_group_info, n_scatterers=8, d_min=2.5,
             anomalous_flag=False, verbose=0):
  f_calc = random_f_calc(
    space_group_info=space_group_info,
    n_scatterers=n_scatterers,
    d_min=d_min,
    anomalous_flag=anomalous_flag,
    verbose=verbose)
  if (f_calc is None): return
  recycle(f_calc, "f_calc", verbose=verbose)
  recycle(abs(f_calc), "f_obs", verbose=verbose)
  recycle(miller.array(
    miller_set=f_calc,
    data=flex.abs(f_calc.data()),
    sigmas=flex.abs(f_calc.data())/10), "f_obs", verbose=verbose)
  to_mtz(
    miller_array=miller.array(
      miller_set=miller.set(
        crystal_symmetry=f_calc,
        anomalous_flag=anomalous_flag,
        indices=f_calc.indices()),
      data=flex.abs(f_calc.data())),
    mtz_label="f_obs").write("tmp.mtz")
  arrays = mtz.Mtz("tmp.mtz").as_miller_arrays()
  assert len(arrays) == 1
  assert arrays[0].anomalous_flag() == anomalous_flag

def run_call_back(flags, space_group_info):
  for anomalous_flag in (False, True):
    exercise(
      space_group_info,
      anomalous_flag=anomalous_flag,
      verbose=flags.Verbose)

def exercise_read_corrupt():
  for i_trial in xrange(2):
    f = open("tmp.mtz", "wb")
    if (i_trial > 0):
      f.write("\0"*80)
    f.close()
    try: mtz.Mtz("tmp.mtz")
    except RuntimeError, e: assert str(e) == "Mtz read failed"
    else: raise AssertionError("Exception expected.")

def run():
  assert mtz.ccp4_liberr_verbosity(-1) == 0
  assert mtz.ccp4_liberr_verbosity(1) == 1
  assert mtz.ccp4_liberr_verbosity(-1) == 1
  assert mtz.ccp4_liberr_verbosity(0) == 0
  assert mtz.ccp4_liberr_verbosity(-1) == 0
  exercise_read_corrupt()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
