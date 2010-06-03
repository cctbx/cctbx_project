import iotbx
from cctbx import maptbx
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import libtbx.load_env
import sys, os

def exercise_with_tst_input_map():
  file_name = libtbx.env.under_dist(
    module_name="iotbx",
    path="ccp4_map/tst_input.map")
  m = iotbx.ccp4_map.map_reader(file_name=file_name)
  assert approx_equal(m.header_min, -0.422722190619)
  assert approx_equal(m.header_max, 0.335603952408)
  assert approx_equal(m.header_mean, 0)
  assert approx_equal(m.header_rms, 0.140116646886)
  assert m.unit_cell_grid == (16, 8, 16)
  assert approx_equal(m.unit_cell_parameters, (
    82.095001220703125, 37.453998565673828, 69.636001586914062,
    90.0, 101.47599792480469, 90.0))
  assert m.space_group_number == 5
  assert m.data.origin() == (0, 0, 0)
  assert m.data.all() == (16, 8, 16)
  assert not m.data.is_padded()

def exercise(args):
  exercise_with_tst_input_map()
  for file_name in args:
    print file_name
    m = iotbx.ccp4_map.map_reader(file_name=file_name)
    print "header_min: ", m.header_min
    print "header_max: ", m.header_max
    print "header_mean:", m.header_mean
    print "header_rms: ", m.header_rms
    print "unit cell grid:", m.unit_cell_grid
    print "unit cell parameters:", m.unit_cell_parameters
    print "space group number:  ", m.space_group_number
    print "map origin:", m.data.origin()
    print "map grid:  ", m.data.all()
    map_stats = maptbx.statistics(m.data)
    assert approx_equal(map_stats.min(), m.header_min)
    assert approx_equal(map_stats.max(), m.header_max)
    assert approx_equal(map_stats.mean(), m.header_mean)
    if (m.header_rms != 0):
      assert approx_equal(map_stats.sigma(), m.header_rms)
    print

def exercise_writer () :
  file_name = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/wizards/partial_refine_001_map_coeffs.mtz",
    test=os.path.isfile)
  if file_name is None :
    print "Can't find map coefficients file, skipping."
    return
  from iotbx import file_reader
  mtz_in = file_reader.any_file(file_name, force_type="hkl").file_object
  miller_arrays = mtz_in.as_miller_arrays()
  map_coeffs = miller_arrays[0]
  fft_map = map_coeffs.fft_map(resolution_factor=1/3.0)
  fft_map.apply_sigma_scaling()
  real_map = fft_map.real_map_unpadded()
  iotbx.ccp4_map.write_ccp4_map(
    file_name="2mFo-DFc.map",
    map_data=real_map,
    unit_cell=map_coeffs.unit_cell(),
    space_group_number=map_coeffs.space_group().type().number())
  m = iotbx.ccp4_map.map_reader(file_name="2mFo-DFc.map")
  from scitbx.array_family import flex
  mmm = flex.double(list(real_map)).min_max_mean()
  assert approx_equal(m.unit_cell_parameters,
                      map_coeffs.unit_cell().parameters())
  assert approx_equal(mmm.min, m.header_min)
  assert approx_equal(mmm.max, m.header_max)
  assert approx_equal(mmm.mean, m.header_mean)

def run(args):
  def have_ext():
    for node in os.listdir(libtbx.env.under_build(path="lib")):
      if (node.startswith("iotbx_ccp4_map_ext")):
        return True
    return False
  if (not have_ext()): # XXX backward compatibility 2008-09-30
    print "Skipping iotbx_ccp4_map tests: extension not available"
  else:
    import iotbx.ccp4_map
    exercise(args=args)
    exercise_writer()
  print format_cpu_times()

if (__name__ == "__main__"):
  run(sys.argv[1:])
