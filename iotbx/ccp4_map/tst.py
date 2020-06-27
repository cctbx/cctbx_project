from __future__ import absolute_import, division, print_function
import iotbx
from cctbx import maptbx
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import libtbx.load_env
from six.moves import cStringIO as StringIO
import sys
from iotbx import mrcfile
from six.moves import zip


def add_list(list_a,list_b):
  assert len(list_a)==len(list_b)
  new_list=[]
  for a,b in zip(list_a,list_b):
    new_list.append(a+b)
  return new_list

def subtract_list(list_a,list_b):
  assert len(list_a)==len(list_b)
  new_list=[]
  for a,b in zip(list_a,list_b):
    new_list.append(a-b)
  return new_list

def exercise_with_tst_input_map(use_mrcfile=None,file_name=None):
  if not file_name:
    file_name = libtbx.env.under_dist(
      module_name="iotbx",
      path="ccp4_map/tst_input.map")

  print("\nTesting read of input map with axis order 3 1 2 "+\
         "\n and use_mrcfile=",use_mrcfile)
  if use_mrcfile:
    m = mrcfile.map_reader(file_name=file_name,verbose=True)
    for label in m.labels: print(label)
  else:
    m = iotbx.ccp4_map.map_reader(file_name=file_name)
  assert approx_equal(m.header_min, -0.422722190619)
  assert approx_equal(m.header_max, 0.335603952408)
  assert approx_equal(m.header_mean, 0)
  assert approx_equal(m.header_rms, 0.140116646886)
  assert m.unit_cell_grid == (16, 8, 16)
  assert approx_equal(m.unit_cell().parameters(), (
    82.095001220703125, 37.453998565673828, 69.636001586914062,
    90.0, 101.47599792480469, 90.0))
  assert m.unit_cell_crystal_symmetry().space_group_number()==5
  assert m.map_data().origin() == (0, 0, 0)
  assert m.map_data().all() == (16, 8, 16)
  assert approx_equal(m.pixel_sizes(),(5.130937576293945, 4.6817498207092285, 4.352250099182129))
  assert not m.map_data().is_padded()
  out = StringIO()
  m.show_summary(out=out)
  assert ("map cell grid: (16, 8, 16)" in out.getvalue())
  uc = m.unit_cell_crystal_symmetry().unit_cell()
  assert approx_equal(m.unit_cell().parameters(), m.unit_cell().parameters())
  assert approx_equal(m.grid_unit_cell().parameters(),
    (5.13094, 4.68175, 4.35225, 90, 101.476, 90))

  # Read and write map with offset

  print("\nReading and writing map with origin not at zero, use_mrcfile=",use_mrcfile)

  from scitbx.array_family.flex import grid
  from scitbx.array_family import flex
  from cctbx import uctbx, sgtbx

  real_map=m.map_data().as_double()
  grid_start=(5,5,5)
  grid_end=(9,10,11)

  if use_mrcfile:
    iotbx.mrcfile.write_ccp4_map(
      file_name="shifted_map.mrc",
      unit_cell=uctbx.unit_cell(m.unit_cell().parameters()),
      space_group=sgtbx.space_group_info("P1").group(),
      gridding_first=grid_start,
      gridding_last=grid_end,
      map_data=real_map,
      labels=flex.std_string(["iotbx.ccp4_map.tst"]),
      verbose=True)
    m = mrcfile.map_reader(file_name='shifted_map.mrc',verbose=True)

  else:
    iotbx.ccp4_map.write_ccp4_map(
      file_name="shifted_map.ccp4",
      unit_cell=uctbx.unit_cell(m.unit_cell().parameters()),
      space_group=sgtbx.space_group_info("P1").group(),
      gridding_first=grid_start,
      gridding_last=grid_end,
      map_data=real_map,
      labels=flex.std_string(["iotbx.ccp4_map.tst"]))
    m = iotbx.ccp4_map.map_reader(file_name='shifted_map.ccp4')

  print(m.map_data().origin(),m.map_data().all())
  print("GRID:",m.unit_cell_grid)
  assert m.unit_cell_grid == (16, 8, 16)
  assert approx_equal(m.unit_cell().parameters(), (
    82.095001220703125, 37.453998565673828, 69.636001586914062,
    90.0, 101.47599792480469, 90.0))
  assert m.unit_cell_crystal_symmetry().space_group_number()== 1
  print(m.map_data().origin(),m.map_data().all())
  assert m.map_data().origin() == (5,5,5)
  assert m.map_data().all() == (5,6,7)

def exercise_with_tst_input_map_2(use_mrcfile=None,file_name=None):
  if not file_name:
    file_name = libtbx.env.under_dist(
      module_name="iotbx",
      path="ccp4_map/ccp4_20_21_22_to_30_40_50.ccp4")

  print("\nTesting read of input map with offset origin and axis order 3 1 2 "+\
         "\n and use_mrcfile=",use_mrcfile)
  if use_mrcfile:
    m = mrcfile.map_reader(file_name=file_name,verbose=True)
  else:
    m = iotbx.ccp4_map.map_reader(file_name=file_name)
  print(m.unit_cell_grid,m.map_data().origin(),m.map_data().all(),m.unit_cell().parameters(), end=' ')
  assert m.unit_cell_grid == (55, 59, 61)
  assert approx_equal(m.unit_cell().parameters(),
    (24.53101921081543, 24.61971664428711, 26.457056045532227,
      90.0, 90.0, 90.0))
  assert m.map_data().origin() == (20, 21, 22)
  assert m.map_data().all() == (11, 20, 29)


def exercise_crystal_symmetry_from_ccp4_map(use_mrcfile=None):
  file_name = libtbx.env.under_dist(
    module_name="iotbx",
    path="ccp4_map/tst_input.map")
  if use_mrcfile:
    from iotbx.mrcfile import crystal_symmetry_from_ccp4_map
  else:
    from iotbx.ccp4_map import crystal_symmetry_from_ccp4_map
  cs = crystal_symmetry_from_ccp4_map.extract_from(file_name=file_name)
  print(cs.unit_cell().parameters(),cs.space_group().info())


def exercise(args,use_mrcfile=None):
  exercise_with_tst_input_map(use_mrcfile=use_mrcfile)
  exercise_with_tst_input_map_2(use_mrcfile=use_mrcfile)
  for file_name in args:
    print("\n",file_name,"use_mrcfile=",use_mrcfile)
    if use_mrcfile:
      m = iotbx.mrcfile.map_reader(file_name=file_name)
    else:
      m = iotbx.ccp4_map.map_reader(file_name=file_name)
    print("header_min: ", m.header_min)
    print("header_max: ", m.header_max)
    print("header_mean:", m.header_mean)
    print("header_rms: ", m.header_rms)
    print("unit cell grid:", m.unit_cell_grid)
    print("unit cell parameters:", m.unit_cell_parameters)
    print("space group number:  ", m.space_group_number)
    print("map origin:", m.map_data().origin())
    print("map grid:  ", m.map_data().all())
    map_stats = maptbx.statistics(m.data)
    assert approx_equal(map_stats.min(), m.header_min)
    assert approx_equal(map_stats.max(), m.header_max)
    assert approx_equal(map_stats.mean(), m.header_mean)
    if (m.header_rms != 0):
      assert approx_equal(map_stats.sigma(), m.header_rms)
    print()

def exercise_read_write_defaults(use_mrcfile=True):
  if not use_mrcfile: return

  from cctbx import uctbx, sgtbx, crystal
  from scitbx.array_family import flex
  mt = flex.mersenne_twister(0)
  nxyz = (4,5,6,)
  grid = flex.grid(nxyz)
  real_map_data = mt.random_double(size=grid.size_1d())
  real_map_data.reshape(grid)
  unit_cell=uctbx.unit_cell((10,10,10,90,90,90))
  space_group=sgtbx.space_group_info("P1").group()
  crystal_symmetry=crystal.symmetry(unit_cell=unit_cell,space_group=space_group)
  iotbx.mrcfile.write_ccp4_map(
    file_name="simple.mrc",
    crystal_symmetry=crystal_symmetry,
    map_data=real_map_data)
  input_real_map = iotbx.mrcfile.map_reader(file_name="simple.mrc")

  # check unit cell and space group:
  cs=input_real_map.crystal_symmetry()
  assert cs.is_similar_symmetry(crystal_symmetry)

  # Check writing map with offset
  origin_shift=(5,6,7)
  iotbx.mrcfile.write_ccp4_map(
    file_name="offset.mrc",
    crystal_symmetry=crystal_symmetry,
    origin_shift_grid_units=origin_shift,
    map_data=real_map_data)
  input_real_map = iotbx.mrcfile.map_reader(file_name="offset.mrc")
  map_data=input_real_map.map_data()
  assert map_data.origin()==origin_shift

def exercise_writer(use_mrcfile=None,output_axis_order=[3,2,1]):
  from cctbx import uctbx, sgtbx
  from scitbx.array_family import flex
  mt = flex.mersenne_twister(0)
  nxyz = (4,5,6,)
  grid = flex.grid(nxyz)
  real_map_data = mt.random_double(size=grid.size_1d())
  real_map_data.reshape(grid)
  unit_cell=uctbx.unit_cell((10,10,10,90,90,90))

  if use_mrcfile:
    from iotbx.mrcfile import create_output_labels
    labels=create_output_labels(
       program_name='test',
       limitations=['extract_unique','no_wrapping_outside_cell'],
       output_labels=['test label'],
       )
    iotbx.mrcfile.write_ccp4_map(
      file_name="four_five_six.mrc",
      unit_cell=unit_cell,
      space_group=sgtbx.space_group_info("P1").group(),
      map_data=real_map_data,
      labels=labels,
      output_axis_order=output_axis_order)
    input_real_map = iotbx.mrcfile.map_reader(file_name="four_five_six.mrc")
    for x in input_real_map.labels:
      print("LABEL: ",x)
    assert str(input_real_map.labels).find('extract_unique')>-1
  else:
    iotbx.ccp4_map.write_ccp4_map(
      file_name="four_five_six.map",
      unit_cell=unit_cell,
      space_group=sgtbx.space_group_info("P1").group(),
      map_data=real_map_data,
      labels=flex.std_string(["iotbx.ccp4_map.tst"]))
    input_real_map = iotbx.ccp4_map.map_reader(file_name="four_five_six.map")

  if use_mrcfile:
    assert not input_real_map.wrapping_from_input_file()


  input_map_data=input_real_map.map_data()
  real_map_mmm = real_map_data.as_1d().min_max_mean()
  input_map_mmm = input_map_data.as_double().as_1d().min_max_mean()
  cc=flex.linear_correlation(real_map_data.as_1d(),input_map_data.as_double().as_1d()).coefficient()
  assert cc > 0.999
  print("\nMRCFILE with 4x5x6 map and axis order %s %s" %(output_axis_order,cc))

  assert approx_equal(input_real_map.unit_cell().parameters(), unit_cell.parameters())
  assert approx_equal(real_map_mmm.min, input_real_map.header_min,eps=0.001)
  assert approx_equal(real_map_mmm.min, input_map_mmm.min,eps=0.001)


  # random small maps of different sizes
  for nxyz in flex.nested_loop((2,1,1),(4,4,4)):
    mt = flex.mersenne_twister(0)
    grid = flex.grid(nxyz)
    real_map = mt.random_double(size=grid.size_1d())
    real_map=real_map-0.5
    real_map.reshape(grid)
    if use_mrcfile:
      iotbx.mrcfile.write_ccp4_map(
        file_name="random.mrc",
        unit_cell=uctbx.unit_cell((1,1,1,90,90,90)),
        space_group=sgtbx.space_group_info("P1").group(),
        gridding_first=(0,0,0),
        gridding_last=tuple(grid.last(False)),
        map_data=real_map,
        labels=flex.std_string(["iotbx.ccp4_map.tst"]))
      m = iotbx.mrcfile.map_reader(file_name="random.mrc")
    else:
      iotbx.ccp4_map.write_ccp4_map(
        file_name="random.map",
        unit_cell=uctbx.unit_cell((1,1,1,90,90,90)),
        space_group=sgtbx.space_group_info("P1").group(),
        gridding_first=(0,0,0),
        gridding_last=tuple(grid.last(False)),
        map_data=real_map,
        labels=flex.std_string(["iotbx.ccp4_map.tst"]))
      m = iotbx.ccp4_map.map_reader(file_name="random.map")

    mmm = flex.double(list(real_map)).min_max_mean()
    m1=real_map.as_1d()
    m2=m.map_data().as_double().as_1d()
    cc=flex.linear_correlation(m1,m2).coefficient()
    assert cc > 0.999
    assert approx_equal(m.unit_cell().parameters(), (1,1,1,90,90,90))
    assert approx_equal(mmm.min, m.header_min)
    assert approx_equal(mmm.max, m.header_max)
    #
    # write unit_cell_grid explicitly to map
    iotbx.ccp4_map.write_ccp4_map(
      file_name="random_b.map",
      unit_cell=uctbx.unit_cell((1,1,1,90,90,90)),
      space_group=sgtbx.space_group_info("P1").group(),
      unit_cell_grid=real_map.all(),
      map_data=real_map,
      labels=flex.std_string(["iotbx.ccp4_map.tst"]))
    m = iotbx.ccp4_map.map_reader(file_name="random_b.map")
    m1=real_map.as_1d()
    m2=m.map_data().as_double().as_1d()
    cc=flex.linear_correlation(m1,m2).coefficient()
    assert cc > 0.999

    mmm = flex.double(list(real_map)).min_max_mean()
    assert approx_equal(m.unit_cell_parameters, (1,1,1,90,90,90))
    assert approx_equal(mmm.min, m.header_min)
    assert approx_equal(mmm.max, m.header_max)
    #

    #
    gridding_first = (0,0,0)
    gridding_last=tuple(grid.last(False))
    map_box = maptbx.copy(real_map, gridding_first, gridding_last)
    map_box.reshape(flex.grid(map_box.all()))

    if use_mrcfile:
      iotbx.mrcfile.write_ccp4_map(
        file_name="random_box.mrc",
        unit_cell=uctbx.unit_cell((1,1,1,90,90,90)),
        space_group=sgtbx.space_group_info("P1").group(),
        map_data=map_box,
        labels=flex.std_string(["iotbx.ccp4_map.tst"]))
    else:
      iotbx.ccp4_map.write_ccp4_map(
        file_name="random_box.map",
        unit_cell=uctbx.unit_cell((1,1,1,90,90,90)),
        space_group=sgtbx.space_group_info("P1").group(),
        map_data=map_box,
        labels=flex.std_string(["iotbx.ccp4_map.tst"]))
  print("OK")

def run(args):
  import iotbx.ccp4_map
  exercise_read_write_defaults()  # only for mrcfile
  for use_mrcfile in [True,False]:
    exercise(args=args,use_mrcfile=use_mrcfile)
    exercise_writer(use_mrcfile=use_mrcfile)
    if use_mrcfile:
      exercise_writer(use_mrcfile=use_mrcfile,output_axis_order=[1,2,3],)
      exercise_writer(use_mrcfile=use_mrcfile,output_axis_order=[2,3,1],)

    exercise_crystal_symmetry_from_ccp4_map(use_mrcfile=use_mrcfile)
  print(format_cpu_times())

if (__name__ == "__main__"):
  run(sys.argv[1:])
