from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx import miller
from cctbx import crystal
from cctbx import maptbx
from libtbx.test_utils import approx_equal
import boost_adaptbx.boost.python as bp
import math, time, random
asu_map_ext = bp.import_ext("cctbx_asymmetric_map_ext")

def relative_difference(a, b):
  return flex.sum(flex.abs(a-b))/max(1.e-30, flex.sum(flex.abs(b)))

def points_and_mask(crystal_symmetry, n_real, n_seeds, radius=2):
  """
  Random blobs of grid points, expanded by space group symmetry. Returns the
  sorted 1d grid indices and the 0/1 map.
  """
  nx,ny,nz = n_real
  ops = crystal_symmetry.space_group().all_ops()
  selected = set()
  for i_seed in range(n_seeds):
    c = [random.randrange(n) for n in n_real]
    for dx in range(-radius, radius+1):
      for dy in range(-radius, radius+1):
        for dz in range(-radius, radius+1):
          if(dx*dx+dy*dy+dz*dz > radius*radius): continue
          site = ((c[0]+dx)/nx, (c[1]+dy)/ny, (c[2]+dz)/nz)
          for op in ops:
            s = op * site
            j = [s[0]*nx, s[1]*ny, s[2]*nz]
            ji = [int(round(v)) for v in j]
            assert max([abs(a-b) for a,b in zip(j,ji)]) < 1.e-6 # grid fits
            selected.add(((ji[0]%nx)*ny + ji[1]%ny)*nz + ji[2]%nz)
  points = flex.size_t(sorted(selected))
  mask = flex.double(flex.grid(n_real), 0)
  mask.as_1d().set_selected(points, 1)
  return points, mask

def exercise_vs_fft(space_group_symbol, unit_cell, d_min=4, step=0.6):
  """
  LoFT equals the Fourier transform of the region mask, for unique and for
  anomalous sets of indices, with and without space group symmetry.
  """
  cs = crystal.symmetry(unit_cell, space_group_symbol)
  cg = maptbx.crystal_gridding(
    unit_cell        = cs.unit_cell(),
    space_group_info = cs.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = step)
  n_real = cg.n_real()
  points, mask = points_and_mask(crystal_symmetry=cs, n_real=n_real,
    n_seeds=5)
  assert points.size() > 100
  n_grid = n_real[0]*n_real[1]*n_real[2]
  v = cs.unit_cell().volume()/n_grid
  asu = asu_map_ext.asymmetric_map(cs.space_group().type(), mask).data()
  for anomalous_flag in [False, True]:
    ms = miller.build_set(crystal_symmetry=cs, anomalous_flag=anomalous_flag,
      d_min=d_min)
    ms = ms.array(data = flex.complex_double(ms.indices().size()))
    f_fft = ms.structure_factors_from_asu_map(
      asu_map_data = asu, n_real = n_real).data()
    f_loft = maptbx.loft(
      points           = points,
      n_real           = n_real,
      miller_indices   = ms.indices(),
      volume_per_point = v)
    assert f_loft.size() == ms.indices().size()
    assert relative_difference(f_loft, f_fft) < 1.e-9, [space_group_symbol,
      anomalous_flag, relative_difference(f_loft, f_fft)]
  # F(000) is the region volume
  f000 = maptbx.loft(points=points, n_real=n_real,
    miller_indices=flex.miller_index([(0,0,0)]), volume_per_point=v)[0]
  assert approx_equal(f000.real, points.size()*v, 1.e-9)
  assert approx_equal(f000.imag, 0, 1.e-9)
  # Order of the points does not matter
  ms = miller.build_set(crystal_symmetry=cs, anomalous_flag=True, d_min=d_min)
  perm = flex.random_permutation(points.size())
  f1 = maptbx.loft(points=points, n_real=n_real,
    miller_indices=ms.indices(), volume_per_point=v)
  f2 = maptbx.loft(points=points.select(perm), n_real=n_real,
    miller_indices=ms.indices(), volume_per_point=v)
  assert relative_difference(f2, f1) < 1.e-12
  # P1 expansion of the unique indices (a hemisphere, mixed signs) against
  # a plain P1 FFT of the real map
  ms_p1 = miller.build_set(crystal_symmetry=cs, anomalous_flag=False,
    d_min=d_min).expand_to_p1()
  ms_p1 = ms_p1.array(data = flex.complex_double(ms_p1.indices().size()))
  f_fft = ms_p1.structure_factors_from_map(map=mask, use_scale=True).data()
  f_loft = maptbx.loft(points=points, n_real=n_real,
    miller_indices=ms_p1.indices(), volume_per_point=v)
  assert relative_difference(f_loft, f_fft) < 1.e-9

def exercise_trivial():
  n_real = (8,10,12)
  hkl = flex.miller_index([(0,0,0), (1,-2,3), (-5,7,-11), (100,-200,300)])
  # no points
  f = maptbx.loft(points=flex.size_t(), n_real=n_real, miller_indices=hkl,
    volume_per_point=1)
  assert f.all_eq(0)
  # one point: a pure phase factor, also for indices beyond the grid size
  jx,jy,jz = 3,7,5
  p = flex.size_t([(jx*n_real[1] + jy)*n_real[2] + jz])
  f = maptbx.loft(points=p, n_real=n_real, miller_indices=hkl,
    volume_per_point=2.5)
  for fi, (h,k,l) in zip(f, hkl):
    arg = 2*math.pi*(h*jx/n_real[0] + k*jy/n_real[1] + l*jz/n_real[2])
    assert approx_equal(fi, 2.5*complex(math.cos(arg), math.sin(arg)), 1.e-9)
  # no reflections
  f = maptbx.loft(points=p, n_real=n_real,
    miller_indices=flex.miller_index(), volume_per_point=1)
  assert f.size() == 0
  # point outside the grid
  try:
    maptbx.loft(points=flex.size_t([8*10*12]), n_real=n_real,
      miller_indices=hkl, volume_per_point=1)
  except RuntimeError: pass
  else: raise RuntimeError("Exception expected.")

def exercise_sphere(d_min=3, R=3):
  """
  Spherical region: agreement with the analytical transform of a sphere,
  and time against FFT of the whole cell.
  """
  cs = crystal.symmetry((60, 72, 96, 90, 90, 90), "P1")
  a,b,c = cs.unit_cell().parameters()[:3]
  center = (11.0, 13.0, 15.0)
  n_real = (int(a/0.6), int(b/0.6), int(c/0.6))
  nx,ny,nz = n_real
  sx,sy,sz = a/nx, b/ny, c/nz
  selected = []
  m = int(R/min(sx,sy,sz))+2
  ic = [int(round(center[0]/sx)), int(round(center[1]/sy)),
        int(round(center[2]/sz))]
  for i in range(ic[0]-m, ic[0]+m+1):
    for j in range(ic[1]-m, ic[1]+m+1):
      for k in range(ic[2]-m, ic[2]+m+1):
        d2 = (center[0]-i*sx)**2 + (center[1]-j*sy)**2 + (center[2]-k*sz)**2
        if(d2 <= R*R): selected.append(((i%nx)*ny + j%ny)*nz + k%nz)
  points = flex.size_t(sorted(selected))
  mask = flex.double(flex.grid(n_real), 0)
  mask.as_1d().set_selected(points, 1)
  ms = miller.build_set(crystal_symmetry=cs, anomalous_flag=False, d_min=d_min)
  ms = ms.array(data = flex.complex_double(ms.indices().size()))
  v = cs.unit_cell().volume()/(nx*ny*nz)
  t0 = time.time()
  f_loft = ms.array(data = maptbx.loft(points=points, n_real=n_real,
    miller_indices=ms.indices(), volume_per_point=v))
  t_loft = time.time()-t0
  t0 = time.time()
  f_fft = ms.structure_factors_from_map(map=mask, use_scale=True)
  t_fft = time.time()-t0
  print("points: %d of %d, reflections: %d"%(points.size(), nx*ny*nz,
    ms.indices().size()))
  print("time LoFT: %.4f FFT: %.4f ratio FFT/LoFT: %.1f"%(t_loft, t_fft,
    t_fft/max(t_loft, 1.e-9)))
  assert relative_difference(f_loft.data(), f_fft.data()) < 1.e-9
  # analytical
  xf = cs.unit_cell().fractionalize(center)
  s = 1./ms.d_spacings().data()
  f_an = flex.complex_double()
  vol = 4*math.pi*R**3/3
  for hkl, si in zip(ms.indices(), s):
    arg = 2*math.pi*si*R
    g = vol*3*(math.sin(arg) - arg*math.cos(arg))/arg**3
    ph = 2*math.pi*(hkl[0]*xf[0] + hkl[1]*xf[1] + hkl[2]*xf[2])
    f_an.append(g*complex(math.cos(ph), math.sin(ph)))
  f_an = ms.array(data = f_an)
  cc = f_loft.map_correlation(other = f_an)
  print("map CC with the analytical sphere transform: %.4f"%cc)
  assert cc > 0.99

def exercise_centring(symbols=None, verbose=False):
  """
  With a space group given, only one of the points related by lattice
  centring translations is summed over (all of them contribute the same).
  The result is the same as without the space group, for all 230 space
  groups and some other settings: asymmetric-unit and anomalous sets of
  indices, and indices that include systematic absences.
  """
  from cctbx import sgtbx
  if(symbols is None):
    symbols = [sgtbx.space_group_info(number = i).type().lookup_symbol()
               for i in range(1, 231)]
    symbols += ["H 3", "R 3 :R", "A 1 2 1", "I 1 2 1", "B 2 21 2",
                "F 2 2 2", "I 41/a :1", "C 1 2/c 1", "H 3 2", "P 1 1 21",
                "F d -3 m :1", "I 41 2 2", "A b a 2"]
  n_centred = 0
  time_plain, time_sg = 0, 0
  for symbol in symbols:
    sgi = sgtbx.space_group_info(symbol = symbol)
    sg = sgi.group()
    cs = crystal.symmetry(
      unit_cell        = sgi.any_compatible_unit_cell(volume = 20000),
      space_group_info = sgi)
    cg = maptbx.crystal_gridding(
      unit_cell        = cs.unit_cell(),
      space_group_info = sgi,
      symmetry_flags   = maptbx.use_space_group_symmetry,
      step             = 1.0)
    n_real = cg.n_real()
    points, mask = points_and_mask(crystal_symmetry = cs, n_real = n_real,
      n_seeds = 3)
    v = cs.unit_cell().volume()/(n_real[0]*n_real[1]*n_real[2])
    if(sg.n_ltr() > 1): n_centred += 1
    sets = [
      miller.build_set(crystal_symmetry = cs, anomalous_flag = False,
        d_min = 3.5).indices(),
      miller.build_set(crystal_symmetry = cs, anomalous_flag = True,
        d_min = 3.5).indices(),
      # all indices of the lattice, systematic absences included
      miller.build_set(crystal_symmetry = crystal.symmetry(
        unit_cell = cs.unit_cell(), space_group_symbol = "P 1"),
        anomalous_flag = True, d_min = 5).indices()]
    for i_set, indices in enumerate(sets):
      assert indices.size() > 10, [symbol, indices.size()]
      t0 = time.time()
      f_plain = maptbx.loft(points = points, n_real = n_real,
        miller_indices = indices, volume_per_point = v)
      time_plain += time.time()-t0
      t0 = time.time()
      f_sg = maptbx.loft(points = points, n_real = n_real,
        miller_indices = indices, volume_per_point = v, space_group = sg)
      time_sg += time.time()-t0
      scale = flex.max(flex.abs(f_plain))
      assert scale > 0
      diff = flex.max(flex.abs(f_sg - f_plain))/scale
      assert diff < 1.e-10, [symbol, i_set, diff]
      if(i_set == 2 and sg.n_ltr() > 1):
        absent = sg.is_sys_absent(indices)
        # reflections absent because of the centring are exactly zero
        centring_only = sgtbx.space_group()
        for t in list(sg.ltr())[1:]: centring_only.expand_ltr(t)
        absent_c = centring_only.is_sys_absent(indices)
        assert absent_c.count(True) > 0
        assert flex.abs(f_sg.select(absent_c)).all_eq(0)
        assert flex.max(flex.abs(f_plain.select(absent_c)))/scale < 1.e-10
    # order of points does not matter
    perm = flex.random_permutation(points.size())
    f2 = maptbx.loft(points = points.select(perm), n_real = n_real,
      miller_indices = sets[0], volume_per_point = v, space_group = sg)
    f1 = maptbx.loft(points = points, n_real = n_real,
      miller_indices = sets[0], volume_per_point = v, space_group = sg)
    assert relative_difference(f2, f1) < 1.e-12
    # a set of points that does not obey the centring is rejected
    if(sg.n_ltr() > 1):
      try:
        maptbx.loft(points = points[1:], n_real = n_real,
          miller_indices = sets[0], volume_per_point = v, space_group = sg)
      except RuntimeError as e: assert "centring" in str(e), str(e)
      else: raise RuntimeError("Exception expected.")
    if(verbose): print(symbol, sg.n_ltr(), points.size(), n_real)
  if(len(symbols) > 200): assert n_centred > 80, n_centred
  print("space groups: %d (centred: %d), time without/with space group: "
    "%.2f %.2f"%(len(symbols), n_centred, time_plain, time_sg))
  # no points, no reflections
  sg = sgtbx.space_group_info("C 2 2 21").group()
  assert maptbx.loft(points = flex.size_t(), n_real = (8,8,8),
    miller_indices = flex.miller_index([(1,1,0)]), volume_per_point = 1,
    space_group = sg).all_eq(0)
  assert maptbx.loft(points = flex.size_t([0, 4*8*8+4*8]), n_real = (8,8,8),
    miller_indices = flex.miller_index(), volume_per_point = 1,
    space_group = sg).size() == 0

def run():
  random.seed(0)
  flex.set_random_seed(0)
  exercise_trivial()
  exercise_vs_fft("P 1",       (21.3, 25.7, 31.1, 83, 97, 105))
  exercise_vs_fft("C 1 2 1",   (30.5, 26.2, 28.9, 90, 103, 90))
  exercise_vs_fft("P 21 21 21",(24.4, 31.9, 27.3, 90, 90, 90))
  exercise_vs_fft("P 63 2 2",  (33.1, 33.1, 41.7, 90, 90, 120))
  exercise_sphere()
  exercise_centring()

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %.2f"%(time.time()-t0))
  print("OK")
