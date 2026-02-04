from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx.maptbx.bcr import qmap
from cctbx import maptbx, adptbx, xray, crystal
import time, math, sys
from libtbx.test_utils import approx_equal
from itertools import combinations
from cctbx.maptbx.bcr import bcr
from libtbx import group_args
import bisect
import libtbx

def get_xrs_and_cg(sites_cart, b_iso, occ, uc_params, n_real, table = "wk1995"):
  cs = crystal.symmetry(uc_params, "P 1")
  sp = crystal.special_position_settings(cs)
  scatterers = flex.xray_scatterer()
  for site_cart in sites_cart:
    site_frac = cs.unit_cell().fractionalize(site_cart)
    scatterer = xray.scatterer(
      "c", site=site_frac, u=adptbx.b_as_u(b_iso), occupancy=occ)
    scatterers.append(scatterer)
  xrs = xray.structure(sp, scatterers)
  #
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell        = cs.unit_cell(),
    space_group_info = cs.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    pre_determined_n_real = n_real,
    )
  xrs.scattering_type_registry(
    table = table,
    d_min = 1,
    types_without_a_scattering_contribution=["?"])
  return xrs, crystal_gridding

def get_fft_map(xrs, cg, d_min):
  f000 = xrs.f_000()
  fc = xrs.structure_factors(
    d_min          = d_min,
    algorithm      = "direct",
    cos_sin_table  = False,
    anomalous_flag = False).f_calc()
  fft_map = fc.fft_map(
    crystal_gridding = cg,
    f_000            = xrs.f_000()
    )
  fft_map.apply_volume_scaling()
  m1 = fft_map.real_map_unpadded()
  return fc.generate_bijvoet_mates(), m1

def get_vrm_map(xrs, cg, d_min, RadFact, RadAdd):
  o = qmap.compute(
    xray_structure = xrs,
    n_real         = cg.n_real(),
    resolution     = d_min,
    RadFact        = RadFact,
    RadAdd         = RadAdd,
    resolutions    = None,
    use_exp_table  = False, # XXX <<<<
    debug          = False,
    show_BCR       = False)
  m2 = o.map_data()
  return m2

def get_formula_map(fc, mask):
  nx,ny,nz = mask.all()
  m3 = flex.double(flex.grid(mask.all()), 0)
  sx,sy,sz = 1./nx, 1./ny, 1./nz
  for i in range(nx):
    for j in range(ny):
      for k in range(nz):
        if mask[i,j,k] > 0.1:
          xf,yf,zf = sx*i,sy*j,sz*k
          mv = fc.direct_summation_at_point(site_frac=[xf,yf,zf]).real
          m3[i,j,k] = mv
  return m3

def get_formula12_or_89_map(n_real, xrs, r_image, radii=None, image=None, approx=None,
                            m12c=None, m12f=None, m89c=None, m89f=None):
  assert [m12c, m12f, m89c, m89f].count(True) == 1
  o = qmap.load_table(element="C", table="wk1995")
  oim = maptbx.atom_curves(scattering_type="C", scattering_table="wk1995")
  d = o["1.0"] ### XXX HARD WIRED VALUE
  B = d["B"]
  C = d["C"]
  R = d["R"]
  nx,ny,nz = n_real
  a,b,c = xrs.unit_cell().parameters()[:3]
  result = flex.double(flex.grid(n_real), 0)
  sx,sy,sz = a/nx,b/ny,c/nz
  mx,my,mz = round(r_image/sx), round(r_image/sy), round(r_image/sz)
  for site_cart in xrs.sites_cart():
    xc,yc,zc = site_cart
    ic,jc,kc = round(xc/sx), round(yc/sy), round(zc/sz)
    for ii in range(-mx, mx):
      for jj in range(-my, my):
        for kk in range(-mz, mz):
          i,j,k = ic+ii,jc+jj,kc+kk
          x,y,z = sx*i, sy*j, sz*k
          r = math.sqrt((xc-x)**2 + (yc-y)**2 + (zc-z)**2)
          if r < 0.001: r = 0.0
          if r > r_image: continue
          if m12c:
            mv = interpolate(x_values=radii, y_values=approx, x=r)
          elif m12f:
            mv = bcr.curve(B=B, C=C, R=R, radii=[r,], b_iso=0)[0]
          elif m89c:
            mv = interpolate(x_values=radii, y_values=image, x=r)
          elif m89f:
            mv = oim.image(d_min=1,b_iso=0,radii=flex.double([r,]),
              compute_derivatives=False, fast=False).image_values[0]
          else: assert 0
          result[i,j,k] += mv
  return result

def get_mask(xrs, n_real, rad):
  return maptbx.mask(
      xray_structure              = xrs,
      n_real                      = n_real,
      mask_value_inside_molecule  = 1,
      mask_value_outside_molecule = 0,
      solvent_radius              = 0,
      atom_radius                 = rad,
      wrapping                    = False)

def get_errors(map1, map2, rho_image_0, inside=None, mask=None):
  if mask is not None:
    assert inside in [True,False]
    s = mask.as_1d()>0.1
    if not inside: s = ~s
    m1, m2 = map1.as_1d(),  map2.as_1d()
    m1, m2 = m1.select(s), m2.select(s)
  else:
    m1 = map1.deep_copy()
    m2 = map2.deep_copy()
  #
  diff       = flex.abs(m1 - m2)
  e_abs_mean = flex.mean( diff  )
  e_abs_max  = flex.max(  diff  )
  e_rel_mean = flex.mean( diff / rho_image_0 )
  e_rel_max  = flex.max(  diff / rho_image_0 )
  string = "e_abs_mean: %6.4f e_abs_max: %6.4f e_rel_max: %6.4f"%(
    e_abs_mean,e_abs_max,   e_rel_max)
  return group_args(string  = string)

def interpolate(x_values, y_values, x):
    """
    Interpolates to find the corresponding y value for a given x using linear
    interpolation.

    Parameters:
    - x_values: List of x values (grid points).
    - y_values: List of corresponding y values (function values).
    - x: The x value for which to find the corresponding y value.

    Returns:
    - The interpolated y value corresponding to x.
    """
    # Check if x_values are ordered
    #if any(x_values[i] > x_values[i + 1] for i in range(len(x_values) - 1)):
    #  raise ValueError("x_values must be in ascending order.")
    # Check if x is outside the interpolation bounds
    if x < x_values[0] or x > x_values[-1]:
      raise ValueError("x is out of interpolation bounds:", x)
    # Use binary search to find the right interval
    i = bisect.bisect_right(x_values, x) - 1  # Get index of the right interval
    # If x is exactly the last point in x_values
    if i == len(x_values) - 1:
      return y_values[i]  # Return the corresponding y directly
    # Perform linear interpolation
    x0, x1 = x_values[i], x_values[i + 1]
    y0, y1 = y_values[i], y_values[i + 1]
    # Linear interpolation formula
    y = y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    return y

def get_image_AU():
  image_AU = flex.double()
  radii_AU = flex.double()
  path=libtbx.env.find_in_repositories("cctbx/maptbx/bcr")
  with open(path+"/CNOPSFeAu-6_1A_dist21.tab", "r") as fo:
    start_reading = False
    for l in fo.readlines():
      if '     N     dist       C              N' in l:
        start_reading = True
        continue
      if start_reading:
        l=l.strip().split()
        radii_AU.append(float(l[1]))
        image_AU.append(float(l[2]))
  return radii_AU, image_AU

def get_image():
  #ff_S = [6.372157,  5.154568, 1.473732,  1.635073,  1.209372, 0.154722,
  #        1.514347, 22.092527, 0.061373, 55.445175,  0.646925, 0.000000]
  ff_C  = [2.657506, 1.078079, 1.490909,  -4.24107,  0.713791, 4.297983,
          14.780758, 0.776775, 42.086842, -0.000294, 0.239535, 0.000000]
  image, radii = maptbx.atom_image_fast(
    ff_packed = ff_C,
    d_min     = 1.0,
    n_grid    = 11000,
    dist_max  = 11.0)
  return image, radii

def get_approx(radii):
  t = qmap.load_table(element="C", table="wk1995")
  d = t["1.0"]
  B = d["B"]
  C = d["C"]
  R = d["R"]
  return bcr.curve(B=B, C=C, R=R, radii=radii, b_iso=0)

def run(d_min=1.0, r_image = 3.0, r_atom = 2.5):
  #
  # READ IMAGE FROM CNOPSFeAu-6_1A_dist21.tab
  #
  radii_AU, image_AU = get_image_AU()
  #
  # COMPUTE IMAGE WITH CCTBX, FORMULAS (8-9)
  #
  image, radii = get_image()
  rho_image_0 = image[0]
  #
  # ASSERT BOTH IMAGES ARE IDENTICAL
  #
  assert image_AU.size() == 11001
  assert approx_equal(image_AU, image)
  assert approx_equal(radii_AU, radii)
  #
  # COMPUTE BCR APPROXIMATION, FORMULAS (1-2)
  #
  approx = get_approx(radii = radii)
  assert approx.size() == image.size()
  #
  # WRITE ALL CURVES
  #
  #file_name = "S_images_approx_1A.dat"
  #with open(file_name, "w") as fo:
  #  print(" radii    imags(AU)     image(cctbx)  approx", file=fo)
  #  for r, i1, i2, a in zip(radii, image_AU, image, approx):
  #    print("%7.4f %13.8f %13.8f %13.8f"%(r, i1, i2, a), file=fo)
  #print()
  #eo = get_errors(map1=image, map2=approx)
  #print("2D image vs approx (file: %s): "%file_name, eo.string)
  #
  # 3D
  #
  uc_params = (64, 72, 80, 90, 90, 90)
  n_real    = [128, 144, 160]
  #
  print()
  print("unit cell:", uc_params)
  print("NXYZ:", n_real)
  print()
  # TWO-ATOM TEST
  two0   = [[32,36,40], [32,   36,   40],]
  two1   = [[32,36,40], [32+1, 36+1, 40+1],]
  two2   = [[32,36,40], [32+2, 36+2, 40+2],]
  two3   = [[32,36,40], [32+3, 36+3, 40+3],]
  two4   = [[32,36,40], [32+4, 36+4, 40+4],]
  two5   = [[32,36,40], [32+5, 36+5, 40+5],]
  cases = [two0, two1, two2, two3, two4, two5]
  #
  for m, it in enumerate(cases):
    sites_cart, b, occ = it, 0, 1
    xrs, cg = get_xrs_and_cg(
      sites_cart = sites_cart,
      b_iso      = b,
      occ        = occ,
      uc_params  = uc_params,
      n_real     = n_real)
    mask = get_mask(xrs=xrs, n_real=cg.n_real(), rad=r_atom)

    fc, m_fft = get_fft_map(xrs=xrs, cg=cg, d_min=d_min)
    #m_vrm     = get_vrm_map(xrs=xrs, cg=cg, d_min=d_min, RadFact=3.0, RadAdd=1.5)
    #m_12c     = get_formula12_or_89_map(n_real=mask.all(), xrs=xrs,
    #            r_image=r_image, radii=radii, approx=approx, m12c=True)
    #m_12f     = get_formula12_or_89_map(n_real=mask.all(), xrs=xrs,
    #            r_image=r_image, m12f=True)
    m_89c     = get_formula12_or_89_map(n_real=mask.all(), xrs=xrs,
                r_image=r_image, radii=radii, image=image, m89c=True)
    #m_89f     = get_formula12_or_89_map(n_real=mask.all(), xrs=xrs, r_image=r_image, m89f=True)
    #m_dir     = get_formula_map(fc=fc, mask=mask11)

    options = [["fft", m_fft],
               #["vrm", m_vrm],
               #["12f", m_12f],
               #["89f", m_89f],
               #["12c", m_12c],
               ["89c", m_89c],
               #["dir", m_dir]
               ]

    for pair in list(combinations(options, 2)):
      one, two = pair
      eo_in  = get_errors(map1=one[1], map2=two[1], rho_image_0=rho_image_0, mask=mask, inside=True)
      eo_out = get_errors(map1=one[1], map2=two[1], rho_image_0=rho_image_0, mask=mask, inside=False)
      print("m: %2d: %s vs %s: in: %s out: %s"%(m, one[0], two[0], eo_in.string, eo_out.string))
      sys.stdout.flush()
    #

if (__name__ == "__main__"):
  start = time.perf_counter()
  for r_image in [3, 5, 10]:
    print("rmax:", r_image)
    run(r_image = r_image)
    print()
  print("Time:", time.perf_counter()-start)
