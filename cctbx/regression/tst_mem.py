from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import iotbx.pdb
import cctbx.maptbx.mem as mem
from cctbx import miller
from six.moves import zip

def map_1d(xrs, map_data, n_steps=300, step_size=0.01):
  atom_center = xrs.scatterers()[0].site
  a_cell = xrs.unit_cell().parameters()[0]
  dist = flex.double()
  rho = flex.double()
  x=-10
  while x<=10:
    x += step_size
    point = x/a_cell, 0, 0
    density_value_at_point_e = map_data.eight_point_interpolation(point)
    dist.append(x)
    rho.append(density_value_at_point_e)
  return dist, rho

def scale(x, y):
  assert x.size() == y.size()
  x = flex.abs(x)
  y = flex.abs(y)
  return flex.sum(x*y)/flex.sum(y*y)

def r_factor(x, r1, r2, eps=1.e-6):
  sel  = x > 0-eps
  sel &= x < 0+eps
  assert sel.count(True) == 1
  r1_ = r1.select(sel)
  r2_ = r2.select(sel)
  return flex.abs( r1_-r2_ )[0]

def run():
  pdb_str="""
CRYST1    5.000    5.000    5.000  90.00  90.00  90.00 P 1
HETATM    1  C    C      1       0.000   0.000   0.000  1.00  5.00           C
END
"""
  xrs = iotbx.pdb.input(source_info=None, lines=pdb_str).xray_structure_simple()
  ###
  fc = xrs.structure_factors(d_min = 0.26, algorithm = "direct").f_calc()
  fft_map_ref = fc.fft_map(resolution_factor = 0.1)
  fft_map_ref.apply_sigma_scaling()
  map_data_ref = fft_map_ref.real_map_unpadded()
  r_r, rho_r = map_1d(xrs=xrs,
    map_data=map_data_ref.deep_copy()/flex.max(map_data_ref))
  F_0 = xrs.structure_factors(d_min=1.4999999, algorithm="direct").f_calc()
  #
  fft_map_0 = miller.fft_map(
    crystal_gridding     = fft_map_ref,
    fourier_coefficients = F_0)
  fft_map_0.apply_sigma_scaling()
  map_data_0 = fft_map_0.real_map_unpadded()
  r_s, rho_s = map_1d(xrs=xrs, map_data=map_data_0)
  sc = scale(rho_r, rho_s)
  rho_s = rho_s*sc
  #
  R = r_factor(x=r_s, r1=rho_r, r2=rho_s)
  cc = flex.linear_correlation(x = rho_r, y = rho_s).coefficient()
  print("Initial dist: %6.3f"%R, cc)
  ### Exercise fixed lam
  lam = 0.74
  for start_map in ["flat", "lde", "min_shifted"]:
    m = mem.run(f=F_0, f_000=xrs.f_000(), lam=lam, resolution_factor=0.1,
      verbose=False, start_map=start_map, max_iterations=1300,
      detect_convergence=False)
    r_m, rho_m = map_1d(xrs=xrs, map_data=m.rho)
    sc = scale(rho_r, rho_m)
    rho_m = rho_m*sc
    R = r_factor(x=r_s, r1=rho_r, r2=rho_m)
    cc = flex.linear_correlation(x = rho_r, y = rho_m).coefficient()
    print(m.show(verbose=False), "*", "dist: %6.3f"%R, cc)
    assert R < 0.0004, R
    assert cc > 0.999, R
  ### Exercise auto-incremented lam, 1
  m = mem.run(f=F_0, f_000=xrs.f_000(), lam=0.01, resolution_factor=0.1,
    lambda_increment_factor = 1.01, beta=0.5, verbose=True,
    start_map="min_shifted", max_iterations=800, detect_convergence=True)
  r_m, rho_m = map_1d(xrs=xrs, map_data=m.rho)
  sc = scale(rho_r, rho_m)
  rho_m = rho_m*sc
  R = r_factor(x=r_s, r1=rho_r, r2=rho_m)
  cc = flex.linear_correlation(x = rho_r, y = rho_m).coefficient()
  print(m.show(verbose=False), "*", "dist: %6.3f"%R, cc)
  assert R < 0.004, R
  assert cc > 0.999, R
  ### Exercise auto-incremented lam, 2
  m = mem.run(f=F_0, f_000=xrs.f_000(), lam=0.01, resolution_factor=0.1,
    lambda_increment_factor = 1.01, beta=0.5, xray_structure=xrs, verbose=True,
    start_map="min_shifted", max_iterations=800, detect_convergence=True)
  r_m, rho_m = map_1d(xrs=xrs, map_data=m.rho)
  sc = scale(rho_r, rho_m)
  rho_m = rho_m*sc
  R = r_factor(x=r_s, r1=rho_r, r2=rho_m)
  cc = flex.linear_correlation(x = rho_r, y = rho_m).coefficient()
  print(m.show(verbose=False), "*", "dist: %6.3f"%R, cc)
  assert R < 0.05, R
  assert cc > 0.998, R
  #
  of = open("gp","w")
  assert rho_r.size() == rho_s.size() == rho_m.size()
  for r, r1, r2, r3 in zip(r_r, rho_r, rho_s, rho_m):
    print("%15.12f %15.12f %15.12f %15.12f"%(r, r1, r2, r3), file=of)
  of.close()
  #
  m.write_mtz_file()
  #
  print("gnuplot command to show results:")
  print("""
plot [-3:3] [-0.1:1.1] "gp" using 1:2 with lines title "exact" lw 3, "gp" using 1:3 with lines title "synthesis" lw 3 lc rgb "black", "gp" using 1:4 with lines title "MEM" lw 3
""")

if __name__ == "__main__" :
  run()
