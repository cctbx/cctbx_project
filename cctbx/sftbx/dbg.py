import sys, os, math
from cctbx.misc import python_utils
from cctbx_boost.arraytbx import shared
from cctbx_boost import sgtbx
from cctbx_boost import adptbx
from cctbx_boost import sftbx
from cctbx_boost import fftbx
from cctbx import xutils
from cctbx.development import debug_utils
from cctbx.development import make_cns_input
from cctbx.development import run_shelx

def add_u_extra(xtal, u_extra):
  xtal_mod = xtal.copy_attributes()
  for site in xtal:
    site.set_Uiso(site.Uiso() + u_extra)
    xtal_mod.add_site(site)
  return xtal_mod

def zero_out_fpfdp(xtal):
  new_sites = shared.XrayScatterer()
  for site in xtal.Sites:
    site.set_fpfdp(0j)
    new_sites.append(site)
  xtal.Sites = new_sites

def show_joined_sets(h1, h2, js):
  for i,j in js.pairs():
    print h1[i], h2[j]
  print "singles 1:"
  for i in js.singles(0):
    print h1[i]
  print "singles 2:"
  for j in js.singles(1):
    print h2[j]

def exercise(SgInfo,
             d_min=3.,
             grid_resolution_factor = 1./3,
             adp=0,
             force_complex=0,
             random_f_prime=0,
             random_f_double_prime=0,
             use_cns=0,
             use_shelx=0):
  elements = ("N", "C", "C", "O", "N", "C", "C", "O")
  if (random_f_prime): random_f_prime = grid_resolution_factor * d_min
  xtal = debug_utils.random_structure(
    SgInfo, elements,
    volume_per_atom=50.,
    min_distance=1.5,
    general_positions_only=0,
    random_f_prime_d_min=random_f_prime,
    random_f_double_prime=random_f_double_prime,
    anisotropic_displacement_parameters=adp)
  if (0):
    zero_out_fpfdp(xtal)
  if (0):
    i_select = 3
    elements = [elements[i_select]]
    site = xtal.Sites[i_select]
    site.set_Coordinates((0,0,0))
    site.set_Uiso(0)
    print site.CAASF().Label(), site.w()
    new_sites = shared.XrayScatterer()
    new_sites.append(site)
    xtal.Sites = new_sites
  if (0):
    assert SgInfo.SgNumber() == 1
    from cctbx_boost import uctbx
    xtal.UnitCell = uctbx.UnitCell((10,10,10))
  print xtal.UnitCell
  debug_utils.print_sites(xtal)
  friedel_flag = not (random_f_double_prime or force_complex)
  print "friedel_flag:", friedel_flag
  MillerIndices = xutils.build_miller_indices(xtal, friedel_flag, d_min)
  Fcalc = xutils.calculate_structure_factors(MillerIndices, xtal)
  if (0):
    # reset fpfdp to verify that correlation is significantly different from 1.
    zero_out_fpfdp(xtal)
  if (use_cns):
    run_cns(elements, xtal, d_min, grid_resolution_factor, friedel_flag, Fcalc)
    return
  if (use_shelx):
    run_shelx.run_shelx(xtal.SgInfo.BuildLookupSymbol(), xtal, Fcalc)
    return
  max_q = 1. / (d_min**2)
  max_prime = 5
  mandatory_grid_factors = xtal.SgOps.refine_gridding()
  grid_logical = sftbx.determine_grid(
    xtal.UnitCell,
    max_q, grid_resolution_factor, max_prime, mandatory_grid_factors)
  rfft = fftbx.real_to_complex_3d(grid_logical)
  quality_factor = 100
  u_extra = sftbx.calc_u_extra(max_q, grid_resolution_factor, quality_factor)
  if (0):
    u_extra = adptbx.B_as_U(20)
  wing_cutoff = 1.e-3
  exp_table_one_over_step_size = -100
  electron_density_must_be_positive = 1
  sampled_density = sftbx.sampled_model_density(
    xtal.UnitCell, xtal.Sites,
    rfft.Nreal(), rfft.Mreal(),
    u_extra, wing_cutoff, exp_table_one_over_step_size,
    force_complex, electron_density_must_be_positive)
  assert friedel_flag == sampled_density.friedel_flag()
  print "u_extra:", sampled_density.u_extra(),
  print "b_extra:", adptbx.U_as_B(sampled_density.u_extra())
  if (1):
    print "number of passed scatterers:", \
      sampled_density.n_passed_scatterers()
    print "number of contributing scatterers:", \
      sampled_density.n_contributing_scatterers()
    print "number of anomalous scatterers:", \
      sampled_density.n_anomalous_scatterers()
    print "wing_cutoff:", sampled_density.wing_cutoff()
    print "exp_table_one_over_step_size:", \
      sampled_density.exp_table_one_over_step_size()
    print "exp_table_size:", sampled_density.exp_table_size()
    print "max_shell_radii:", sampled_density.max_shell_radii(),
    print "(%.4f, %.4f, %.4f)" % sampled_density.max_shell_radii_frac()
  tags = sftbx.grid_tags(rfft.Nreal())
  sym_flags = sftbx.map_symmetry_flags(1)
  tags.build(xtal.SgInfo, sym_flags)
  sampled_density.apply_symmetry(tags)
  if (friedel_flag):
    map = sampled_density.map_real_as_shared()
    map_stats = shared.statistics(map)
    if (0):
      print "Electron density"
      print "max %.6g" % (map_stats.max())
      print "min %.6g" % (map_stats.min())
      print "mean %.6g" % (map_stats.mean())
      print "sigma %.6g" % (map_stats.sigma())
  else:
    map = sampled_density.map_complex_as_shared()
  if (0):
    map = sftbx.structure_factor_map(
      xtal.SgOps, Fcalc.H, Fcalc.F, rfft.Ncomplex())
    rfft.backward(map)
    map_stats = shared.statistics(map)
    print "True electron density"
    print "max %.6g" % (map_stats.max())
    print "min %.6g" % (map_stats.min())
    print "mean %.6g" % (map_stats.mean())
    print "sigma %.6g" % (map_stats.sigma())
  if (friedel_flag):
    rfft.forward(map)
    map_stats = shared.statistics(map)
    if (0):
      print "Transformed electron density"
      print "max %.6g" % (map_stats.max())
      print "min %.6g" % (map_stats.min())
      print "mean %.6g" % (map_stats.mean())
      print "sigma %.6g" % (map_stats.sigma())
      print "Ncomplex", rfft.Ncomplex()
  else:
    cfft = fftbx.complex_to_complex_3d(rfft.Nreal())
    cfft.forward(map)
  if (0):
    map = sftbx.structure_factor_map(
      xtal.SgOps, Fcalc.H, Fcalc.F, rfft.Ncomplex())
  if (friedel_flag):
    miller_indices, fcal = sftbx.collect_structure_factors(
      xtal.UnitCell, xtal.SgInfo,
      max_q, map, rfft.Ncomplex())
  else:
    miller_indices, fcal = sftbx.collect_structure_factors(
      xtal.UnitCell, xtal.SgInfo,
      max_q, map, cfft.N())
  sampled_density.eliminate_u_extra_and_normalize(miller_indices, fcal)
  if (0):
    u_extra = sampled_density.u_extra()
    xtal_extra = add_u_extra(xtal, u_extra)
    fcalc_extra = xutils.calculate_structure_factors(MillerIndices, xtal_extra)
    show_structure_factor_correlation(
      "before", Fcalc.H, 0, Fcalc.F, fcalc_extra.F)
    sftbx.eliminate_u_extra(
      xtal.UnitCell, u_extra, MillerIndices.H, fcalc_extra.F)
    show_structure_factor_correlation(
      "after", Fcalc.H, 0, Fcalc.F, fcalc_extra.F)
  js = shared.join_sets(MillerIndices.H, miller_indices)
  if (0):
    show_joined_sets(MillerIndices.H, miller_indices, js)
  assert js.pairs().size() + js.singles(0).size() == MillerIndices.H.size()
  assert js.pairs().size() + js.singles(1).size() == miller_indices.size()
  assert js.pairs().size() == MillerIndices.H.size()
  for i in xrange(2):
    assert js.singles(i).size() == 0
  show_structure_factor_correlation(
    "sgtbx_dir/sgtbx_fft", Fcalc.H, js, Fcalc.F, fcal,
    min_corr_ampl=0.99, max_mean_w_phase_error=3.,
    verbose=0)

def write_cns_input(elements, xtal, d_min, grid_resolution_factor,
                    friedel_flag):
  cns_input = make_cns_input.topology(xtal.Sites)
  cns_input += make_cns_input.coordinates(xtal.Sites)
  cns_input += make_cns_input.unit_cell(xtal.UnitCell)
  cns_input += make_cns_input.symmetry(xtal.SgOps)
  l = cns_input.append
  l("""write coordinates end
coordinates orthogonalize end

xray
  fft
    grid=%.12g
  end
end"""
% (grid_resolution_factor,))
  l("""
xray
  anomalous=%s
  generate 100000. %.12g
  mapresolution %.12g
end"""
% (("false", "true")[friedel_flag==0], d_min, d_min))
  cns_input += make_cns_input.predict("f_dir", "direct")
  cns_input += make_cns_input.predict("f_fft", "fft")
  l("""xray
write reflections output=tmp.hkl end
end
""")
  l("stop")
  f = open("tmp.cns", "w")
  for l in cns_input:
    print >> f, l
  f.close()

def show_regression(x, y, label, min_correlation = 0):
  xy_regr = shared.linear_regression(x, y)
  assert xy_regr.is_well_defined()
  print label, "cc:", xy_regr.cc(), "m:", xy_regr.m()
  assert min_correlation == 0 or xy_regr.cc() >= min_correlation

def phase_error(p1, p2):
  d_as_r = math.pi / 180
  return math.acos(math.cos((p1 - p2) * d_as_r)) / d_as_r

class structure_factor_comparison:

  def __init__(self, label, min_corr_ampl=0, max_mean_w_phase_error=0,
               verbose=0):
    python_utils.adopt_init_args(self, locals())
    self.amp1 = shared.double()
    self.amp2 = shared.double()
    self.sum_amp1_minus_amp2_sq = 0
    self.sum_amp1_sq = 0
    self.sum_w_phase_error = 0
    self.sum_w = 0

  def add(self, h, f1, f2):
    a1, p1 = xutils.f_as_ampl_phase(f1)
    a2, p2 = xutils.f_as_ampl_phase(f2)
    if (self.verbose):
      print h
      print " ", a1, p1
      print " ", a2, p2
      print " " * 20, phase_error(p1, p2)
    self.amp1.append(a1)
    self.amp2.append(a2)
    self.sum_amp1_minus_amp2_sq += (a1 - a2)**2
    self.sum_amp1_sq += a1**2
    self.sum_w_phase_error += (a1 + a2) * phase_error(p1, p2)
    self.sum_w += (a1 + a2)

  def report(self):
    if (self.sum_amp1_sq):
      r = self.sum_amp1_minus_amp2_sq / self.sum_amp1_sq
      print "R-factor:", r
    if (self.sum_w):
      self.mean_w_phase_error = self.sum_w_phase_error / self.sum_w
    show_regression(
      self.amp1, self.amp2, self.label + " ampl", self.min_corr_ampl)
    print self.label + (" mean weighted phase error: %.2f" % (
      self.mean_w_phase_error,))
    if (self.max_mean_w_phase_error):
      assert self.mean_w_phase_error <= self.max_mean_w_phase_error

def show_structure_factor_correlation(label, h1, joined_sets, f1, f2,
                                      min_corr_ampl=0,
                                      max_mean_w_phase_error=0,
                                      verbose=0):
  sf_cmp = structure_factor_comparison(
    label, min_corr_ampl, max_mean_w_phase_error, verbose)
  if (joined_sets == 0):
    assert f1.size() == f2.size()
    for i in xrange(f1.size()):
      sf_cmp.add(h1[i], f1[i], f2[i])
  else:
    for i,j in joined_sets.pairs():
      sf_cmp.add(h1[i], f1[i], f2[j])
  sf_cmp.report()

def run_cns(elements, xtal, d_min, grid_resolution_factor,
            friedel_flag=0, fcalc=0):
  from cctbx.macro_mol import cns_input
  write_cns_input(elements, xtal, d_min, grid_resolution_factor, friedel_flag)
  try: os.unlink("tmp.hkl")
  except: pass
  os.system("cns < tmp.cns > tmp.out")
  f = open("tmp.hkl", "r")
  reader = cns_input.CNS_xray_reflection_Reader(f)
  reflection_file = reader.load()
  f.close()
  f_dir_h = reflection_file.reciprocal_space_objects["F_DIR"].H
  f_dir_f = reflection_file.reciprocal_space_objects["F_DIR"].data
  f_fft_h = reflection_file.reciprocal_space_objects["F_FFT"].H
  f_fft_f = reflection_file.reciprocal_space_objects["F_FFT"].data
  assert f_dir_h.size() == f_fft_h.size()
  show_structure_factor_correlation(
    "cns_dir/cns_fft", f_dir_h, 0, f_dir_f, f_fft_f, 0.99)
  if (fcalc):
    assert fcalc.H.size() == f_dir_h.size()
    asym_f_dir = sftbx.map_to_asym_index(
      xtal.SgInfo, friedel_flag, f_dir_h, f_dir_f)
    cns_h = asym_f_dir.asym_miller_indices()
    js = shared.join_sets(fcalc.H, cns_h)
    if (0):
      show_joined_sets(fcalc.H, cns_h, js)
    assert js.pairs().size() == fcalc.H.size()
    show_structure_factor_correlation(
      "sftbx_dir/cns_dir", fcalc.H, js, fcalc.F, asym_f_dir.asym_data_array(),
      0.99, 1.0, verbose=1)

def run():
  Flags = debug_utils.command_line_options(sys.argv[1:], (
    "RandomSeed",
    "AllSpaceGroups",
    "Isotropic",
    "Anisotropic",
    "Dispersive",
    "Anomalous",
    "ForceComplex",
    "cns",
    "shelx",
  ))
  if (not Flags.RandomSeed): debug_utils.set_random_seed(0)
  if (not (Flags.Isotropic or Flags.Anisotropic)):
    Flags.Isotropic = 1
    # XXX Flags.Anisotropic = 1
  symbols_to_stdout = 0
  if (len(sys.argv) > 1 + Flags.n):
    symbols = sys.argv[1:]
  else:
    symbols = debug_utils.get_test_space_group_symbols(Flags.AllSpaceGroups)
    symbols_to_stdout = 1
  for RawSgSymbol in symbols:
    if (RawSgSymbol.startswith("--")): continue
    SgSymbols = sgtbx.SpaceGroupSymbols(RawSgSymbol)
    SgInfo = sgtbx.SpaceGroup(SgSymbols).Info()
    LookupSymbol = SgInfo.BuildLookupSymbol()
    sys.stdout.flush()
    print >> sys.stderr, LookupSymbol
    sys.stderr.flush()
    if (symbols_to_stdout):
      print LookupSymbol
      sys.stdout.flush()
    exercise(
      SgInfo,
      adp=Flags.Anisotropic,
      force_complex=Flags.ForceComplex,
      random_f_prime=Flags.Dispersive,
      random_f_double_prime=Flags.Anomalous,
      use_cns=Flags.cns,
      use_shelx=Flags.shelx)
    sys.stdout.flush()

if (__name__ == "__main__"):
  run()
  t = os.times()
  print "u+s,u,s:", t[0] + t[1], t[0], t[1]
