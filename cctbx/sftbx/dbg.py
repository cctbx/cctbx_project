import sys, os
from cctbx_boost.arraytbx import shared
from cctbx_boost import sgtbx
from cctbx_boost import adptbx
from cctbx_boost import sftbx
from cctbx_boost import fftbx
from cctbx import xutils
from cctbx.development import debug_utils
from cctbx.development import make_cns_input

def add_u_extra(xtal, u_extra):
  xtal_mod = xtal.copy_attributes()
  for site in xtal:
    site.set_Uiso(site.Uiso() + u_extra)
    xtal_mod.add_site(site)
  return xtal_mod

def print_structure_factors(SgInfo,
                            adp=0,
                            d_min=3.,
                            grid_resolution_factor = 1./3,
                            use_cns=0):
  elements = ("N", "C", "C", "O", "N", "C", "C", "O")
  xtal = debug_utils.random_structure(
    SgInfo, elements,
    volume_per_atom=50.,
    min_distance=1.5,
    general_positions_only=0,
    anisotropic_displacement_parameters=adp)
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
  MillerIndices = xutils.build_miller_indices(xtal, d_min)
  Fcalc = xutils.calculate_structure_factors(MillerIndices, xtal)
  if (use_cns):
    run_cns(elements, xtal, d_min, grid_resolution_factor)
    return
  max_q = 1. / (d_min**2)
  max_prime = 5
  mandatory_grid_factors = xtal.SgOps.refine_gridding()
  grid_logical = sftbx.determine_grid(
    xtal.UnitCell,
    max_q, grid_resolution_factor, max_prime, mandatory_grid_factors)
  fft = fftbx.real_to_complex_3d(grid_logical)
  quality_factor = 100
  u_extra = sftbx.calc_u_extra(max_q, grid_resolution_factor, quality_factor)
  if (0):
    u_extra = adptbx.B_as_U(20)
  wing_cutoff = 1.e-3
  exp_table_one_over_step_size = -100
  sampled_density = sftbx.sampled_model_density(
    xtal.UnitCell, xtal.Sites,
    fft.Nreal(), fft.Mreal(),
    u_extra, wing_cutoff, exp_table_one_over_step_size)
  print "u_extra:", sampled_density.u_extra(),
  print "b_extra:", adptbx.U_as_B(sampled_density.u_extra())
  if (1):
    print "wing_cutoff:", sampled_density.wing_cutoff()
    print "exp_table_one_over_step_size:", \
      sampled_density.exp_table_one_over_step_size()
    print "exp_table_size:", sampled_density.exp_table_size()
  tags = sftbx.grid_tags(fft.Nreal())
  sym_flags = sftbx.map_symmetry_flags(1)
  tags.build(xtal.SgInfo, sym_flags)
  sampled_density.apply_symmetry(tags)
  map = sampled_density.map_as_shared()
  map_stats = shared.statistics(map)
  if (0):
    print "Electron density"
    print "max %.6g" % (map_stats.max())
    print "min %.6g" % (map_stats.min())
    print "mean %.6g" % (map_stats.mean())
    print "sigma %.6g" % (map_stats.sigma())
  if (0):
    map = sftbx.structure_factor_map(
      xtal.SgOps, Fcalc.H, Fcalc.F, fft.Ncomplex())
    fft.backward(map)
    map_stats = shared.statistics(map)
    print "True electron density"
    print "max %.6g" % (map_stats.max())
    print "min %.6g" % (map_stats.min())
    print "mean %.6g" % (map_stats.mean())
    print "sigma %.6g" % (map_stats.sigma())
  fft.forward(map)
  map_stats = shared.statistics(map)
  if (0):
    print "Transformed electron density"
    print "max %.6g" % (map_stats.max())
    print "min %.6g" % (map_stats.min())
    print "mean %.6g" % (map_stats.mean())
    print "sigma %.6g" % (map_stats.sigma())
    print "Ncomplex", fft.Ncomplex()
  if (0):
    map = sftbx.structure_factor_map(
      xtal.SgOps, Fcalc.H, Fcalc.F, fft.Ncomplex())
  miller_indices, fcal = sftbx.collect_structure_factors(
    xtal.UnitCell, xtal.SgInfo,
    max_q, map, fft.Ncomplex())
  sampled_density.eliminate_u_extra_and_normalize(miller_indices, fcal)
  if (0):
    u_extra = sampled_density.u_extra()
    xtal_extra = add_u_extra(xtal, u_extra)
    fcalc_extra = xutils.calculate_structure_factors(MillerIndices, xtal_extra)
    show_structure_factor_correlation("before", Fcalc.F, fcalc_extra.F)
    sftbx.eliminate_u_extra(
      xtal.UnitCell, u_extra, MillerIndices.H, fcalc_extra.F)
    show_structure_factor_correlation("after", Fcalc.F, fcalc_extra.F)
  js = shared.join_sets(MillerIndices.H, miller_indices)
  if (0):
    for i,j in js.pairs():
      print MillerIndices.H[i], miller_indices[j]
    print "singles 1:"
    for i in js.singles(0):
      print MillerIndices.H[i]
    print "singles 2:"
    for i in js.singles(1):
      print miller_indices[i]
  assert js.pairs().size() + js.singles(0).size() == MillerIndices.H.size()
  assert js.pairs().size() + js.singles(1).size() == miller_indices.size()
  assert js.pairs().size() == MillerIndices.H.size()
  for i in xrange(2):
    assert js.singles(i).size() == 0
  x = shared.double()
  y = shared.double()
  for i,j in js.pairs():
    assert MillerIndices.H[i] == miller_indices[j]
    if (0):
      print MillerIndices.H[i],
      print "(" + debug_utils.format_structure_factor(Fcalc.F[i]) + ")",
      print "(" + debug_utils.format_structure_factor(fcal[j]) + ")"
    x.append(abs(Fcalc.F[i]))
    y.append(abs(fcal[j]))
  xy_regr = shared.linear_regression(x, y)
  assert xy_regr.is_well_defined()
  print "cc:", xy_regr.cc(), "m:", xy_regr.m()

def write_cns_input(elements, xtal, d_min, grid_resolution_factor):
  cns_input = make_cns_input.topology(elements)
  cns_input += make_cns_input.coordinates(xtal.Sites)
  cns_input += make_cns_input.unit_cell(xtal.UnitCell)
  cns_input += make_cns_input.symmetry(xtal.SgOps)
  l = cns_input.append
  l("""write coordinates end
coordinates orthogonalize end

xray
  @@CNS_XRAYLIB:scatter.lib
  fft
    grid=%.12g
  end
end"""
% (grid_resolution_factor,))
  l("""
xray
  ANOMalous=FALSe
  generate 100000. %.12g
  mapresolution %.12g
end"""
% (d_min, d_min))
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

def show_structure_factor_correlation(label, f1, f2,
                                      min_corr_ampl = 0, min_corr_phases = 0):
  assert f1.size() == f2.size()
  a1 = shared.double()
  p1 = shared.double()
  a2 = shared.double()
  p2 = shared.double()
  for i in xrange(f1.size()):
    a, p = xutils.f_as_ampl_phase(f1[i])
    a1.append(a)
    p1.append(p)
    a, p = xutils.f_as_ampl_phase(f2[i])
    a2.append(a)
    p2.append(p)
  show_regression(a1, a2, label + " ampl", min_corr_ampl)
  show_regression(p1, p2, label + " phases", min_corr_phases)

def run_cns(elements, xtal, d_min, grid_resolution_factor, fcalc = 0):
  from cctbx.macro_mol import cns_input
  write_cns_input(elements, xtal, d_min, grid_resolution_factor)
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
  ampl_dir = shared.double()
  phase_dir = shared.double()
  ampl_fft = shared.double()
  phase_fft = shared.double()
  for i in xrange(f_dir_h.size()):
    assert f_dir_h[i] == f_fft_h[i]
    a, p = xutils.f_as_ampl_phase(f_dir_f[i])
    ampl_dir.append(a)
    phase_dir.append(p)
    a, p = xutils.f_as_ampl_phase(f_fft_f[i])
    ampl_fft.append(a)
    phase_fft.append(p)
  show_regression(ampl_dir, ampl_fft, "ampl dir/fft", 0.99)
  show_regression(phase_dir, phase_fft, "phase dir/fft")
  if (fcalc):
    # XXX this does not work:
    # XXX need to map to common asymmetric unit
    assert 0
    js = shared.join_sets(f_dir_h, fcalc.H)
    assert js.pairs().size() == f_dir_h.size()
    ampl_fcalc = shared.double()
    phase_fcalc = shared.double()
    for i,j in js.pairs():
      a, p = xutils.f_as_ampl_phase(fcalc.F[j])
      ampl_fcalc.append(a)
      phase_fcalc.append(p)
    show_regression(ampl_dir, ampl_fcalc, "ampl dir/fcalc", 0.99)
    show_regression(phase_dir, phase_fcalc, "phase dir/fcalc", 0.99)

def run():
  Flags = debug_utils.command_line_options(sys.argv[1:], (
    "RandomSeed",
    "AllSpaceGroups",
    "Isotropic",
    "Anisotropic",
    "cns",
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
    if (Flags.Isotropic):
      print_structure_factors(SgInfo, adp=0, use_cns=Flags.cns)
    if (Flags.Anisotropic):
      print_structure_factors(SgInfo, adp=1)
    sys.stdout.flush()

if (__name__ == "__main__"):
  run()
  t = os.times()
  print "u+s,u,s:", t[0] + t[1], t[0], t[1]
