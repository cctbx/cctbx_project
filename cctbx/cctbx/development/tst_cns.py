import sys, os, random
from cctbx_boost.arraytbx import shared
from cctbx_boost import sgtbx
from cctbx import xutils
from cctbx.development import debug_utils
from cctbx.development import make_cns_input
from cctbx.macro_mol import cns_input

def make_miller_lookup_dict(miller_indices): # XXX push to C++
  result = {}
  for i in xrange(miller_indices.size()):
    result[miller_indices[i]] = i
  return result

def generate_random_hl(miller_index_set, coeff_range=100):
  is_centric = miller_index_set.SgOps.isCentric
  hl = shared.hendrickson_lattman()
  for h in miller_index_set.H:
    if (is_centric(h)):
      hl.append((0,0,0,0))
    else:
      hl.append(
        [2 * coeff_range * random.random() - coeff_range for i in xrange(4)])
  return hl

def are_similar_hl(coeff1, coeff2, tolerance=1.e-2):
  for i in xrange(4):
    if (abs(coeff1[i] - coeff2[i]) > tolerance): return 0
  return 1

def verify(sg_miller_index_set, sg_fcalc_set, sg_hl,
           p1_miller_indices, p1_fcalc, p1_hl):
  sgops = sg_miller_index_set.SgOps
  asu = sgtbx.ReciprocalSpaceASU(sg_miller_index_set.SgInfo)
  lookup_dict = make_miller_lookup_dict(sg_miller_index_set.H)
  for p1_i in xrange(p1_miller_indices.size()):
    h_asu = sgtbx.Miller_AsymIndex(sgops, asu, p1_miller_indices[p1_i])
    h_eq = h_asu.one_column(sg_miller_index_set.friedel_flag)
    fcalc_asu = h_eq.complex_eq(p1_fcalc[p1_i])
    hl_asu = h_eq.hl_eq(p1_hl[p1_i])
    sg_i = lookup_dict[h_eq.H()]
    assert abs(sg_fcalc_set.F[sg_i] - fcalc_asu) < 1.e-2
    if (not are_similar_hl(sg_hl[sg_i], hl_asu)):
      print "Error:", sg_miller_index_set.SgInfo.BuildLookupSymbol()
      print sg_miller_index_set.H[sg_i]
      print "i:", sg_hl[sg_i]
      print "o:", hl_asu
      return # XXX
      raise AssertionError

def write_cns_input(miller_index_set, fcalc_set, hl):
  assert miller_index_set.H.size() == hl.size()
  cns_input = make_cns_input.unit_cell(miller_index_set.UnitCell)
  cns_input += make_cns_input.symmetry(miller_index_set.SgOps)
  l = cns_input.append
  l("xray")
  if (miller_index_set.friedel_flag):
    l("  anomalous=false")
  else:
    l("  anomalous=true")
  l("  declare name=fcalc domain=reciprocal type=complex end")
  l("  declare name=pa domain=reciprocal type=real end")
  l("  declare name=pb domain=reciprocal type=real end")
  l("  declare name=pc domain=reciprocal type=real end")
  l("  declare name=pd domain=reciprocal type=real end")
  l("  group type hl")
  l("    object pa")
  l("    object pb")
  l("    object pc")
  l("    object pd")
  l("  end")
  l("  reflection")
  for i in xrange(miller_index_set.H.size()):
    l(  ("    index %d %d %d" % miller_index_set.H[i])
      + (" fcalc=%.6g %.6g\n" % xutils.f_as_ampl_phase(fcalc_set.F[i], 1))
      + (" pa=%.6g pb=%.6g pc=%.6g pd=%.6g" % hl[i]))
  l("  end")
  l("  expand")
  l("  write reflections output=\"tmp.hkl\" end")
  l("end")
  l("stop")
  f = open("tmp.cns", "w")
  for l in cns_input:
    print >> f, l
  f.close()

def exercise(SgInfo, d_min=2., friedel_flag=0, verbose=0):
  elements = ("N", "C", "C", "O")
  xtal = debug_utils.random_structure(
    SgInfo, elements,
    volume_per_atom=50.,
    min_distance=1.5,
    general_positions_only=0)
  sg_miller_index_set = xutils.build_miller_indices(xtal, friedel_flag, d_min)
  sg_fcalc_set = xutils.calculate_structure_factors(sg_miller_index_set, xtal)
  sg_hl = generate_random_hl(sg_miller_index_set)
  write_cns_input(sg_miller_index_set, sg_fcalc_set, sg_hl)
  try: os.unlink("tmp.hkl")
  except: pass
  os.system("cns < tmp.cns > tmp.out")
  f = open("tmp.hkl", "r")
  reader = cns_input.CNS_xray_reflection_Reader(f)
  reflection_file = reader.load()
  f.close()
  if (0 or verbose):
    print reflection_file
  assert reflection_file.anomalous == (not sg_miller_index_set.friedel_flag)
  p1_miller_indices, p1_hl = reflection_file.join_hl_group()
  p1_fcalc_rso = reflection_file.reciprocal_space_objects["FCALC"]
  assert not shared.join_sets(p1_miller_indices, p1_fcalc_rso.H).have_singles()
  verify(sg_miller_index_set, sg_fcalc_set, sg_hl,
    p1_miller_indices, p1_fcalc_rso.data, p1_hl)

def run():
  Flags = debug_utils.command_line_options(sys.argv[1:], (
    "RandomSeed",
    "AllSpaceGroups",
    "Anomalous",
  ))
  if (not Flags.RandomSeed): debug_utils.set_random_seed(0)
  symbols_to_stdout = 0
  if (len(sys.argv) > 1 + Flags.n):
    symbols = sys.argv[1:]
  else:
    symbols = debug_utils.get_test_space_group_symbols(Flags.AllSpaceGroups)
    #symbols_to_stdout = 1
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
    exercise(SgInfo, friedel_flag=not Flags.Anomalous)
    sys.stdout.flush()

if (__name__ == "__main__"):
  run()
  t = os.times()
  print "u+s,u,s:", t[0] + t[1], t[0], t[1]
