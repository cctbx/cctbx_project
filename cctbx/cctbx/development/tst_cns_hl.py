import sys, os, random
from cctbx.array_family import flex
from cctbx import sgtbx
from cctbx import miller
from cctbx.development import debug_utils
from cctbx.development import make_cns_input
from cctbx.development import random_structure
from cctbx.macro_mol import cns_input
from scitbx.python_utils.complex_math import abs_arg

def generate_random_hl(miller_set, coeff_range=100):
  is_centric = miller_set.space_group().is_centric
  hl = flex.hendrickson_lattman()
  for h in miller_set.indices():
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

def verify(sg_fcalc_array, sg_hl,
           p1_miller_indices, p1_fcalc, p1_hl):
  space_group = sg_fcalc_array.space_group()
  asu = sg_fcalc_array.space_group_info().reciprocal_space_asu()
  lookup_dict = miller.make_lookup_dict(sg_fcalc_array.indices())
  for p1_i, p1_h in p1_miller_indices.items():
    h_asu = miller.asym_index(space_group, asu, p1_h)
    h_eq = h_asu.one_column(sg_fcalc_array.anomalous_flag())
    fcalc_asu = h_eq.complex_eq(p1_fcalc[p1_i])
    hl_asu = h_eq.hendrickson_lattman_eq(p1_hl[p1_i])
    sg_i = lookup_dict[h_eq.h()]
    assert abs(sg_fcalc_array.data()[sg_i] - fcalc_asu) < 1.e-2
    if (not are_similar_hl(sg_hl[sg_i], hl_asu)):
      print "Error:", sg_fcalc_array.space_group_info()
      print sg_fcalc_array.indices()[sg_i]
      print "i:", sg_hl[sg_i]
      print "o:", hl_asu
      if (0): return
      raise AssertionError

def write_cns_input(fcalc_array, hl):
  assert fcalc_array.data().size() == hl.size()
  cns_input = make_cns_input.xray_unit_cell(fcalc_array.unit_cell())
  cns_input += make_cns_input.xray_symmetry(fcalc_array.space_group())
  cns_input += make_cns_input.xray_anomalous(fcalc_array.anomalous_flag())
  l = cns_input.append
  l("xray")
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
  for i,h in fcalc_array.indices().items():
    l(  ("    index %d %d %d" % h)
      + (" fcalc=%.6g %.6g\n" % abs_arg(fcalc_array.data()[i], deg=True))
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

def exercise(space_group_info, anomalous_flag=False, d_min=2., verbose=0):
  sg_fcalc_array = random_structure.xray_structure(
    space_group_info,
    elements=("N", "C", "C", "O"),
    random_u_iso=True,
    random_occupancy=True
    ).structure_factors_direct(
      anomalous_flag=anomalous_flag, d_min=d_min).f_calc_array()
  sg_hl = generate_random_hl(sg_fcalc_array)
  write_cns_input(sg_fcalc_array, sg_hl)
  try: os.unlink("tmp.hkl")
  except: pass
  os.system("cns < tmp.cns > tmp.out")
  f = open("tmp.hkl", "r")
  reflection_file = cns_input.cns_reflection_file(f)
  f.close()
  if (0 or verbose):
    print reflection_file.show_summary()
  assert reflection_file.anomalous == sg_fcalc_array.anomalous_flag()
  p1_miller_indices, p1_hl = reflection_file.join_hl_group()
  p1_fcalc_rso = reflection_file.reciprocal_space_objects["FCALC"]
  assert not miller.match_indices(
    p1_miller_indices, p1_fcalc_rso.indices).have_singles()
  verify(sg_fcalc_array, sg_hl, p1_miller_indices, p1_fcalc_rso.data, p1_hl)

def run_call_back(flags, space_group_info):
  for anomalous_flag in (False, True):
    exercise(space_group_info, anomalous_flag, verbose=flags.Verbose)

def run():
  make_cns_input.check_cns_availability()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
