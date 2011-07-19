from cctbx.array_family import flex
from cctbx import miller
from cctbx.development import make_cns_input
from cctbx.development import random_structure
from cctbx.regression.tst_miller import generate_random_hl
from iotbx.cns import reflection_reader
from libtbx.complex_math import abs_arg
from libtbx.test_utils import approx_equal
from libtbx import easy_run
import sys, os

def verify(sg_fcalc, sg_hl, sg_cns, p1_cns):
  sg_phase_integrals = miller.array(
    miller_set=miller.set(
      crystal_symmetry=sg_fcalc,
      indices=sg_cns.miller_indices,
      anomalous_flag=sg_fcalc.anomalous_flag()),
    data=sg_cns.hl).phase_integrals()
  for h, cns_pi, miller_pi in zip(sg_cns.miller_indices,
                                  sg_cns.pi.data,
                                  sg_phase_integrals.data()):
    if (abs(cns_pi - miller_pi) > 1.e-2):
      print "Error:", h, cns_pi, miller_pi
      if (0): return
      raise AssertionError
  #
  p1_phase_integrals = miller.array(
    miller_set=miller.set(
      crystal_symmetry=sg_fcalc.cell_equivalent_p1(),
      indices=p1_cns.miller_indices,
      anomalous_flag=sg_fcalc.anomalous_flag()),
    data=p1_cns.hl).phase_integrals()
  for h, cns_pi, miller_pi in zip(p1_cns.miller_indices,
                                  p1_cns.pi.data,
                                  p1_phase_integrals.data()):
    if (abs(cns_pi - miller_pi) > 1.e-2):
      print "Error:", h, cns_pi, miller_pi
      if (0): return
      raise AssertionError
  #
  space_group = sg_fcalc.space_group()
  asu = sg_fcalc.space_group_info().reciprocal_space_asu()
  lookup_dict = miller.make_lookup_dict(sg_fcalc.indices())
  for p1_i, h in enumerate(p1_cns.miller_indices):
    h_asu = miller.asym_index(space_group, asu, h)
    h_eq = h_asu.one_column(sg_fcalc.anomalous_flag())
    fcalc_asu = h_eq.complex_eq(p1_cns.fcalc.data[p1_i])
    hl_asu = h_eq.hendrickson_lattman_eq(p1_cns.hl[p1_i])
    sg_i = lookup_dict[h_eq.h()]
    assert abs(sg_fcalc.data()[sg_i] - fcalc_asu) < 1.e-2
    if (not approx_equal(sg_hl[sg_i], hl_asu, eps=1.e-2)):
      print "Error:", sg_fcalc.space_group_info()
      print sg_fcalc.indices()[sg_i]
      print "i:", sg_hl[sg_i]
      print "o:", hl_asu
      if (0): return
      raise AssertionError

def write_cns_input(fcalc_array, hl, test_merge=False):
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
  for i,h in enumerate(fcalc_array.indices()):
    l(  ("    index %d %d %d" % h)
      + (" fcalc=%.6g %.6g\n" % abs_arg(fcalc_array.data()[i], deg=True))
      + ("    pa=%.6g pb=%.6g pc=%.6g pd=%.6g" % hl[i]))
  l("  end")
  if (not test_merge):
    l("  declare name=pi domain=reciprocal type=complex end")
    l("  do (pi=get_fom[phistep=5,CEN360=false](pa,pb,pc,pd)) (all)")
    l("  write reflections output=\"tmp_sg.hkl\" end")
    l("  expand")
    l("  do (pi=get_fom[phistep=5,CEN360=false](pa,pb,pc,pd)) (all)")
    l("  write reflections output=\"tmp_p1.hkl\" end")
  else:
    l("  anomalous=false")
    l("  write reflections output=\"tmp_merged.hkl\" end")
  l("end")
  l("stop")
  f = open("tmp.cns", "w")
  for l in cns_input:
    print >> f, l
  f.close()

class read_reflection_arrays(object):

  def __init__(self, file_name, anomalous_flag, verbose):
    reflection_file = reflection_reader.cns_reflection_file(open(file_name))
    if (0 or verbose):
      print reflection_file.show_summary()
    assert reflection_file.anomalous == anomalous_flag
    names, self.miller_indices, self.hl = reflection_file.join_hl_group()
    self.fcalc = reflection_file.reciprocal_space_objects["FCALC"]
    self.pi = reflection_file.reciprocal_space_objects["PI"]
    assert not miller.match_indices(
      self.miller_indices, self.fcalc.indices).have_singles()
    assert not miller.match_indices(
      self.miller_indices, self.pi.indices).have_singles()

def exercise(space_group_info, anomalous_flag=False, d_min=2., verbose=0):
  sg_fcalc = random_structure.xray_structure(
    space_group_info,
    elements=("N", "C", "C", "O"),
    random_f_double_prime=anomalous_flag,
    random_u_iso=True,
    random_occupancy=True
    ).structure_factors(
      anomalous_flag=anomalous_flag, d_min=d_min, algorithm="direct").f_calc()
  sg_hl = generate_random_hl(sg_fcalc)
  write_cns_input(sg_fcalc, sg_hl.data())
  try: os.unlink("tmp_sg.hkl")
  except OSError: pass
  try: os.unlink("tmp_p1.hkl")
  except OSError: pass
  easy_run.fully_buffered(command="cns < tmp.cns > tmp.out") \
    .raise_if_errors_or_output()
  sg_cns = read_reflection_arrays("tmp_sg.hkl", anomalous_flag, verbose)
  p1_cns = read_reflection_arrays("tmp_p1.hkl", anomalous_flag, verbose)
  verify(sg_fcalc, sg_hl.data(), sg_cns, p1_cns)
  if (anomalous_flag):
    hl_merged = sg_hl.average_bijvoet_mates()
    fc_merged = sg_fcalc.average_bijvoet_mates()
    write_cns_input(sg_fcalc, sg_hl.data(), test_merge=True)
    try: os.unlink("tmp_merged.hkl")
    except OSError: pass
    easy_run.fully_buffered(command="cns < tmp.cns > tmp.out") \
      .raise_if_errors_or_output()
    reflection_file = reflection_reader.cns_reflection_file(
      open("tmp_merged.hkl"))
    if (not sg_fcalc.space_group().is_centric()):
      fc_merged_cns = reflection_file.reciprocal_space_objects["FCALC"]
      fc_merged_cns = fc_merged.customized_copy(
        indices=fc_merged_cns.indices,
        data=fc_merged_cns.data).map_to_asu().common_set(fc_merged)
      assert fc_merged_cns.indices().all_eq(fc_merged.indices())
      fc_merged_a = fc_merged.select_acentric()
      fc_merged_cns_a = fc_merged_cns.select_acentric()
      for part in [flex.real, flex.imag]:
        cc = flex.linear_correlation(
          part(fc_merged_a.data()),
          part(fc_merged_cns_a.data())).coefficient()
        if (cc < 1-1.e-6):
          print "FAILURE acentrics", sg_fcalc.space_group_info()
          if (0): return
          raise AssertionError
    names, miller_indices, hl = reflection_file.join_hl_group()
    assert names == ["PA", "PB", "PC", "PD"]
    hl_merged_cns = hl_merged.customized_copy(indices=miller_indices, data=hl)\
      .map_to_asu().common_set(hl_merged)
    assert hl_merged_cns.indices().all_eq(hl_merged.indices())
    for h,a,b in zip(hl_merged.indices(),
                     hl_merged.data(),
                     hl_merged_cns.data()):
      if (not approx_equal(a, b, eps=5.e-3)):
        print h
        print "cctbx:", a
        print "  cns:", b
        if (0): return
        raise AssertionError

def run_call_back(flags, space_group_info):
  for anomalous_flag in (False, True):
    exercise(space_group_info, anomalous_flag, verbose=flags.Verbose)

if (__name__ == "__main__"):
  make_cns_input.tst_run_requiring_cns(
    args=sys.argv[1:], call_back=run_call_back)
