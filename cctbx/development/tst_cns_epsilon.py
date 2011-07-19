from cctbx import miller
from cctbx import crystal
from cctbx.development import make_cns_input
from iotbx.cns import reflection_reader
from libtbx import easy_run
import sys, os

def verify(crystal_symmetry, anomalous_flag, reflection_file):
  assert reflection_file.anomalous == anomalous_flag
  cns_m = reflection_file.reciprocal_space_objects["CNS_M"]
  cns_e = reflection_file.reciprocal_space_objects["CNS_E"]
  cns_c = reflection_file.reciprocal_space_objects["CNS_C"]
  cns_a = reflection_file.reciprocal_space_objects["CNS_A"]
  cns_p = reflection_file.reciprocal_space_objects["CNS_P"]
  for cns_x in (cns_e, cns_c, cns_a, cns_p):
    assert cns_x.indices.id() == cns_m.indices.id()
  space_group = crystal_symmetry.space_group()
  for i,h in enumerate(cns_m.indices):
    m_i = cns_m.data[i]
    e_i = cns_e.data[i]
    c_i = cns_c.data[i]
    a_i = cns_a.data[i]
    p_i = cns_p.data[i]

    if (p_i < 0 and p_i != -1.): p_i += 180.

    assert (c_i == 0) == (a_i != 0)
    assert (c_i == 0 and p_i == -1.) or (c_i != 0 and p_i != -1.), \
           'c_i = %d, p_i = %g' % (c_i, p_i)

    m = space_group.multiplicity(h, anomalous_flag)
    e = space_group.epsilon(h)
    c = space_group.is_centric(h)
    p = space_group.phase_restriction(h).ht_angle(True)

    if (c or anomalous_flag):
      assert e == space_group.order_p() // m
    else:
      assert e == (2 * space_group.order_p()) // m

    try:
      assert m_i == m, 'multiplicity mismatch'
      assert e_i == e, 'epsilon mismatch'
      assert c_i == c, 'centric mismatch'
      assert p_i == p, 'restricted phase mismatch'
    except AssertionError, exc:
      print crystal_symmetry.space_group_info()
      print 'index=', h
      print 'm:', m_i, m
      print 'e:', e_i, e
      print 'c:', c_i, c
      print 'p:', p_i, p
      raise AssertionError, exc

    assert (not space_group.is_sys_absent(h))
    assert (e == space_group.order_p() // space_group.multiplicity(h, True))
    eq = miller.sym_equiv_indices(space_group, h)
    assert (eq.multiplicity(anomalous_flag) == m)
    assert (eq.epsilon() == e)
    assert (eq.is_centric() == c)
    assert (eq.phase_restriction().ht_angle(True) == p)

def write_cns_input(crystal_symmetry, anomalous_flag, d_min):
  cns_input = make_cns_input.xray_unit_cell(crystal_symmetry.unit_cell())
  cns_input += make_cns_input.xray_symmetry(crystal_symmetry.space_group())
  cns_input += make_cns_input.xray_anomalous(anomalous_flag)
  cns_input += make_cns_input.xray_generate(10000, d_min)
  l = cns_input.append
  l("xray")
  l("  declare name=cns_m domain=reciprocal type=integer end")
  l("  declare name=cns_e domain=reciprocal type=integer end")
  l("  declare name=cns_c domain=reciprocal type=integer end")
  l("  declare name=cns_a domain=reciprocal type=integer end")
  l("  declare name=cns_p domain=reciprocal type=real    end")
  l("  do (cns_m = mult) (all)")
  l("  do (cns_e = epsilon) (all)")
  l("  do (cns_c = 0) (all)")
  l("  do (cns_c = 1) (centric)")
  l("  do (cns_a = 0) (all)")
  l("  do (cns_a = 1) (acentric)")
  l("  do (cns_p = -1) (all)")
  l("  do (cns_p = centric_phase) (centric)")
  l("  write reflections output=\"tmp.hkl\" end")
  l("end")
  l("stop")
  f = open("tmp.cns", "w")
  for l in cns_input:
    print >> f, l
  f.close()

def exercise(space_group_info, anomalous_flag=False, d_min=2., verbose=0):
  crystal_symmetry = crystal.symmetry(
    space_group_info.any_compatible_unit_cell(1000),
    space_group_info=space_group_info)
  write_cns_input(crystal_symmetry, anomalous_flag, d_min)
  try: os.unlink("tmp.hkl")
  except KeyboardInterrupt: raise
  except Exception: pass
  easy_run.fully_buffered(command="cns < tmp.cns > tmp.out") \
    .raise_if_errors_or_output()
  f = open("tmp.hkl", "r")
  reflection_file = reflection_reader.cns_reflection_file(f)
  f.close()
  if (0 or verbose):
    print reflection_file.show_summary()
  verify(crystal_symmetry, anomalous_flag, reflection_file)

def run_call_back(flags, space_group_info):
  for anomalous_flag in (False, True):
    exercise(space_group_info, anomalous_flag, verbose=flags.Verbose)

if (__name__ == "__main__"):
  make_cns_input.tst_run_requiring_cns(
    args=sys.argv[1:], call_back=run_call_back)
