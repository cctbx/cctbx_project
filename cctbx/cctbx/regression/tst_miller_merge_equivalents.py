from cctbx import miller
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
import random
import sys

def exercise(space_group_info, anomalous_flags,
             n_scatterers=8, d_min=2, verbose=0):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["const"]*n_scatterers)
  f = abs(structure.structure_factors(
    d_min=d_min, anomalous_flag=anomalous_flags).f_calc())
  fs = miller.array(miller_set=f, data=f.data(), sigmas=flex.sqrt(f.data()))
  for a in (f, fs):
    m = a.merge_equivalents()
    j = m.array().adopt_set(a)
    assert flex.linear_correlation(j.data(),
                                   a.data()).coefficient() > 1-1.e-6
    if (a.sigmas() is not None):
      assert flex.linear_correlation(j.sigmas(),
                                     a.sigmas()).coefficient() > 1-1.e-6
  redundancies = flex.size_t()
  for i in fs.indices().indices():
    redundancies.append(random.randrange(5)+1)
  space_group = space_group_info.group()
  r_indices = flex.miller_index()
  r_data = flex.double()
  r_sigmas = flex.double()
  for i,n in redundancies.items():
    h = fs.indices()[i]
    h_eq = miller.sym_equiv_indices(space_group, h).indices()
    for j in xrange(n):
      r_indices.append(h_eq[random.randrange(len(h_eq))].h())
      r_data.append(fs.data()[i])
      r_sigmas.append(fs.sigmas()[i])
  r = miller.array(
    miller_set=miller.set(
      crystal_symmetry=fs,
      indices=r_indices,
      anomalous_flag=fs.anomalous_flag()),
    data=r_data,
    sigmas=r_sigmas)
  noise = flex.double()
  for i in r.indices().indices():
    noise.append(random.random())
  r = r.sort(by_value=noise)
  m = r.merge_equivalents()
  j = m.array().adopt_set(fs)
  assert flex.linear_correlation(
    j.data(),
    fs.data()).coefficient() > 1-1.e-6
  assert flex.linear_correlation(
    j.sigmas()*flex.sqrt(redundancies.as_double()),
    fs.sigmas()).coefficient() > 1-1.e-6

def run_call_back(flags, space_group_info):
  for anomalous_flag in (00000, 0001):
    exercise(space_group_info, anomalous_flag)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print "OK"

if (__name__ == "__main__"):
  run()
