from cctbx import mintbx
from cctbx import miller
from cctbx import adptbx
from scitbx import lbfgs
from cctbx.array_family import flex
from cctbx.development import debug_utils
from cctbx.development import random_structure
from scitbx.python_utils.misc import adopt_init_args
import sys

class k_b_scaling_minimizer:

  def __init__(self, unit_cell, miller_indices, multiplicities,
               data_reference, data_scaled,
               k_initial, b_initial,
               refine_k, refine_b,
               min_iterations=50, max_calls=1000):
    adopt_init_args(self, locals())
    self.anisotropic = hasattr(self.b_initial, "__len__")
    self.k_min = 1 # refine correction factor for k_initial
    self.b_min = self.b_initial
    self.x = self.pack(self.k_min, self.b_min)
    self.n = self.x.size()
    self.minimizer = lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs.termination_parameters(
        min_iterations=min_iterations,
        max_calls=max_calls))
    self()
    del self.x
    self.k_min *= self.k_initial

  def pack(self, k, b):
    v = []
    if (self.refine_k): v.append(k)
    if (self.refine_b):
      if (self.anisotropic): v += list(b)
      else:                  v.append(b)
    return flex.double(tuple(v))

  def unpack_x(self):
    i = 0
    if (self.refine_k):
      self.k_min = self.x[i]
      i += 1
    if (self.refine_b):
      if (self.anisotropic):
        self.b_min = tuple(self.x)[i:]
      else:
        self.b_min = self.x[i]

  def __call__(self):
    self.unpack_x()
    tg = mintbx.k_b_scaling_target_and_gradients(
      self.unit_cell, self.miller_indices, self.multiplicities,
      self.data_reference, self.data_scaled,
      self.k_initial * self.k_min, self.b_min,
      self.refine_k, self.refine_b)
    self.f = tg.target()
    if (self.anisotropic):
      self.g = self.pack(tg.gradient_k(), tg.gradients_b_cif())
    else:
      self.g = self.pack(tg.gradient_k(), tg.gradient_b_iso())
    return self.x, self.f, self.g

def exercise(space_group_info, anomalous_flag, d_min=2., verbose=0):
  structure_factors = random_structure.xray_structure(
    space_group_info,
    elements=("N", "C", "C", "O", "N", "C", "C", "O"),
    volume_per_atom=50.,
    min_distance=1.5,
    anisotropic_flag=True,
    random_f_prime_d_min=1.0,
    random_f_double_prime=anomalous_flag,
    random_u_iso=True,
    random_occupancy=True
    ).structure_factors(
        anomalous_flag=anomalous_flag, d_min=d_min, algorithm="direct")
  if (0 or verbose):
    structure_factors.xray_structure().show_summary()
  f_ref = abs(structure_factors.f_calc())
  multiplicities = f_ref.multiplicities()
  for anisotropic_flag in (False, True):
    f_sca = miller.array(miller_set=f_ref, data=flex.double())
    k_sim = 100
    b_iso = 13
    b_cif = [5,10,15,20,25,30]
    if (anisotropic_flag):
      u_star = adptbx.u_cif_as_u_star(
        f_ref.unit_cell(), adptbx.b_as_u(b_cif))
    for i,h in enumerate(f_ref.indices()):
      if (anisotropic_flag):
        dw = adptbx.debye_waller_factor_u_star(h, u_star)
      else:
        dw = adptbx.debye_waller_factor_b_iso(f_ref.unit_cell(),h,b_iso)
      f_sca.data().push_back(k_sim * f_ref.data()[i] * dw)
      if (0 or verbose):
        print h, f_ref.data()[i], f_sca.data()[i]
    k_min = 1
    if (anisotropic_flag):
      b_min = [0,0,0,0,0,0]
    else:
      b_min = 0
    for p in xrange(10):
      for refine_k, refine_b in ((1,0), (0,1), (1,0), (0,1), (1,0), (1,1)):
        minimized = k_b_scaling_minimizer(
          f_ref.unit_cell(),
          f_ref.indices(),
          multiplicities.data(),
          f_ref.data(),
          f_sca.data(),
          k_min, b_min,
          refine_k, refine_b)
        k_min = minimized.k_min
        b_min = minimized.b_min
    if (0 or verbose):
      print "anisotropic_flag:", anisotropic_flag
      print "k_min:", minimized.k_min
      print "b_min:", minimized.b_min
      print "target:", minimized.f,
      print "after %d iteration(s)" % (minimized.minimizer.iter(),)
    assert abs(minimized.k_min - k_sim) < 1.e-4
    if (anisotropic_flag):
      for i in xrange(6):
        assert abs(minimized.b_min[i] - b_cif[i]) < 1.e-4
    else:
      assert abs(minimized.b_min - b_iso) < 1.e-4
  if (0 or verbose):
    print

def run_call_back(flags, space_group_info):
  for anomalous_flag in (False, True):
    exercise(space_group_info, anomalous_flag, verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
