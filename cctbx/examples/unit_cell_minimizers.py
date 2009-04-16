from cctbx import uctbx
from cctbx.eltbx import wavelengths
from cctbx.array_family import flex
import scitbx.lbfgs
import scitbx.minimizers
import libtbx.utils

def residual(
      two_thetas_obs, miller_indices, wavelength, unit_cell):
  two_thetas_calc = unit_cell.two_theta(miller_indices, wavelength, deg=True)
  return flex.sum(flex.pow2(two_thetas_obs - two_thetas_calc))

def gradients(
      two_thetas_obs, miller_indices, wavelength, unit_cell, eps=1.e-6):
  result = flex.double()
  for i in xrange(6):
    rs = []
    for signed_eps in [eps, -eps]:
      params_eps = list(unit_cell.parameters())
      params_eps[i] += signed_eps
      rs.append(
        residual(
          two_thetas_obs, miller_indices, wavelength,
          uctbx.unit_cell(params_eps)))
    result.append((rs[0]-rs[1])/(2*eps))
  return result

def hessian(
      two_thetas_obs, miller_indices, wavelength, unit_cell, eps=1.e-6):
  result = flex.double()
  for i in xrange(6):
    rs = []
    for signed_eps in [eps, -eps]:
      params_eps = list(unit_cell.parameters())
      params_eps[i] += signed_eps
      rs.append(
        gradients(
          two_thetas_obs, miller_indices, wavelength,
          uctbx.unit_cell(params_eps)))
    result.extend((rs[0]-rs[1])/(2*eps))
  result.reshape(flex.grid(6,6))
  u = result.matrix_symmetric_as_packed_u(relative_epsilon=eps*10)
  return u.matrix_packed_u_as_symmetric()

class refinery:

  def __init__(self,
        two_thetas_obs, miller_indices, wavelength, unit_cell,
        mode):
    self.two_thetas_obs = two_thetas_obs
    self.miller_indices = miller_indices
    self.wavelength = wavelength
    #
    self.number_of_function_evaluations = 0
    self.progress = []
    if (mode < 3):
      self.x = flex.double(unit_cell.parameters())
      self.diag_mode = [None, "once", "always"][mode]
      scitbx.lbfgs.run(target_evaluator=self)
      self.plot_legend = "lbfgs_"+str(self.diag_mode)
    else:
      assert mode == 3
      m = scitbx.minimizers.newton_more_thuente_1994(
        function=self, x0=flex.double(unit_cell.parameters()))
      m.show_statistics()
      self.x = m.x_star
      self.plot_legend = "newton"

  def unit_cell(self):
    return uctbx.unit_cell(iter(self.x))

  def compute_functional_gradients_diag(self, compute_diag=True):
    unit_cell = self.unit_cell()
    f = residual(
      self.two_thetas_obs, self.miller_indices, self.wavelength, unit_cell)
    if (len(self.progress) == 0):
      self.progress.append((0, f))
    self.number_of_function_evaluations += 1
    g = gradients(
      self.two_thetas_obs, self.miller_indices, self.wavelength, unit_cell)
    print "functional: %12.6g" % f, "gradient norm: %12.6g" % g.norm()
    if (not compute_diag):
      return f, g
    h = hessian(
      self.two_thetas_obs, self.miller_indices, self.wavelength, unit_cell)
    d = h.matrix_diagonal()
    return f, g, 1/d

  def compute_functional_and_gradients(self):
    return self.compute_functional_gradients_diag(compute_diag=False)

  def callback_after_step(self, minimizer=None, x=None):
    if (minimizer is not None):
      functional = self.functional(x=self.x)
      print "STEP lbfgs", functional
    else:
      functional = self.functional(x=x)
      print "STEP newton_more_thuente_1994", functional
    self.progress.append((self.number_of_function_evaluations, functional))

  def functional(self, x):
    result = residual(
      self.two_thetas_obs, self.miller_indices, self.wavelength,
      unit_cell=uctbx.unit_cell(iter(x)))
    if (len(self.progress) == 0):
      self.progress.append((0, result))
    return result

  def gradients(self, x):
    self.number_of_function_evaluations += 1
    return gradients(
      self.two_thetas_obs, self.miller_indices, self.wavelength,
      unit_cell=uctbx.unit_cell(iter(x)))

  def hessian(self, x):
    return hessian(
      self.two_thetas_obs, self.miller_indices, self.wavelength,
      unit_cell=uctbx.unit_cell(iter(x)))

def show_fit(two_thetas_obs, miller_indices, wavelength, unit_cell):
  two_thetas_calc = unit_cell.two_theta(miller_indices, wavelength, deg=True)
  for h,o,c in zip(miller_indices, two_thetas_obs, two_thetas_calc):
    print "(%2d, %2d, %2d)" % h, "%6.2f - %6.2f = %6.2f" % (o, c, o-c)
  print

two_theta_and_index_list = """\
  8.81   0  1  1
 12.23   0  0  2
 12.71   0  2  0
 12.97   1  1  0
 13.79   0  1  2
 14.11   0  2  1
 14.35   1  1  1
 16.68   1  0  2
 17.03   1  2  0
 17.67   0  2  2
 17.86   1  1  2
 19.47   0  1  3
 21.03   1  2  2
 22.26   1  3  0
 22.41   0  2  3
 22.56   1  1  3
 22.72   2  0  0
 23.10   1  3  1
 24.40   2  1  1
 24.60   0  0  4
 25.17   1  2  3
 25.43   0  1  4
 25.87   2  0  2
 26.11   2  2  0
 26.32   0  4  1
 26.66   2  1  2
 26.84   2  2  1
 27.15   1  0  4
 27.78   0  2  4
 27.90   1  1  4
 28.44   0  4  2
 28.72   1  4  1
 28.92   2  2  2
 29.02   1  3  3
 30.08   2  1  3
 30.49   2  3  1
 30.69   1  4  2
 31.34   0  3  4
 31.56   0  1  5
 32.12   2  2  3
""".splitlines()

def run():
  two_thetas_obs = flex.double()
  miller_indices = flex.miller_index()
  for line in two_theta_and_index_list:
    fields = line.split()
    assert len(fields) == 4
    two_thetas_obs.append(float(fields[0]))
    miller_indices.append([int(s) for s in fields[1:]])

  wavelength = wavelengths.characteristic("CU").as_angstrom()
  unit_cell_start = uctbx.unit_cell((10,10,10,90,90,90))
  show_fit(
    two_thetas_obs, miller_indices, wavelength, unit_cell_start)

  refined_accu = []
  for mode in range(4):
    refined = refinery(
      two_thetas_obs, miller_indices, wavelength, unit_cell_start, mode=mode)
    refined_accu.append(refined)
    print

  p = open("tmp.xy", "w")
  print >> p, "@with g0"
  for i,r in enumerate(refined_accu):
    print >> p, '@ s%d legend "%s"' % (i, r.plot_legend)
    print >> p, '@ s%d symbol 1' % i
  for refined in refined_accu:
    for x,y in refined.progress:
      if (x > 15): break
      print >> p, x,y
    print >> p, "&"
  del p

  if (0):
    show_fit(
      two_thetas_obs, miller_indices, wavelength, refined.unit_cell())

  print refined.unit_cell()
  print
  print libtbx.utils.format_cpu_times()

if (__name__ == "__main__"):
  run()
