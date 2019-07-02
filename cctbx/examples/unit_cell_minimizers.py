from __future__ import absolute_import, division, print_function
from cctbx import uctbx
from cctbx.eltbx import wavelengths
from cctbx.array_family import flex
import scitbx.lbfgs
import scitbx.lbfgsb
import scitbx.minimizers
import libtbx.utils
import platform
import sys
from six.moves import range
from six.moves import zip

def residual(
      two_thetas_obs, miller_indices, wavelength, unit_cell):
  two_thetas_calc = unit_cell.two_theta(miller_indices, wavelength, deg=True)
  return flex.sum(flex.pow2(two_thetas_obs - two_thetas_calc))

def gradients(
      two_thetas_obs, miller_indices, wavelength, unit_cell, eps=1.e-6):
  result = flex.double()
  for i in range(6):
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
  for i in range(6):
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
    self.mode = mode
    self.number_of_gradient_evaluations = 0
    self.functionals = []
    if (mode < 2):
      self.plot_legend = "%d:newton_%s" % (mode, ["full", "diag"][mode])
      m = scitbx.minimizers.newton_more_thuente_1994(
        function=self, x0=flex.double(unit_cell.parameters()))
      m.show_statistics()
      self.x = m.x_star
    elif (mode < 7):
      diagco, use_hessian, lbfgs_impl_switch = [
        (0, 0, 0),
        (0, 0, 1),
        (1, 1, 0),
        (1, 1, 1),
        (2, 1, 1)][mode-2]
      self.plot_legend = "%d:lbfgs_d=%d_u=%d_l=%d" % (
        mode, diagco, use_hessian, lbfgs_impl_switch)
      print("plot_legend:", self.plot_legend)
      self.x = self.run_lbfgs_raw(
        unit_cell=unit_cell,
        diagco=diagco,
        use_hessian=use_hessian,
        lbfgs_impl_switch=lbfgs_impl_switch)
    elif (mode < 8):
      self.plot_legend = "%d:lbfgsb" % mode
      print("plot_legend:", self.plot_legend)
      self.x = self.run_lbfgsb(unit_cell=unit_cell)
    else:
      raise AssertionError("bad mode=%d" % mode)

  def unit_cell(self):
    return uctbx.unit_cell(iter(self.x))

  def functional(self, x):
    if (0):
      print("functional(): x =", list(x))
    if (flex.min(x[:3]) < 1):
      print("FunctionalException: small length")
      raise scitbx.minimizers.FunctionalException
    if (flex.min(x[3:]) < 50):
      print("FunctionalException: small angle")
      raise scitbx.minimizers.FunctionalException
    try:
      result = residual(
        self.two_thetas_obs, self.miller_indices, self.wavelength,
        unit_cell=uctbx.unit_cell(iter(x)))
    except KeyboardInterrupt: raise
    except Exception as e:
      print("FunctionalException:", str(e))
      raise scitbx.minimizers.FunctionalException
    if (len(self.functionals) != 0 and result > 2 * self.functionals[0]):
      print("FunctionalException: greater than 2 * initial")
      raise scitbx.minimizers.FunctionalException
    if (0):
      print("functional result:", result)
    self.functionals.append(result)
    return result

  def gradients(self, x):
    self.number_of_gradient_evaluations += 1
    return gradients(
      self.two_thetas_obs, self.miller_indices, self.wavelength,
      unit_cell=uctbx.unit_cell(iter(x)))

  def hessian(self, x):
    result = hessian(
      self.two_thetas_obs, self.miller_indices, self.wavelength,
      unit_cell=uctbx.unit_cell(iter(x)))
    if (self.mode == 1):
      d = result.matrix_diagonal()
      result = flex.double(flex.grid(6,6), 0)
      result.matrix_diagonal_set_in_place(diagonal=d)
    if (0):
      from scitbx import matrix
      print("hessian:")
      print(matrix.sqr(result))
    return result

  def run_lbfgs_raw(self,
        unit_cell,
        diagco,
        use_hessian,
        lbfgs_impl_switch):
    assert use_hessian in [0,1,2]
    n = 6
    m = 5
    x = flex.double(unit_cell.parameters())
    x_last = x.deep_copy()
    diag = flex.double(n, 0)
    iprint = [1, 0]
    eps = 1.0e-5
    xtol = 1.0e-16
    size_w = n*(2*m+1)+2*m
    w = flex.double(size_w)
    def diag_from_hessian():
      h = self.hessian(x)
      d = h.matrix_diagonal()
      assert d.all_gt(0)
      return 1/d
    iflag = 0
    while True:
      assert iflag in [0,1,2,100]
      if (iflag in [0,1]):
        try: f = self.functional(x=x)
        except scitbx.minimizers.FunctionalException: return x_last
        x_last = x.deep_copy()
        g = self.gradients(x=x)
      if (iflag == 0):
        if (diagco == 0):
          diag = flex.double(n, 0)
        elif (use_hessian == 0):
          diag = flex.double(n, 1)
        else:
          diag = diag_from_hessian()
          if (use_hessian == 1):
            diag0 = diag.deep_copy()
      elif (iflag == 2):
        if (use_hessian == 0):
          diag = flex.double(n, 1)
        elif (use_hessian == 1):
          diag = diag0.deep_copy()
        else:
          diag = diag_from_hessian()
      iflag = [scitbx.lbfgs.raw_reference,
               scitbx.lbfgs.raw][lbfgs_impl_switch](
        n=n, m=m, x=x, f=f, g=g, diagco=diagco, diag=diag,
        iprint=iprint, eps=eps, xtol=xtol, w=w, iflag=iflag)
      if (iflag <= 0): break
    return x

  def run_lbfgsb(self, unit_cell, iprint=1):
    l = flex.double([1,1,1,50,50,50])
    u = flex.double(6, 0)
    nbd = flex.int(6, 1) # all x have lower bound
    minimizer = scitbx.lbfgsb.minimizer(
      n=6,
      m=5,
      l=l,
      u=u,
      nbd=nbd,
      factr=1.0e+7,
      pgtol=1.0e-5,
      iprint=iprint)
    x = flex.double(unit_cell.parameters())
    f = 0
    g = flex.double(6, 0)
    while True:
      if (minimizer.process(x, f, g)):
        f = self.functional(x=x)
        g = self.gradients(x=x)
      elif (minimizer.is_terminated()):
        break
    return x

def show_fit(two_thetas_obs, miller_indices, wavelength, unit_cell):
  two_thetas_calc = unit_cell.two_theta(miller_indices, wavelength, deg=True)
  for h,o,c in zip(miller_indices, two_thetas_obs, two_thetas_calc):
    print("(%2d, %2d, %2d)" % h, "%6.2f - %6.2f = %6.2f" % (o, c, o-c))
  print()

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

def run(args):
  two_thetas_obs = flex.double()
  miller_indices = flex.miller_index()
  for line in two_theta_and_index_list:
    fields = line.split()
    assert len(fields) == 4
    two_thetas_obs.append(float(fields[0]))
    miller_indices.append([int(s) for s in fields[1:]])

  wavelength = wavelengths.characteristic("CU").as_angstrom()
  unit_cell_start = uctbx.unit_cell((10,10,10,80,100,80))
  show_fit(
    two_thetas_obs, miller_indices, wavelength, unit_cell_start)

  refined_accu = []
  assert len(args) in [0,1]
  if (len(args) != 0):
    arg = args[0]
    if (arg.find(",") > 0):
      modes = eval("["+arg+"]")
    elif (arg.find(":") > 0):
      modes = eval("range("+arg.replace(":",",")+"+1)")
    else:
      modes = [eval(arg)]
  else:
    modes = list(range(8))
  print("modes:", modes)
  print()
  for mode in modes:
    refined = refinery(
      two_thetas_obs, miller_indices, wavelength, unit_cell_start, mode=mode)
    refined_accu.append(refined)
    print()

  p = open("tmp.xy", "w")
  print("@with g0", file=p)
  print('@ title "%s"' % "\\n".join([
    platform.platform(),
    platform.node()]), file=p)
  for i,r in enumerate(refined_accu):
    print('@ s%d legend "%s"' % (i, r.plot_legend), file=p)
    print('@ s%d symbol 1' % i, file=p)
  for refined in refined_accu:
    for x,y in enumerate(refined.functionals):
      if (x > 15): break
      print(x,y, file=p)
    print("&", file=p)
  del p

  if (0):
    show_fit(
      two_thetas_obs, miller_indices, wavelength, refined.unit_cell())

  print(refined.unit_cell())
  print()
  print(libtbx.utils.format_cpu_times())

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
