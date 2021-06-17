from __future__ import absolute_import, division, print_function
from cctbx.geometry_restraints import angle, angle_delta_deg
from cctbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
import sys
from six.moves import range

def derivs_fd(a, order, eps=1.e-6):
  result = flex.vec3_double()
  sites0 = a.sites
  sites = [list(site) for site in sites0]
  for i_site in range(3):
    ds = []
    for i_dim in range(3):
      samples = []
      for signed_eps in [eps, -eps]:
        sites[i_site][i_dim] = sites0[i_site][i_dim] + signed_eps
        a_eps = angle(sites=sites, angle_ideal=a.angle_ideal, weight=a.weight)
        if (order == 1):
          samples.append(a_eps.residual())
        elif (order == 2):
          samples.append(a_eps.gradients()[i_site][i_dim])
        else:
          raise RuntimeError
      sites[i_site][i_dim] = sites0[i_site][i_dim]
      ds.append((samples[0]-samples[1])/(2*eps))
    result.append(tuple(ds))
  return result

def compare_derivs(out, ana, fin, expect_failure, eps=1.e-6):
  s = 1 / max(1, flex.max(flex.abs(ana.as_double())))
  print("  ana:", list(ana*s), file=out)
  print("  fin:", list(fin*s), file=out)
  if (not expect_failure):
    assert approx_equal(ana*s, fin*s, eps=eps)
  else:
    # This is to suppress output because it confuses testing functionality.
    # Apparently, it was never working and development was stopped long time ago.
    # To continue development/debugging - uncomment approx_equal
    pass
    # approx_equal(ana*s, fin*s, eps=eps)
  print(file=out)

def check_derivs(out, a, expect_failure=False):
  print("sites:", a.sites, file=out)
  print(file=out)
  gc = a.grads_and_curvs()
  print("grads:", a.sites, file=out)
  g_fd = derivs_fd(a=a, order=1)
  compare_derivs(out=out, ana=gc[:3], fin=g_fd, expect_failure=expect_failure)
  print("curvs:", a.sites, file=out)
  c_fd = derivs_fd(a=a, order=2)
  compare_derivs(out=out, ana=gc[3:], fin=c_fd, expect_failure=expect_failure)
  print(file=out)

def write_plots(method, rot_scale=2, rot_step=1, c_truncate_at=-2e4):
  assert method in ["ana", "fin"]
  plot_file_names = []
  site0 = (-1,0,0)
  site1 = (0,0,0)
  axis = matrix.col((0,0,1))
  for angle_ideal in [180, 120]:
    def init_plots():
      result = []
      for j in range(3):
        result.append([])
      return result
    g_plots = init_plots()
    c_plots = init_plots()
    for rot_deg_sc in range(90*rot_scale, 270*rot_scale+rot_step, rot_step):
      rot_deg = rot_deg_sc / rot_scale
      r = axis.axis_and_angle_as_r3_rotation_matrix(angle=rot_deg, deg=True)
      a = angle(
        sites=[site0, site1, r*site0],
        angle_ideal=angle_ideal,
        weight=1)
      if (method == "ana"):
        gc = a.grads_and_curvs()
      else:
        gc = derivs_fd(a=a, order=1)
        gc.extend(derivs_fd(a=a, order=2))
      for j in range(3):
        g_plots[j].append((rot_deg, gc[2][j]))
      for j in range(3):
        c_plots[j].append((rot_deg, gc[5][j]))
    def write(deriv, plots):
      file_name = "angle_%d_%s_%s.xy" % (angle_ideal, deriv, method)
      plot_file_names.append(file_name)
      with open(file_name, "w") as f:
        print("@with g0", file=f)
        print('@ title "%s ideal=%d method=%s"' % (
          deriv, angle_ideal, method), file=f)
        for j in range(3):
          print('@ s%d legend "%s"' % (j, "xyz"[j]), file=f)
        for plot in plots:
          for x,y in plot:
            if (deriv == "curv" and y < c_truncate_at): y = c_truncate_at
            print(x,y, file=f)
          print("&", file=f)
    write(deriv="grad", plots=g_plots)
    write(deriv="curv", plots=c_plots)
  return plot_file_names

def run(args):
  assert args in [[], ["--verbose"]]
  if (len(args) == 0):
    out = null_out()
  else:
    out = sys.stdout
  #
  mt = flex.mersenne_twister(seed=0)
  #
  for i_trial in range(10):
    l0 = mt.random_double() + 0.5
    l1 = mt.random_double() + 0.5
    l2 = mt.random_double() + 0.5
    angle_model = mt.random_double() * 178 + 1 \
                + 180 * (mt.random_size_t() % 3 - 1)
    v = matrix.col(mt.random_double_point_on_sphere())
    axis = v.ortho()
    site1 = v * l1
    site0 = site1 + v * l0
    r = axis.axis_and_angle_as_r3_rotation_matrix(angle=angle_model, deg=True)
    site2 = site1 + (r * v) * l2
    a = angle(
      sites=[site0, site1, site2],
      angle_ideal=mt.random_double() * 720 - 360,
      weight=mt.random_double() * 10 + 0.1)
    assert approx_equal(min(
      abs(angle_delta_deg(angle_1=a.angle_model, angle_2= angle_model)),
      abs(angle_delta_deg(angle_1=a.angle_model, angle_2=-angle_model))), 0)
    check_derivs(out=out, a=a)
  #
  for site2 in [(0,2.3,0), (0,0,2.5)]:
    perm = flex.size_t([0,1,2])
    while True:
      a = angle(
        sites=tuple(flex.vec3_double(
          [(1.2,0,0), (0,0,0), site2]).select(perm)),
        angle_ideal=mt.random_double() * 720 - 360,
        weight=mt.random_double() * 10 + 0.1)
      check_derivs(out=out, a=a)
      if (not perm.next_permutation()):
        break
  #
  for site0 in [(1,0,0),(0,1,0),(0,0,1),(1,1,1)]:
    perm = flex.size_t([0,1,2])
    while True:
      a = angle(
        sites=tuple(flex.vec3_double(
          [site0, (0,0,0), -matrix.col(site0)]).select(perm)),
        angle_ideal=180,
        weight=1.3)
      check_derivs(out=out, a=a, expect_failure=True)
      if (not perm.next_permutation()):
        break
  #
  plot_file_names = []
  for method in ["ana", "fin"]:
    plot_file_names.extend(write_plots(method=method))
  f = open("angle_xy_as_pdf_commands", "w")
  for file_name in plot_file_names:
    print("ppdf %s > %s" % (file_name, file_name.replace(".xy", ".pdf")), file=f)
  f.close()
  #
  print("OK")

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
