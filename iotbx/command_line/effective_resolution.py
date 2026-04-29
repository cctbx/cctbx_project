"""Calculate effective resolution of a reflection file"""

from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.resolution

from iotbx import reflection_file_reader
from cctbx import maptbx
from cctbx.array_family import flex
from libtbx.utils import Sorry
import sys, math, time
from scitbx import regular_grid_on_unit_sphere
import six
from six.moves import range

def one_d_image_along_axis(n, step, uc_length):
  rho  = flex.double()
  dist = flex.double()
  r = 0
  while r < uc_length/2:
    rho_ = 0
    for n_key, n_value in six.iteritems(n):
      rho_ += n_value*math.cos(2*math.pi*r*n_key/uc_length)
    dist.append(r)
    rho.append(rho_)
    r+=step
  return dist, rho

def second_derivatives(rho, delta):
  rho_2nd = flex.double()
  for i in range(rho.size()):
    if(i>=1 and i<rho.size()-1): tau = (rho[i+1]+rho[i-1]-2*rho[i])/delta**2
    elif(i==0):                  tau = (rho[i+1]+rho[i+1]-2*rho[i])/delta**2
    else:                        tau = (rho[i+0]+rho[i-1]-2*rho[i])/delta**2
    rho_2nd.append(tau)
  result = flex.double()
  for i in range(rho_2nd.size()):
    rho_ave = 0
    span = [-2,-1,0,1,2]
    #span = [-5,-4,-3,-2,-1,0,1,2, 3,4,5]
    for j in span:
      ij = i+j
      if(ij<0):               ij=0
      if(ij>=rho_2nd.size()): ij=rho_2nd.size()-1
      rho_ave += rho_2nd[ij]
    rho_ave = rho_ave/len(span)
    result.append(rho_ave)
  return result

def compute_d_eff(r, rho_2nd):
  v0 = rho_2nd[0]
  for i in range(rho_2nd.size()):
    if(v0*rho_2nd[i]<0):
      return r[i]*2.5 # formulas (9)-(10)
  return None

def compute(miller_array, step_scale=0.0005):
  miller_array.show_comprehensive_summary(prefix="  ")
  step = miller_array.d_min()*step_scale
  #
  ma_p1 = miller_array.expand_to_p1()
  #
  n_h = {}
  n_k = {}
  n_l = {}
  indices = ma_p1.indices()
  for ind in indices:
    h,k,l = ind
    n_h.setdefault(h, flex.int()).append(1)
    n_k.setdefault(k, flex.int()).append(1)
    n_l.setdefault(l, flex.int()).append(1)
  def count(d):
    for k in d.keys():
      d[k] = d[k].size()
    return d
  n_h = count(n_h)
  n_k = count(n_k)
  n_l = count(n_l)
  # resolutions along axes
  a,b,c = miller_array.unit_cell().parameters()[:3]
  x, rho_x = one_d_image_along_axis(n=n_h, step=step, uc_length=a)
  y, rho_y = one_d_image_along_axis(n=n_k, step=step, uc_length=b)
  z, rho_z = one_d_image_along_axis(n=n_l, step=step, uc_length=c)
  # 2nd derivatives
  r2x = second_derivatives(rho=rho_x, delta=step)
  r2y = second_derivatives(rho=rho_y, delta=step)
  r2z = second_derivatives(rho=rho_z, delta=step)
  # effective resolution along axes
  d_eff_a = compute_d_eff(r=x, rho_2nd=r2x)
  d_eff_b = compute_d_eff(r=y, rho_2nd=r2y)
  d_eff_c = compute_d_eff(r=z, rho_2nd=r2z)
  print("  Effective resolution along axes a,b,c: %6.3f %6.3f %6.3f"%(
    d_eff_a, d_eff_b, d_eff_c))
  # all directions
  l = 0.8 * min(d_eff_a/2.5, d_eff_b/2.5, d_eff_c/2.5)
  r = 1.2 * max(d_eff_a/2.5, d_eff_b/2.5, d_eff_c/2.5)
  us = regular_grid_on_unit_sphere.rosca(m=9, hemisphere=True)
  d_effs = flex.double()
  o = maptbx.ft_analytical_1d_point_scatterer_at_origin(N=100000)
  for i, u in enumerate(us):
    o.compute(
      miller_indices=indices,
      step=step,
      left=l,
      right=r,
      u_frac=miller_array.unit_cell().fractionalize(u))
    dist, rho_ = o.distances(), o.rho()
    rho2 = second_derivatives(rho=rho_, delta=step)
    d_eff = compute_d_eff(r=dist, rho_2nd=rho2)
    d_effs.append(d_eff)
  print("  Effective resolution (min,max): %8.3f%8.3f"%(
    flex.min(d_effs), flex.max(d_effs)))

def run(args):
  if(len(args)!=1):
    raise Sorry("Reflection file expected.")
  reflection_file = reflection_file_reader.any_reflection_file(
    file_name = args[0])
  miller_arrays = reflection_file.as_miller_arrays(
    force_symmetry=True,
    merge_equivalents=False)
  if(miller_arrays is None):
    raise Sorry("Warning: unknown file format:", file_name)
  for ma in miller_arrays:
    if(type(ma.data()) == type(flex.double())):
      print("Processing data array with labels:", ma.info().label_string())
      compute(miller_array=ma)
      print()

if (__name__ == "__main__"):
  t0 = time.time()
  run(args=sys.argv[1:])
  print("Time: %8.3f"%(time.time()-t0))
