from __future__ import absolute_import, division, print_function
class shelxl_wght_ls(object):

  def __init__(O, f_obs, i_obs, i_sig, i_calc=None, f_calc=None, wa=0.1, wb=0):
    assert [i_calc, f_calc].count(None) == 1
    if (i_calc is None):
      from cctbx.array_family import flex
      i_calc = flex.norm(f_calc)
    from cctbx import xray
    raw = xray.targets_shelxl_wght_ls_kwt_b_dv(
      f_obs=f_obs,
      i_obs=i_obs,
      i_sig=i_sig,
      ic=i_calc,
      wa=wa,
      wb=wb)
    assert len(raw) == 5
    O.scale_factor, \
    O.weights, \
    O.target, \
    O.i_gradients, \
    O.i_curvatures = raw
    if (f_calc is None):
      O.f_gradients = None
      O.f_hessians = None
    else:
      g = O.i_gradients
      c = O.i_curvatures
      O.f_gradients = 2 * g * f_calc
      a = flex.real(f_calc)
      b = flex.imag(f_calc)
      aa = 2 * g + 4 * a * a * c
      bb = 2 * g + 4 * b * b * c
      ab =         4 * a * b * c
      O.f_hessians = flex.vec3_double(aa, bb, ab)

  def target_work(O):
    return O.target

  def gradients_work(O):
    return O.f_gradients

  def hessians_work(O):
    return O.f_hessians
