from __future__ import division
from cctbx.omz import bfgs
from cctbx import xray
import cctbx.xray.targets
from cctbx.array_family import flex
from scitbx import matrix
import libtbx.phil
from itertools import count
from math import pi, atan2
import sys

def delta_estimation_minus_cos(limit, grad, curv):
  return limit/pi * atan2(pi/limit*grad, curv)

class shift_limit_pair(object):

  __slots__ = ["neg", "pos", "curv_est_factor"]

  def __init__(O, neg, pos=None, curv_est_factor=0.1):
    assert neg > 0
    if (pos is None):
      pos = neg
    else:
      assert pos > 0
    assert curv_est_factor > 0
    assert curv_est_factor <= 1
    O.neg = neg
    O.pos = pos
    O.curv_est_factor = curv_est_factor

  def get(O, grad):
    if (grad < 0): return O.neg
    return O.pos

  def curv_est_shift(O, grad):
    if (grad < 0):
      return O.neg * curv_est_factor
    return -O.pos * curv_est_factor

class dynamic_shift_limit_site(object):
  """neg, pos estimation: xrefine_pot first_negative_after_min(curv)
  neg = pos = d_min / 2
  """

  __slots__ = ["width"]

  def __init__(O, width):
    O.width = width

  def pair(O, x):
    return shift_limit_pair(neg=O.width)

class dynamic_shift_limit_u_iso(object):

  __slots__ = ["d_min"]

  def __init__(O, d_min):
    O.d_min = d_min

  def pair(O, x):
    """neg, pos estimation: xrefine_pot first_negative(curv)
    d_min y_interc slope
      0.5   0.04    1.6
      1     0.06    1.5
      2     0.2     1.4
      3     0.4     1.3
      4     0.6     1.2
      5     1.4     1.1
    """
    u_iso = x
    y_interc = 0.05 * O.d_min**2
    slope = max(1., 1.6 - 0.1 * O.d_min)
    slope_attenuation = 0.5 # XXX to be optimized
    pos = slope_attenuation * slope * max(0., u_iso) + y_interc
    neg = max(0.01, u_iso)
    return shift_limit_pair(neg=neg, pos=pos)

class xfgc_info(object):

  __slots__ = ["x", "funcl", "grads", "curvs", "is_iterate", "approx_quads"]

  def __init__(k, work_obj, is_iterate):
    k.x = work_obj.x.deep_copy()
    k.funcl = work_obj.funcl
    k.grads = work_obj.grads
    k.curvs = work_obj.curvs
    k.is_iterate = is_iterate
    k.approx_quads = None

  def check_strong_wolfe(k, l, pk, ak, c1, c2):
    assert l is not k
    fk = k.funcl
    fl = l.funcl
    gk = k.grads
    gl = l.grads
    decr_cond = (fl <= fk + c1 * ak * gk.dot(pk))
    curv_cond = (abs(gl.dot(pk)) <= c2 * abs(gk.dot(pk)))
    if (0): print "CHECK Wolfe decr curv", decr_cond, curv_cond
    if (1): # verify that curv cond can be reformulated using sk
      sk = l.x - k.x
      curv_cond_sk = (abs(gl.dot(sk)) <= c2 * abs(gk.dot(sk)))
      assert curv_cond_sk == curv_cond
    return decr_cond and curv_cond

def sign0(x):
  if (x < 0): return -1
  if (x > 0): return 1
  return 0

class refinement(object):

  def __init__(O,
        i_obs,
        f_obs,
        xray_structure,
        params,
        reference_structure,
        expected_n_refinable_parameters=None,
        plot_samples_id=None):
    O.params = params
    O.i_obs = i_obs
    O.f_obs = f_obs
    O.xray_structure = xray_structure
    O.calc_weights_if_needed()
    O.reference_structure = reference_structure
    O.grads = None
    O.curvs = None
    O.grads_mean_sq = None
    O.show_rms_info()
    if (O.reference_structure is None):
      O.x_reference = None
    else:
      O.pack_variables(xray_structure=O.reference_structure)
      O.x_reference = O.x
    O.pack_variables()
    print "Number of variables:", O.x.size()
    if (expected_n_refinable_parameters is not None):
      assert O.x.size() == expected_n_refinable_parameters
    #
    O.plot_samples_id = plot_samples_id
    if (len(O.params.plot_samples.stages) != 0):
      ix = O.params.plot_samples.ix
      if (ix is None or ix < 0 or ix >= O.x.size()):
        raise RuntimeError("Out of range: plot_samples.ix = %s" % str(ix))
    #
    O.xfgc_infos = []
    O.update_fgc(is_iterate=True)
    O.pseudo_curvs = None
    O.pseudo_curvs_i_info = None
    O.termination_remark = " UNDEFINED"
    #
    O.plot_samples("initial")
    #
    if (O.params.use_classic_lbfgs):
      O.classic_lbfgs()
    elif (O.params.use_lbfgs_emulation):
      O.lbfgs_emulation()
    else:
      O.developmental_algorithms()
    print "Number of iterations, evaluations: %d %d%s" % (
      O.i_step+1, len(O.xfgc_infos), O.termination_remark)
    #
    O.plot_samples("final")

  def plot_samples(O, stage):
    p = O.params.plot_samples
    if (stage not in p.stages):
      return
    ix = p.ix
    i_sc, x_type = O.x_info[ix]
    xs = O.xray_structure
    site_inp = xs.scatterers()[i_sc].site
    ys = []
    xyv = []
    x_inp = O.x[ix]
    try:
      if (x_type == "u"):
        assert p.u_min < p.u_max
        assert p.u_steps > 0
        for i_step in xrange(p.u_steps+1):
          u = p.u_min + i_step / p.u_steps * (p.u_max - p.u_min)
          O.x[ix] = u
          ls = O.__get_tg()
          y = ls.target_work()
          ys.append(y)
          xyv.append((u,y,O.x[ix]))
      else:
        assert p.x_radius > 0
        assert p.x_steps > 0
        O.x[ix] = x_inp + 1
        O.__unpack_variables()
        dist = xs.unit_cell().distance(xs.scatterers()[i_sc].site, site_inp)
        assert dist != 0
        x_scale = p.x_radius / dist
        for i_step in xrange(-p.x_steps//2, p.x_steps//2+1):
          x = i_step / p.x_steps * 2 * x_scale
          O.x[ix] = x_inp + x
          ls = O.__get_tg()
          y = ls.target_work()
          ys.append(y)
          dist = xs.unit_cell().distance(xs.scatterers()[i_sc].site, site_inp)
          if (i_step < 0): dist *= -1
          xyv.append((dist,y,O.x[ix]))
    finally:
      O.x[ix] = x_inp
      O.__unpack_variables()
    def build_info_str():
      sc = xs.scatterers()[i_sc]
      ss = xs.site_symmetry_table().get(i_sc)
      return "%s|%03d|%03d|%s|%s|occ=%.2f|%s|%s" % (
        O.plot_samples_id,
        ix,
        i_sc,
        sc.label,
        sc.scattering_type,
        sc.occupancy,
        ss.special_op_simplified(),
        x_type)
    info_str = build_info_str()
    base_name_plot_files = "%s_%03d_%s" % (p.file_prefix, ix, stage)
    def write_xyv():
      if (p.file_prefix is None): return
      f = open(base_name_plot_files+".xy", "w")
      print >> f, "# %s" % info_str
      print >> f, "# %s" % str(xs.unit_cell())
      print >> f, "# %s" % str(xs.space_group_info().symbol_and_number())
      for x,y,v in xyv:
        print >> f, x, y, v
    write_xyv()
    from libtbx import pyplot
    fig = pyplot.figure()
    ax = fig.add_subplot(1, 1, 1)
    x,y,v = zip(*xyv)
    ax.plot(x,y, "r-")
    ax.plot([x_inp, x_inp], [min(ys), max(ys)], "k--")
    if (O.x_reference is not None):
      x = O.x_reference[ix]
      ax.plot([x, x], [min(ys), max(ys)], "r--")
    ax.set_title(info_str, fontsize=12)
    ax.axis([xyv[0][0], xyv[-1][0], None, None])
    def write_pdf():
      if (p.file_prefix is None): return
      fig.savefig(base_name_plot_files+".pdf", bbox_inches="tight")
    write_pdf()
    if (p.gui):
      pyplot.show()

  def classic_lbfgs(O):
    import scitbx.lbfgs
    O.i_step = 0
    scitbx.lbfgs.run(target_evaluator=O)

  def compute_functional_and_gradients(O):
    O.update_fgc()
    return O.funcl, O.grads

  def callback_after_step(O, minimizer):
    O.update_fgc(is_iterate=True)
    print "%4d: %s" % (O.i_step+1, O.format_rms_info())
    sys.stdout.flush()
    if (O.grads_mean_sq < O.params.grads_mean_sq_threshold):
      O.termination_remark = ""
      return True
    if (O.i_step+1 == O.params.iteration_limit):
      O.termination_remark = " (iteration limit reached)"
      return True
    O.i_step += 1

  def lbfgs_emulation(O, memory_range=5):
    assert len(O.xfgc_infos) == 1
    for O.i_step in xrange(O.params.iteration_limit):
      if (O.i_step == 0):
        dests = -O.grads
        stp = 1 / O.grads.norm()
      else:
        active_infos = O.get_active_infos()
        assert len(active_infos) > 1
        if (memory_range is not None):
          active_infos = active_infos[-(memory_range+1):]
        memory = O.build_bfgs_memory(active_infos=active_infos)
        assert memory is not None
        k_1 = active_infos[-1]
        k_2 = active_infos[-2]
        gamma = bfgs.h0_scaling(
          sk=k_1.x-k_2.x,
          yk=k_1.grads-k_2.grads)
        hk0 = flex.double(O.x.size(), gamma)
        dests = -bfgs.hg_two_loop_recursion(
          memory=memory, hk0=hk0, gk=O.grads)
        stp = 1
      stp = O.line_search(dests, stp=stp)
      assert stp is not None
      O.update_fgc(is_iterate=True)
      print "%4d: %s" % (O.i_step+1, O.format_rms_info())
      sys.stdout.flush()
      if (O.grads_mean_sq < O.params.grads_mean_sq_threshold):
        O.termination_remark = ""
        break
    else:
      O.termination_remark = " (iteration limit reached)"

  def developmental_algorithms(O):
    for O.i_step in xrange(O.params.iteration_limit):
      O.compute_step()
      s = "%4d: %s" % (O.i_step+1, O.format_rms_info())
      if (O.aq_sel_size is not None):
        s += " aq(%d, %d)" % (O.aq_sel_size, O.aq_n_used)
      print s
      sys.stdout.flush()
      if (O.grads_mean_sq < O.params.grads_mean_sq_threshold):
        O.termination_remark = ""
        break
    else:
      O.termination_remark = " (iteration limit reached)"

  def pack_variables(O, xray_structure=None):
    if (xray_structure is None):
      xray_structure = O.xray_structure
    O.x = flex.double()
    O.x_info = []
    O.gact_indices = flex.size_t()
    O.dynamic_shift_limits = []
    site_limits = [0.15/p for p in xray_structure.unit_cell().parameters()[:3]]
    d_min = O.f_obs.d_min()
    i_all = 0
    sstab = xray_structure.site_symmetry_table()
    for i_sc,sc in enumerate(xray_structure.scatterers()):
      assert sc.flags.use_u_iso()
      assert not sc.flags.use_u_aniso()
      #
      site_symmetry = sstab.get(i_sc)
      if (site_symmetry.is_point_group_1()):
        p = sc.site
        l = site_limits
      else:
        p = site_symmetry.site_constraints().independent_params(
          all_params=sc.site)
        l = site_symmetry.site_constraints().independent_params(
          all_params=site_limits)
      O.x.extend(flex.double(p))
      O.x_info.extend([(i_sc,"xyz"[i]) for i in xrange(len(p))])
      O.dynamic_shift_limits.extend(
        [dynamic_shift_limit_site(width=width) for width in l])
      for i in xrange(len(p)):
        O.gact_indices.append(i_all)
        i_all += 1
      #
      O.x.append(sc.u_iso)
      O.x_info.append((i_sc,"u"))
      O.dynamic_shift_limits.append(dynamic_shift_limit_u_iso(d_min=d_min))
      O.gact_indices.append(i_all)
      i_all += 1
      #
      i_all += 3 # occ, fp, fdp

  def __unpack_variables(O):
    ix = 0
    sstab = O.xray_structure.site_symmetry_table()
    for i_sc,sc in enumerate(O.xray_structure.scatterers()):
      site_symmetry = sstab.get(i_sc)
      if (site_symmetry.is_point_group_1()):
        sc.site = tuple(O.x[ix:ix+3])
        ix += 3
      else:
        constr = site_symmetry.site_constraints()
        np = constr.n_independent_params()
        sc.site = constr.all_params(independent_params=tuple(O.x[ix:ix+np]))
        ix += np
      sc.u_iso = O.x[ix]
      ix += 1
    assert ix == O.x.size()

  def show_dests(O, dests):
    from libtbx.str_utils import format_value
    ix = 0
    sstab = O.xray_structure.site_symmetry_table()
    for i_sc,sc in enumerate(O.xray_structure.scatterers()):
      site_symmetry = sstab.get(i_sc)
      if (site_symmetry.is_point_group_1()):
        vals = list(dests[ix:ix+3])
        ix += 3
      else:
        constr = site_symmetry.site_constraints()
        np = constr.n_independent_params()
        vals = list(dests[ix:ix+np]) + [None] * (3-np)
        ix += np
      vals.append(dests[ix])
      ix += 1
      print " ".join([format_value("%15.6f", v) for v in vals])
    assert ix == O.x.size()

  def get_f_calc(O):
    p = O.params.f_calc_options
    return O.f_obs.structure_factors_from_scatterers(
      xray_structure=O.xray_structure,
      algorithm=p.algorithm,
      cos_sin_table=p.cos_sin_table).f_calc()

  def r1_factor(O):
    from libtbx import Auto
    return O.f_obs.r1_factor(
      other=O.get_f_calc(), scale_factor=Auto, assume_index_matching=True)

  def calc_weights_if_needed(O):
    O.weights = None
    if (O.params.i_obs_weights == "shelxl_wght_once"):
      assert O.params.target_obs_type == "I"
      f_calc = O.get_f_calc()
      i_calc = flex.norm(f_calc.data())
      from cctbx.xray.targets.tst_shelxl_wght_ls import calc_k, calc_w
      k = calc_k(f_obs=O.f_obs.data(), i_calc=i_calc)
      assert O.i_obs.sigmas().all_ge(0.01)
      O.weights = calc_w(
        wa=0.1,
        wb=0,
        i_obs=O.i_obs.data(),
        i_sig=O.i_obs.sigmas(),
        i_calc=i_calc,
        k=k)

  def __get_tg(O):
    O.__unpack_variables()
    if (O.params.target_obs_type == "F"):
      obs = O.f_obs
    else:
      obs = O.i_obs
    if (O.params.target_type == "ls"):
      if (O.params.i_obs_weights in ["unit", "shelxl_wght_once"]):
        return xray.targets_least_squares(
          compute_scale_using_all_data=True,
          obs_type=O.params.target_obs_type,
          obs=obs.data(),
          weights=O.weights,
          r_free_flags=None,
          f_calc=O.get_f_calc().data(),
          derivatives_depth=2,
          scale_factor=O.params.f_calc_scale_factor)
      if (O.params.i_obs_weights != "shelxl_wght"):
        raise RuntimeError(
          "Unknown: i_obs_weights = %s" % O.params.i_obs_weights)
      assert O.params.target_obs_type == "I"
      return xray.targets.shelxl_wght_ls(
        f_obs=O.f_obs.data(),
        i_obs=O.i_obs.data(),
        i_sig=O.i_obs.sigmas(),
        f_calc=O.get_f_calc().data(),
        i_calc=None,
        wa=0.1,
        wb=0)
    elif (O.params.target_type == "cc"):
      return xray.targets_correlation(
        obs_type=O.params.target_obs_type,
        obs=obs.data(),
        weights=O.weights,
        r_free_flags=None,
        f_calc=O.get_f_calc().data(),
        derivatives_depth=2)
    elif (O.params.target_type == "r1"):
      assert O.params.target_obs_type == "F"
      from cctbx.xray.targets import r1
      return r1.target(
        f_obs=O.f_obs.data(),
        f_calc=O.get_f_calc().data())
    raise RuntimeError("Unknown target_type.")

  def update_fgc(O, is_iterate=False):
    if (len(O.xfgc_infos) != 0):
      prev_xfgc = O.xfgc_infos[-1]
      if (prev_xfgc.x.all_eq(O.x)):
        if (not prev_xfgc.is_iterate):
          prev_xfgc.is_iterate = is_iterate
        return
    tg = O.__get_tg()
    assert tg.target_work() is not None
    gact = O.xray_structure.grads_and_curvs_target_simple(
      miller_indices=O.f_obs.indices(),
      da_db=tg.gradients_work(),
      daa_dbb_dab=tg.hessians_work())
    O.funcl = tg.target_work()
    O.grads = gact.grads.select(O.gact_indices)
    O.curvs = gact.curvs.select(O.gact_indices)
    O.grads_mean_sq = flex.mean_sq(O.grads)
    O.xfgc_infos.append(xfgc_info(work_obj=O, is_iterate=is_iterate))

  def get_active_infos(O, i_start=0):
    result = []
    for info in O.xfgc_infos[i_start:]:
      if (info.is_iterate):
        result.append(info)
    return result

  def build_bfgs_memory(O, active_infos):
    result = []
    for iinfo in xrange(len(active_infos)-1):
      k = active_infos[iinfo]
      l = active_infos[iinfo+1]
      m = bfgs.memory_element(s=l.x-k.x, y=l.grads-k.grads)
      if (m.rho is None):
        return None
      result.append(m)
    return result

  def update_dests_using_bfgs_formula(O, dests):
    O.aq_sel_size = -2
    O.aq_n_used = -2
    if (O.params.bfgs_estimate_limit_factor <= 0):
      return
    aq_sel = flex.size_t()
    aq_sel_size_start = 0
    iinfo_active = []
    for iinfo in xrange(len(O.xfgc_infos)-1,-1,-1):
      info = O.xfgc_infos[iinfo]
      if (info.is_iterate):
        if (aq_sel_size_start == 0):
          aq_sel = info.approx_quads
          aq_sel_size_start = aq_sel.size()
          if (aq_sel_size_start < 2):
            return
        else:
          next_aq_sel = aq_sel.intersection(other=info.approx_quads)
          if (    next_aq_sel.size() < aq_sel_size_start * 0.9
              and len(iinfo_active) > 1):
            break
          aq_sel = next_aq_sel
        iinfo_active.append(iinfo)
    iinfo_active.sort()
    O.aq_sel_size = aq_sel.size()
    if (len(iinfo_active) < 2 or O.aq_sel_size < 2):
      return
    O.aq_n_used = -1
    assert iinfo_active[-1] == len(O.xfgc_infos)-1
    curvs = O.xfgc_infos[iinfo_active[-1]].curvs.select(aq_sel)
    assert curvs.all_gt(0)
    hk0 = 1 / curvs
    memory = []
    for iinfo in iinfo_active[:-1]:
      k = O.xfgc_infos[iinfo]
      l = O.xfgc_infos[iinfo+1]
      xk = k.x.select(aq_sel)
      xl = l.x.select(aq_sel)
      gk = k.grads.select(aq_sel)
      gl = l.grads.select(aq_sel)
      m = bfgs.memory_element(s=xl-xk, y=gl-gk)
      gks = gk.dot(m.s)
      gls = gl.dot(m.s)
      wolfe_curv_cond = (gls >= 0.9 * gks)
        # Nocedal & Wright (1999) Equation 3.7b
        # reformulated using sk instead of pk
      if (not wolfe_curv_cond):
        return
      if (m.rho is None):
        print "Warning: rho <= 0"
        return
      memory.append(m)
    aq_dests = -bfgs.hg_two_loop_recursion(
      memory=memory, hk0=hk0, gk=O.xfgc_infos[-1].grads.select(aq_sel))
    O.aq_n_used = 0
    for aq_dest,ix in zip(aq_dests, aq_sel):
      dsl = O.dynamic_shift_limits[ix]
      limit = dsl.pair(x=O.x[ix]).get(grad=O.grads[ix])
      if (abs(aq_dest) <= O.params.bfgs_estimate_limit_factor * limit):
        dests[ix] = aq_dest
        O.aq_n_used += 1

  def approx_curvs(O, shift_limit_factor=0.1):
    x_on_entry = O.x
    shifts = flex.double()
    for ix,dsl,g in zip(count(), O.dynamic_shift_limits, O.grads):
      shift = dsl.pair(x=O.x[ix]).get(grad=g) * shift_limit_factor
      if (g > 0):
        shift *= -1
      assert shift != 0
      shifts.append(shift)
    if (0): print "shifts:", list(shifts)
    grads_z = O.grads
    O.x = x_on_entry + shifts
    O.update_fgc()
    grads_p = O.grads
    O.x = x_on_entry - shifts
    O.update_fgc()
    grads_m = O.grads
    O.x = x_on_entry
    O.update_fgc()
    O.xfgc_infos.pop()
    O.xfgc_infos.pop()
    O.xfgc_infos.pop()
    apprx = (grads_p - grads_m) / (2*shifts)
    if (0):
      print "curvs:", list(O.curvs)
      print "apprx:", list(apprx)
      flex.linear_correlation(O.curvs, apprx).show_summary()
    O.curvs = apprx

  def compute_step(O):
    if (O.params.use_curvs):
      O.compute_step_using_curvs()
    else:
      O.compute_step_just_grads()

  def compute_step_just_grads(O):
    inp_i_info = len(O.xfgc_infos) - 1
    inp_info = O.xfgc_infos[-1]
    limits = flex.double()
    for ix,dsl,g in zip(count(), O.dynamic_shift_limits, inp_info.grads):
      limits.append(dsl.pair(x=O.x[ix]).get(grad=g))
    assert limits.all_gt(0)
    def get_pseudo_curvs():
      ag_max = flex.max(flex.abs(inp_info.grads))
      assert ag_max != 0
      dests = (-inp_info.grads/ag_max) * (limits/2)
      assert flex.abs(dests).all_le(limits/2*(1+1e-6))
      assert (dests > 0).all_eq(inp_info.grads < 0)
      O.pseudo_curvs_i_info = inp_i_info
      return dests
    if (O.pseudo_curvs is None):
      dests = get_pseudo_curvs()
    else:
      active_infos = O.get_active_infos(O.pseudo_curvs_i_info)
      assert len(active_infos) > 1
      memory = O.build_bfgs_memory(active_infos=active_infos)
      if (memory is None):
        O.pseudo_curvs = None
        dests = get_pseudo_curvs()
      else:
        hk0 = 1 / O.pseudo_curvs
        dests = -bfgs.hg_two_loop_recursion(
          memory=memory, hk0=hk0, gk=inp_info.grads)
        madl = flex.max(flex.abs(dests / limits))
        if (madl > 1):
          print "madl:", madl
          dests *= (1/madl)
        assert flex.abs(dests).all_le(limits*(1+1e-6))
    dest_adj = O.line_search(dests, stpmax=2.0)
    print "dest_adj:", dest_adj
    if (dest_adj is not None):
      dests *= dest_adj
    elif (O.pseudo_curvs is not None):
      O.pseudo_curvs = None
      dests = get_pseudo_curvs()
      dest_adj = O.line_search(dests, stpmax=2.0)
      if (dest_adj is not None):
        dests *= dest_adj
    if (O.pseudo_curvs is None):
      assert (dests > 0).all_eq(inp_info.grads < 0)
      assert flex.abs(dests).all_le(limits*(1+1e-6))
      O.pseudo_curvs = -inp_info.grads / dests
      assert O.pseudo_curvs.all_gt(0)
    O.x = inp_info.x + dests
    O.update_fgc(is_iterate=True)
    O.aq_sel_size = None
    O.aq_n_used = None

  def compute_step_using_curvs(O):
    if (len(O.xfgc_infos) > 1 and O.params.use_gradient_flips):
      prev = O.xfgc_infos[-2]
    else:
      prev = None
    dests = flex.double()
    approx_quads = flex.size_t()
    if (O.params.try_approx_curvs):
      O.approx_curvs()
    for ix,dsl,g,c in zip(count(), O.dynamic_shift_limits, O.grads, O.curvs):
      limit = dsl.pair(x=O.x[ix]).get(grad=g)
      dest = None
      if (prev is not None):
        prev_g = prev.grads[ix]
        if (sign0(g) != sign0(prev_g)):
          x = O.x[ix]
          prev_x = prev.x[ix]
          xm = (g*prev_x - prev_g*x) / (g - prev_g)
          dest = xm - x
          if   (dest >  limit): dest =  limit
          elif (dest < -limit): dest = -limit
      if (dest is None):
        if (c > 0
              and O.params.approx_quad_limit_factor > 0
              and abs(g) < O.params.approx_quad_limit_factor * limit * c):
          dest = -g / c
          approx_quads.append(ix)
        else:
          dest = -delta_estimation_minus_cos(
            limit=limit, grad=g, curv=c)
      dests.append(dest)
    O.xfgc_infos[-1].approx_quads = approx_quads
    O.update_dests_using_bfgs_formula(dests)
    O.x_before_line_search = O.xfgc_infos[-1].x
    if (O.params.use_line_search):
      dest_adj = O.line_search(dests, stpmax=1.0)
      print "dest_adj:", dest_adj
      if (dest_adj is not None and dest_adj < 1):
        dests *= dest_adj
    if (O.params.show_dests): O.show_dests(dests)
    O.x = O.x_before_line_search + dests
    O.update_fgc(is_iterate=True)

  def line_search(O, dests, stp=1, stpmax=None):
    import scitbx.math
    k = O.xfgc_infos[-1]
    line_search = scitbx.math.line_search_more_thuente_1994()
    line_search.ftol = 1e-4
    line_search.gtol = 0.9
    if (stpmax is not None):
      assert stp <= stpmax
      line_search.stpmax = stpmax
    try:
      line_search.start(
        x=O.x,
        functional=O.funcl,
        gradients=O.grads,
        search_direction=dests,
        initial_estimate_of_satisfactory_step_length=stp)
    except RuntimeError, e:
      if (str(e) != "Search direction not descent."):
        raise
      return None
    assert line_search.info_code == -1
    O.update_fgc()
    while (line_search.info_code == -1):
      chk = k.check_strong_wolfe(
        l=O.xfgc_infos[-1],
        pk=dests,
        ak=line_search.stp,
        c1=line_search.ftol,
        c2=line_search.gtol)
      line_search.next(x=O.x, functional=O.funcl, gradients=O.grads)
      if (chk):
        assert line_search.info_code != -1
      elif (line_search.info_code == 5):
        assert line_search.info_meaning \
            == "The step is at the upper bound stpmax."
        return None
      elif (line_search.info_code == 6):
        assert line_search.info_meaning.startswith(
          "Rounding errors prevent further progress.")
        return None
      else:
        assert line_search.info_code == -1, (
          line_search.info_code, line_search.info_meaning)
      O.update_fgc()
    return line_search.stp

  def get_rms_info(O):
    if (O.reference_structure is None):
      return None
    xs = O.xray_structure
    rs = O.reference_structure
    xf = xs.sites_frac()
    rf = rs.sites_frac()
    cd = xf - rf
    # TODO: use scattering power as weights, move to method of xray.structure
    ave_csh = matrix.col(cd.mean())
    ave_csh_perp = matrix.col(xs.space_group_info()
      .subtract_continuous_allowed_origin_shifts(translation_frac=ave_csh))
    caosh_corr = ave_csh_perp - ave_csh
    ad = cd + caosh_corr
    ud = xs.scatterers().extract_u_iso() \
       - rs.scatterers().extract_u_iso()
    omx = xs.unit_cell().orthogonalization_matrix()
    O.crmsd = (omx * cd).rms_length()
    O.armsd = (omx * ad).rms_length()
    O.urmsd = flex.mean_sq(ud)**0.5
    if (O.params.show_distances_to_reference_structure):
      for sc, a, u in zip(xs.scatterers(), omx * ad, ud):
        print "    %-10s" % sc.label, \
          " ".join(["%6.3f" % v for v in a]), \
          "%6.3f" % u
    return (O.crmsd, O.armsd, O.urmsd)

  def format_rms_info(O):
    s = ""
    info = O.get_rms_info()
    if (info is not None):
      for r in info:
        s += " %5.3f" % r
    if (O.grads_mean_sq is not None):
      s += " f=%8.2e |g|=%8.2e" % (O.funcl, O.grads_mean_sq)
    s += " r1=%.3f" % O.r1_factor()
    return s

  def show_rms_info(O):
    s = O.format_rms_info()
    if (len(s) != 0):
       print " cRMSD aRMSD uRMSD"
       print s
       sys.stdout.flush()

def run_refinement(
      structure_ideal,
      structure_shake,
      params,
      i_obs=None,
      f_obs=None):
  assert (i_obs is None) == (f_obs is None)
  print "Ideal structure:"
  structure_ideal.show_summary().show_scatterers()
  print
  print "Modified structure:"
  structure_shake.show_summary().show_scatterers()
  print
  print "rms difference:", \
    structure_ideal.rms_difference(other=structure_shake)
  print
  sdt = params.show_distances_threshold
  if (sdt > 0):
    print "structure_shake inter-atomic distances:"
    structure_shake.show_distances(distance_cutoff=sdt)
    print
  if (f_obs is None):
    i_obs = structure_ideal.structure_factors(
      anomalous_flag=False,
      d_min=1,
      algorithm="direct",
      cos_sin_table=False).f_calc().intensities()
    f_obs = i_obs.array(data=flex.sqrt(i_obs.data()))
  return refinement(
    i_obs=i_obs,
    f_obs=f_obs,
    xray_structure=structure_shake,
    params=params,
    reference_structure=structure_ideal)

def get_master_phil(
      iteration_limit=50,
      show_distances_threshold=0,
      grads_mean_sq_threshold=1e-6,
      f_calc_options_algorithm="*direct fft",
      additional_phil_string=""):
  return libtbx.phil.parse("""
    general_positions_only = True
      .type = bool
    shake_sites_rmsd = 0.5
      .type = float
    shake_adp_spread = 20
      .type = float
    show_distances_threshold = %(show_distances_threshold)s
      .type = float
    target_type = *ls cc r1
      .type = choice
    target_obs_type = *F I
      .type = choice
    i_obs_weights = *unit shelxl_wght shelxl_wght_once
      .type = choice
    iteration_limit = %(iteration_limit)s
      .type = int
    grads_mean_sq_threshold = %(grads_mean_sq_threshold)s
      .type = float
    use_classic_lbfgs = False
      .type = bool
    use_lbfgs_emulation = False
      .type = bool
    use_curvs = True
      .type = bool
    approx_quad_limit_factor = 1.0
      .type = float
    bfgs_estimate_limit_factor = 1.0
      .type = float
    use_line_search = False
      .type = bool
    use_gradient_flips = True
      .type = bool
    try_approx_curvs = False
      .type = bool
    f_calc_scale_factor = 0
      .type = float
    show_dests = False
      .type = bool
    show_distances_to_reference_structure = False
      .type = bool
    f_calc_options {
      algorithm = %(f_calc_options_algorithm)s
        .type = choice
      cos_sin_table = False
        .type = bool
    }
    plot_samples {
      stages = initial final
        .type = choice(multi=True)
      ix = None
        .type = int
      x_radius = 5
        .type = float
      x_steps = 100
        .type = int
      u_min = -0.05
        .type = float
      u_max = 0.45
        .type = float
      u_steps = 100
        .type = int
      gui = True
        .type = bool
      file_prefix = None
        .type = str
    }
""" % vars() + additional_phil_string)
