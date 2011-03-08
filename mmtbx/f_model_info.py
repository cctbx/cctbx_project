from libtbx import adopt_init_args
from libtbx.str_utils import format_value
import sys, re, string

def _scale_helper(num, den, selection=None, num_num=False):
  from cctbx.array_family import flex
  if (selection is not None):
    num = num.select(selection)
    den = den.select(selection)
  if (den.size() == 0):
    raise RuntimeError("No data for scale calculation.")
  denom = flex.sum(den*den)
  if (denom == 0):
    raise RuntimeError("Zero denominator in scale calculation.")
  if (num_num): return flex.sum(num*num) / denom
  return flex.sum(num*den) / denom

def n_as_s(format, value):
  if (value is None):
    return format_value(format=format, value=value)
  if (isinstance(value, (int, float))):
    return (format % value).strip()
  return [(format % v).strip() for v in value]

class resolution_bin(object):
  def __init__(self,
               i_bin         = None,
               d_range       = None,
               completeness  = None,
               alpha_work    = None,
               beta_work     = None,
               r_work        = None,
               r_free        = None,
               target_work   = None,
               target_free   = None,
               n_work        = None,
               n_free        = None,
               mean_f_obs    = None,
               fom_work      = None,
               scale_k1_work = None,
               pher_work     = None,
               pher_free     = None,
               sigmaa        = None):
    adopt_init_args(self, locals())

class info(object):
  def __init__(self,
               fmodel,
               free_reflections_per_bin = 140,
               max_number_of_bins = 30):
    from cctbx.array_family import flex
    mp = fmodel.mask_params
    self.target_name = fmodel.target_name
    if(self.target_name == "twin_lsq_f"):
      self.twin_fraction = fmodel.twin_fraction
      self.twin_law = fmodel.twin_law
    else:
      self.twin_fraction = None
      self.twin_law = None
    self.r_work = fmodel.r_work()
    self.r_free = fmodel.r_free()
    self.r_all = fmodel.r_all()
    self.target_work = fmodel.target_w()
    self.target_free = fmodel.target_t()
    self.overall_scale_k1 = fmodel.scale_k1()
    self.number_of_test_reflections = fmodel.f_calc_t().data().size()
    self.number_of_work_reflections = fmodel.f_calc_w().data().size()
    self.number_of_reflections = fmodel.f_obs().data().size()
    self.k_sol = fmodel.shell_k_sols()
    self.b_sol = fmodel.b_sol()
    self.b_cart = fmodel.b_cart()
    self.b_iso = fmodel.b_iso()
    self.mask_solvent_radius = mp.solvent_radius
    self.mask_shrink_radius = mp.shrink_truncation_radius
    self.mask_grid_step_factor = mp.grid_step_factor
    self.ml_phase_error = flex.mean(fmodel.phase_errors())
    self.ml_coordinate_error = fmodel.model_error_ml()
    if hasattr(fmodel, "sigmaa") :
      try :
        self.sigmaa = fmodel.sigmaa().sigmaa() # miller array
      except RuntimeError, e :
        self.sigmaa = None
    else :
      self.sigmaa = None
    self.d_max, self.d_min = fmodel.f_obs().resolution_range()
    self.completeness_in_range = fmodel.f_obs().completeness(d_max = self.d_max)
    self.completeness_d_min_inf = fmodel.f_obs().completeness()
    f_obs_6 = fmodel.f_obs().resolution_filter(d_min = 6)
    self.completeness_6_inf = f_obs_6.completeness()
    self.min_f_obs_over_sigma = fmodel.f_obs().min_f_over_sigma(
      return_none_if_zero_sigmas=True)
    self.sf_algorithm = fmodel.sfg_params.algorithm
    alpha_w, beta_w = fmodel.alpha_beta_w()
    self.alpha_work_min, self.alpha_work_max, self.alpha_work_mean = \
      alpha_w.data().min_max_mean().as_tuple()
    self.beta_work_min, self.beta_work_max, self.beta_work_mean = \
      beta_w.data().min_max_mean().as_tuple()
    self.fom_work_min, self.fom_work_max, self.fom_work_mean = \
      fmodel.figures_of_merit_work().min_max_mean().as_tuple()
    self.pher_work_min, self.pher_work_max, self.pher_work_mean = \
      fmodel.phase_errors_work().min_max_mean().as_tuple()
    self.pher_free_min, self.pher_free_max, self.pher_free_mean = \
      fmodel.phase_errors_test().min_max_mean().as_tuple()
    self.bins = self.statistics_in_resolution_bins(
      fmodel = fmodel,
      free_reflections_per_bin = free_reflections_per_bin,
      max_number_of_bins = max_number_of_bins)

  def statistics_in_resolution_bins(self, fmodel, free_reflections_per_bin,
                                    max_number_of_bins):
    from mmtbx import bulk_solvent
    from cctbx.array_family import flex
    if(self.target_name == "twin_lsq_f"):
      return fmodel.statistics_in_resolution_bins()
    result = []
    target_functor = fmodel.target_functor()
    target_result = target_functor(compute_gradients=False)
    tpr = target_result.target_per_reflection()
    if(tpr.size() != 0):
      tpr_w = tpr.select(fmodel.active_arrays.work_sel)
      tpr_t = tpr.select(fmodel.active_arrays.free_sel)
    fo_t = fmodel.f_obs_free()
    fc_t = fmodel.f_model_scaled_with_k1_t()
    fo_w = fmodel.f_obs_work()
    fc_w = fmodel.f_model_scaled_with_k1_w()
    alpha_w, beta_w = fmodel.alpha_beta_w()
    alpha_t, beta_t = fmodel.alpha_beta_t()
    pher_w = fmodel.phase_errors_work()
    pher_t = fmodel.phase_errors_test()
    fom = fmodel.figures_of_merit_work()
    fmodel.f_obs().setup_binner(n_bins=fmodel.determine_n_bins(
      free_reflections_per_bin=free_reflections_per_bin,
      max_n_bins=max_number_of_bins))
    fo_t.use_binning_of(fmodel.f_obs())
    fc_t.use_binning_of(fo_t)
    fo_w.use_binning_of(fo_t)
    fc_w.use_binning_of(fo_t)
    alpha_w.use_binning_of(fo_t)
    alpha_t.use_binning_of(fo_t)
    beta_w.use_binning_of(fo_t)
    beta_t.use_binning_of(fo_t)
    if hasattr(fmodel, "sigmaa") :
      try :
        sigmaa = fmodel.sigmaa().sigmaa()
      except RuntimeError, e :
        sigmaa = None
      else :
        sigmaa.use_binning_of(fo_t)
    else :
      sigmaa = None
    for i_bin in fo_t.binner().range_used():
      sel_t = fo_t.binner().selection(i_bin)
      sel_w = fo_w.binner().selection(i_bin)
      sel_all = fmodel.f_obs().binner().selection(i_bin)
      sel_fo_all = fmodel.f_obs().select(sel_all)
      sel_fo_t = fo_t.select(sel_t)
      sel_fc_t = fc_t.select(sel_t)
      sel_fo_w = fo_w.select(sel_w)
      sel_fc_w = fc_w.select(sel_w)
      if (tpr.size() == 0):
        sel_tpr_w = None
        sel_tpr_t = None
      else:
        denom_w = sel_fo_w.data().size()
        denom_t = sel_fo_t.data().size()
        if(denom_w != 0):
           sel_tpr_w = flex.sum(tpr_w.select(sel_w))/denom_w
        else:
           sel_tpr_w = flex.sum(tpr_w.select(sel_w))
        if(denom_t != 0):
           sel_tpr_t = flex.sum(tpr_t.select(sel_t))/denom_t
        else:
           sel_tpr_t = flex.sum(tpr_t.select(sel_t))
      d_max_,d_min_ = sel_fo_all.d_max_min()
      completeness = sel_fo_all.completeness(d_max = d_max_)
      d_range = fo_t.binner().bin_legend(
        i_bin=i_bin, show_bin_number=False, show_counts=False)
      s_fo_w_d = sel_fo_w.data()
      s_fc_w_d = sel_fc_w.data()
      assert s_fo_w_d.size() == s_fc_w_d.size()
      s_fc_w_d_a = flex.abs(s_fc_w_d)
      sigmaa_bin = None
      if (sigmaa is not None) :
        sigmaa_bin = flex.mean_default(sigmaa.select(sel_all).data(), None)
      if(s_fo_w_d.size() > 0):
        bin = resolution_bin(
          i_bin        = i_bin,
          d_range      = d_range,
          completeness = completeness,
          alpha_work   = flex.mean_default(alpha_w.select(sel_w).data(),None),
          beta_work    = flex.mean_default(beta_w.select(sel_w).data(),None),
          r_work       = bulk_solvent.r_factor(s_fo_w_d, s_fc_w_d, 1),
          r_free       = bulk_solvent.r_factor(sel_fo_t.data(), sel_fc_t.data(), 1),
          target_work  = sel_tpr_w,
          target_free  = sel_tpr_t,
          n_work       = sel_fo_w.data().size(),
          n_free       = sel_fo_t.data().size(),
          scale_k1_work= _scale_helper(num=s_fo_w_d, den=s_fc_w_d_a),
          mean_f_obs   = flex.mean_default(sel_fo_all.data(),None),
          fom_work     = flex.mean_default(fom.select(sel_w),None),
          pher_work    = flex.mean_default(pher_w.select(sel_w),None),
          pher_free    = flex.mean_default(pher_t.select(sel_t),None),
          sigmaa       = sigmaa_bin)
        result.append(bin)
    return result

  def show_rwork_rfree_number_completeness(self, prefix="", title=None, out = None):
    if(out is None): out = sys.stdout
    if(title is not None):
      print >> out, prefix+title
    print >> out,\
      prefix+" BIN  RESOLUTION RANGE  COMPL.    NWORK NFREE   RWORK  RFREE"
    fmt = " %s %s    %s %s %s  %s %s"
    for bin in self.bins:
      print >> out,prefix+fmt%(
        format_value("%3d", bin.i_bin),
        format_value("%-17s", bin.d_range),
        format_value("%4.2f", bin.completeness),
        format_value("%8d", bin.n_work),
        format_value("%5d", bin.n_free),
        format_value("%6.4f", bin.r_work),
        format_value("%6.4f", bin.r_free))

  def show_remark_3(self, out = None):
    from cctbx import sgtbx
    if(out is None): out = sys.stdout
    pr = "REMARK   3  "
    print >> out,pr+"REFINEMENT TARGET : %s"%self.target_name.upper()
    print >> out,pr
    print >> out,pr+"DATA USED IN REFINEMENT."
    print >> out,pr+" RESOLUTION RANGE HIGH (ANGSTROMS) : %s"%format_value("%-8.3f", self.d_min)
    print >> out,pr+" RESOLUTION RANGE LOW  (ANGSTROMS) : %s"%format_value("%-8.3f", self.d_max)
    print >> out,pr+" MIN(FOBS/SIGMA_FOBS)              : %s"%format_value("%-6.2f", self.min_f_obs_over_sigma)
    print >> out,pr+" COMPLETENESS FOR RANGE        (%s) : %-6.2f"%\
      ("%", self.completeness_in_range*100.0)
    print >> out,pr+" NUMBER OF REFLECTIONS             : %-10d"%self.number_of_reflections
    print >> out,pr
    print >> out,pr+"FIT TO DATA USED IN REFINEMENT."
    print >> out,pr+" R VALUE     (WORKING + TEST SET) : %s"%format_value("%-6.4f",self.r_all)
    print >> out,pr+" R VALUE            (WORKING SET) : %s"%format_value("%-6.4f", self.r_work)
    print >> out,pr+" FREE R VALUE                     : %s"%format_value("%-6.4f", self.r_free)
    print >> out,pr+" FREE R VALUE TEST SET SIZE   (%s) : %-6.2f"%("%",
      float(self.number_of_test_reflections)/self.number_of_reflections*100.)
    print >> out,pr+" FREE R VALUE TEST SET COUNT      : %-10d"%self.number_of_test_reflections
    print >> out,pr
    self.show_rwork_rfree_number_completeness(prefix = pr,
      title = "FIT TO DATA USED IN REFINEMENT (IN BINS).", out = out)
    print >> out,pr
    print >> out,pr+"BULK SOLVENT MODELLING."
    print >> out,pr+" METHOD USED        : FLAT BULK SOLVENT MODEL"
    print >> out,pr+" SOLVENT RADIUS     : %s"%format_value("%-8.2f", self.mask_solvent_radius)
    print >> out,pr+" SHRINKAGE RADIUS   : %s"%format_value("%-8.2f", self.mask_shrink_radius)
    print >> out,pr+" GRID STEP FACTOR   : %s"%format_value("%-8.2f", self.mask_grid_step_factor)
    if( len(self.k_sol) == 1 ):
      print >> out,pr+" K_SOL              : %s"%format_value("%-8.3f", self.k_sol[0])
    else:
      # TODO: pretty print
      print >> out,pr," K_SOLS    : ", self.k_sol
    print >> out,pr+" B_SOL              : %s"%format_value("%-8.3f", self.b_sol)
    print >> out,pr
    if(self.twin_fraction is not None):
      print >> out,pr+"TWINNING INFORMATION."
      print >> out,pr+" FRACTION: %s"%format_value("%-8.3f", self.twin_fraction)
      print >> out,pr+" OPERATOR: %s"%\
        format_value("%-s", sgtbx.change_of_basis_op(self.twin_law).as_hkl())
    print >> out,pr+"ERROR ESTIMATES."
    print >> out,pr+" COORDINATE ERROR (MAXIMUM-LIKELIHOOD BASED)     : %s"%\
      format_value("%-8.2f", self.ml_coordinate_error)
    print >> out,pr+" PHASE ERROR (DEGREES, MAXIMUM-LIKELIHOOD BASED) : %s"%\
      format_value("%-8.2f", self.ml_phase_error)
    print >> out,pr
    print >> out,pr+"OVERALL SCALE FACTORS."
    print >> out,pr+" SCALE = SUM(|F_OBS|*|F_MODEL|)/SUM(|F_MODEL|**2) : %s"%\
      format_value("%-12.4f", self.overall_scale_k1)
    print >> out,pr+" ANISOTROPIC SCALE MATRIX ELEMENTS (IN CARTESIAN BASIS)."
    print >> out,pr+"  B11 : %s"%format_value("%-15.4f", self.b_cart[0])
    print >> out,pr+"  B22 : %s"%format_value("%-15.4f", self.b_cart[1])
    print >> out,pr+"  B33 : %s"%format_value("%-15.4f", self.b_cart[2])
    print >> out,pr+"  B12 : %s"%format_value("%-15.4f", self.b_cart[3])
    print >> out,pr+"  B13 : %s"%format_value("%-15.4f", self.b_cart[4])
    print >> out,pr+"  B23 : %s"%format_value("%-15.4f", self.b_cart[5])
    print >> out,pr
    print >> out,pr+"R FACTOR FORMULA."
    print >> out,pr+" R = SUM(||F_OBS|-SCALE*|F_MODEL||)/SUM(|F_OBS|)"
    print >> out,pr
    print >> out,pr+"TOTAL MODEL STRUCTURE FACTOR (F_MODEL)."
    print >> out,pr+" F_MODEL = FB_CART * (F_CALC_ATOMS + F_BULK)"
    print >> out,pr+"  F_BULK = K_SOL * EXP(-B_SOL * S**2 / 4) * F_MASK"
    print >> out,pr+"  F_CALC_ATOMS = ATOMIC MODEL STRUCTURE FACTORS"
    print >> out,pr+"  FB_CART = EXP(-H(t) * A(-1) * B * A(-1t) * H)"
    print >> out,pr+"   A = orthogonalization matrix, H = MILLER INDEX"
    print >> out,pr+"   (t) = TRANSPOSE, (-1) = INVERSE"
    print >> out,pr
    print >> out,pr+"STRUCTURE FACTORS CALCULATION ALGORITHM : %-s"%\
      self.sf_algorithm.upper()
    out.flush()

  def show_targets(self, out = None, text = ""):
    if(out is None): out = sys.stdout
    part1 = "|-"+text
    part2 = "-|"
    n = 79 - len(part1+part2)
    print >> out, part1 + "-"*n + part2
    part3 = "| target_work(%s) = %s  r_work = %s  r_free = %s" % (
      self.target_name,
      format_value(format="%.6g",  value = self.target_work),
      format_value(format="%6.4f", value = self.r_work),
      format_value(format="%6.4f", value = self.r_free))
    n = 78 - len(str(part3)+"|")
    print >> out, part3, " "*n +"|"
    print >> out, "|" +"-"*77+"|"
    out.flush()

  def _header_resolutions_nreflections(self, header, out):
    if(header is None): header = ""
    line1 = "(resolution: "
    line2 = format_value("%6.2f",self.d_min).strip()
    line3 = format_value("%6.2f",self.d_max).strip()
    line4 = " - "
    line5 = " A; n_refl. = "
    line6 = format_value("%d",self.number_of_reflections).strip()
    tl = header+"-"+line1+line2+line4+line3+line5+line6+")"
    line_len = len("|-"+"|"+tl)
    fill_len = 80-line_len-1
    print >> out, "|-"+tl+"-"*(fill_len)+"|"
    out.flush()

  def _rfactors_and_bulk_solvent_and_scale_params(self, out):
    out.flush()
    r_work = format_value("%6.4f",self.r_work).strip()
    r_free = format_value("%6.4f",self.r_free).strip()
    scale  = format_value("%6.3f",self.overall_scale_k1).strip()
    if( len(self.k_sol)==1 ):
      k_sol  = format_value("%4.2f",self.k_sol[0]).strip()
    else:
      # TODO: pretty print
      k_sol = str(self.k_sol)
    b_sol  = format_value("%6.2f",self.b_sol).strip()
    b0,b1,b2,b3,b4,b5 = n_as_s("%7.2f",self.b_cart)
    b_iso  = format_value("%7.2f",self.b_iso).strip()
    line = "| r_work= "+r_work+"   r_free= "+r_free+"   ksol= "+k_sol+\
           "   Bsol= "+b_sol+"   scale= "+scale
    np = 79 - (len(line) + 1)
    if(np < 0): np = 0
    print >> out, line + " "*np + "|"
    print >> out, "| "+"  "*38+"|"
    print >> out, "| overall anisotropic scale matrix (Cartesian basis; B11,B22,B33,B12,B13,B23):|"
    c = ","
    line4 = "| ("+b0+c+b1+c+b2+c+b3+c+b4+c+b5+"); trace/3= "+b_iso
    np = 79 - (len(line4) + 1)
    line4 = line4 + " "*np + "|"
    print >> out, line4
    out.flush()

  def show_rfactors_targets_scales_overall(self, header = None, out=None):
    from cctbx import sgtbx
    if(out is None): out = sys.stdout
    out.flush()
    p = " "
    self._header_resolutions_nreflections(header=header, out=out)
    print >> out, "| "+"  "*38+"|"
    self._rfactors_and_bulk_solvent_and_scale_params(out=out)
    err = format_value("%6.2f",self.ml_coordinate_error)
    print >> out, "| "+"  "*38+"|"
    line6="| maximum likelihood estimate for coordinate error: "+err+" A"
    np = 79 - (len(line6) + 1)
    line6 = line6 + " "*np + "|"
    print >> out, line6
    line7="| x-ray target function (%s) for work reflections: %s"% (
      self.target_name, n_as_s("%15.6f",self.target_work))
    np = 79 - (len(line7) + 1)
    line7 = line7 + " "*np + "|"
    print >> out, line7
    if(self.twin_fraction is not None):
      line8="| twin fraction: "+format_value("%-4.2f",self.twin_fraction)+\
        "  twin operator: "+\
        format_value("%-s",sgtbx.change_of_basis_op(self.twin_law).as_hkl())
      np = 79 - (len(line8) + 1)
      line8 = line8 + " "*np + "|"
      print >> out, line8
    print >> out, "|"+"-"*77+"|"
    out.flush()

  def show_rfactors_targets_in_bins(self, out = None):
    if(out is None): out = sys.stdout
    print >> out, "|"+"-"*77+"|"
    print >> out, "| Bin     Resolution   Compl.  No. Refl.    R-factors          Targets        |"
    print >> out, "|number     range              work test   work   test        work        test|"
    for bin in self.bins:
      print >> out, "|%3d: %-17s %4.2f %6d %4d %s %s %s %s|" % (
        bin.i_bin,
        bin.d_range,
        bin.completeness,
        bin.n_work,
        bin.n_free,
        format_value("%6.4f",  bin.r_work),
        format_value("%6.4f",  bin.r_free),
        format_value("%11.5g", bin.target_work),
        format_value("%11.5g", bin.target_free))
    print >> out, "|"+"-"*77+"|"
    out.flush()

  def show_fom_pher_alpha_beta_in_bins(self, out = None):
    if(out is None): out = sys.stdout
    print >> out, "|"+"-"*77+"|"
    print >> out, "|R-free likelihood based estimates for figures of merit, absolute phase error,|"
    print >> out, "|and distribution parameters alpha and beta (Acta Cryst. (1995). A51, 880-887)|"
    print >> out, "|"+" "*77+"|"
    print >> out, "| Bin     Resolution      No. Refl.   FOM  Phase Scale    Alpha        Beta   |"
    print >> out, "|  #        range        work  test        error factor                       |"
    for bin in self.bins:
      print >> out, "|%3d: %-17s%6d%6d%s%s%s%s%s|" % (
        bin.i_bin,
        bin.d_range,
        bin.n_work,
        bin.n_free,
        format_value("%6.2f",  bin.fom_work),
        format_value("%7.2f",  bin.pher_work),
        format_value("%7.2f",  bin.scale_k1_work),
        format_value("%9.2f",  bin.alpha_work),
        format_value("%14.2f", bin.beta_work))
    print >>out, "|alpha:            min =%s max =%s mean =%s|"%(
      format_value("%12.2f", self.alpha_work_min),
      format_value("%16.2f", self.alpha_work_max),
      format_value("%13.2f", self.alpha_work_mean))
    print >>out, "|beta:             min =%s max =%s mean =%s|"%(
      format_value("%12.2f", self.beta_work_min),
      format_value("%16.2f", self.beta_work_max),
      format_value("%13.2f", self.beta_work_mean))
    print >>out, "|figures of merit: min =%s max =%s mean =%s|"%(
      format_value("%12.2f", self.fom_work_min),
      format_value("%16.2f", self.fom_work_max),
      format_value("%13.2f", self.fom_work_mean))
    print >>out, "|phase err.(work): min =%s max =%s mean =%s|"%(
      format_value("%12.2f", self.pher_work_min),
      format_value("%16.2f", self.pher_work_max),
      format_value("%13.2f", self.pher_work_mean))
    print >>out, "|phase err.(test): min =%s max =%s mean =%s|"%(
      format_value("%12.2f", self.pher_free_min),
      format_value("%16.2f", self.pher_free_max),
      format_value("%13.2f", self.pher_free_mean))
    print >> out, "|"+"-"*77+"|"
    out.flush()

  def show_all(self, header = "", out = None):
    if(out is None): out = sys.stdout
    self.show_rfactors_targets_scales_overall(header = header, out = out)
    print >> out
    self.show_rfactors_targets_in_bins(out = out)
    print >> out
    self.show_fom_pher_alpha_beta_in_bins(out = out)

  # re-arrange binned statistics for phenix GUI (or logfile)
  def export_bins_table_data (self, title="Statistics by resolution bin") :
    return export_bins_table_data(self.bins, title)

def show_histogram(data, n_slots, log):
  from cctbx.array_family import flex
  hm = flex.histogram(data = data, n_slots = n_slots)
  lc_1 = hm.data_min()
  s_1 = enumerate(hm.slots())
  for (i_1,n_1) in s_1:
    hc_1 = hm.data_min() + hm.slot_width() * (i_1+1)
    print >> log, "%10.3f - %-10.3f : %d" % (lc_1, hc_1, n_1)
    lc_1 = hc_1

def export_bins_table_data (bins, title="Statistics by resolution bin") :
  from iotbx import data_plots
  table_stats = ["r_work", "r_free", "completeness", "fom_work",
                 "pher_free", "scale_k1_work"]
  labels = ["Resolution", "R-work", "R-free", "Completeness", "FOM",
                   "Phase error", "Scale factor"]
  graph_names = ["R-work/R-free vs. resolution",
                 "Completeness vs. resolution",
                 "Figure of merit vs. resolution",
                 "Phase error vs. resolution",
                 "Scale factor vs. resolution"]
  graph_columns = [[0,1,2], [0,3], [0,4], [0,5], [0,6]]
  if hasattr(bins[0], "sigmaa") and (bins[0].sigmaa is not None) :
    table_stats.append("sigmaa")
    labels.append("SigmaA")
    graph_names.append("SigmaA vs. resolution")
    graph_columns.append([0,7])
  data_rows = []
  for bin in bins :
    bin_stats = []
    (min_res_str, max_res_str) = re.sub("\s*", "", bin.d_range).split("-")
    (min_res, max_res) = (string.atof(min_res_str), string.atof(max_res_str))
    bin_stats.append(1 / (max_res ** 2))
    for stat_attr_name in table_stats :
      bin_stats.append(getattr(bin, stat_attr_name))
    data_rows.append(bin_stats)
  data = [[row[i] for row in data_rows] for i in xrange(len(data_rows[0]))]
  t = data_plots.table_data(
    title=title,
    column_labels=labels,
    graph_names=graph_names,
    graph_columns=graph_columns,
    data=data,
    x_is_inverse_d_min=True)
  return t
