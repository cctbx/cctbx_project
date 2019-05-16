from __future__ import absolute_import, division, print_function
from mmtbx.scaling import absolute_scaling
import mmtbx.scaling
from iotbx import data_plots
from scitbx.array_family import flex
from scitbx.math import scale_curves
from scitbx import simplex
from scitbx.math import chebyshev_polynome
from libtbx.utils import Sorry
from libtbx import table_utils
import sys,math
from six.moves import zip
from six.moves import range

low_lim = 0.00142857142857
high_lim = 0.914857142857

mean_coefs = flex.double([-0.43300151026105321, 0.33427752644857256, -0.36429668743046412, 0.32015342663362861, -0.33187593112421665, 0.25808435074324798, -0.21097660027094089, 0.19245214093632027, -0.16848242519036311, 0.15291501268526811, -0.12910039774102444, 0.099417816273142834, -0.088720990827720211, 0.095701367393211653, -0.10629070210548312, 0.10981060473882812, -0.098166077509690655, 0.081269422162739482, -0.070505309537945482, 0.062534977493702376, -0.052309784988656668, 0.043125649942558686, -0.033948257000056978, 0.025057256942873911, -0.017846896900315476, 0.013015472828176888, -0.0090867510653793015, 0.005271961028977065, -0.00066906935937178775, -0.001734805197078687, 0.0019473391336692384, -0.00077699241439945844, -0.00012166448304738755, 0.0022859384260305684, -0.0041703347977807507, 0.0048616662471392688, -0.0049936320989681648, 0.0062946025455139906, -0.0054911210539810868, 0.0040222680071224075, -0.0033087237096459769, 0.0042634379859494151, -0.0037060560347156168, 0.0026770515762505058, -0.0020954947717182685, 0.0035512064084911323, -0.0028854530832642875, 0.002100343979509825, -0.0014536705634179688, 0.0024292695349115174])

std_coefs = flex.double([0.22609807236783072, -0.051083004382385722, 0.10050223345603099, -0.059797469000968342, 0.078452075293119358, -0.061376061912225756, 0.046001019180730129, -0.046818252753277688, 0.037535878728343949, -0.031883497361025873, 0.031132854775465228, -0.026248228833806508, 0.025229855893282804, -0.022987539515026915, 0.018603638709078982, -0.020688685663116515, 0.021490882895477355, -0.019155463126466928, 0.018694555361791723, -0.017919220545523508, 0.01688432732243788, -0.016177982330096936, 0.013523772618558827, -0.011497460798395623, 0.010090183879313928, -0.0077311745570573494, 0.0069765868372828593, -0.0085904927919333608, 0.0079389398144072369, -0.0063106616713442193, 0.0072030470015979342, -0.0082688707324504833, 0.0075456582719407002, -0.0078122483377159966, 0.007131121698384397, -0.004898714984268495, 0.0045473543279292298, -0.0055478491205527948, 0.0041818036356804219, -0.0032442174724577502, 0.0035282617908206485, -0.0026738719276938735, 0.0012942832126333331, -0.001864418991069073, 0.001979588643470585, -0.001413729012970848, 0.00074827896319899767, -0.00089235624086152622, 0.00061639083311362331, -0.0007922443411235876])

class relative_wilson(mmtbx.scaling.xtriage_analysis):
  def __init__(self,
      miller_obs,
      miller_calc,
      min_d_star_sq=0.0,
      max_d_star_sq=2.0,
      n_points=2000,
      level=6.0):
    assert miller_obs.indices().all_eq(miller_calc.indices())
    if (miller_obs.is_xray_amplitude_array()):
      miller_obs = miller_obs.f_as_f_sq()
    if (miller_calc.is_xray_amplitude_array()):
      miller_calc = miller_calc.f_as_f_sq()
    self.obs  = miller_obs.deep_copy()
    self.calc = miller_calc.deep_copy()
    self.mind = min_d_star_sq
    self.maxd = max_d_star_sq
    self.m    = n_points
    self.n    = 2
    self.level = level

    norma_obs  = absolute_scaling.kernel_normalisation(
      miller_array=self.obs,
      auto_kernel=True,
      n_bins=45,
      n_term=17)
    norma_calc = absolute_scaling.kernel_normalisation(
      miller_array=self.calc,
      auto_kernel=True,
      n_bins=45,
      n_term=17)

    obs_d_star_sq  = norma_obs.d_star_sq_array
    calc_d_star_sq = norma_calc.d_star_sq_array
    sel_calc_obs = norma_calc.bin_selection.select(norma_obs.bin_selection)
    sel_obs_calc = norma_obs.bin_selection.select(norma_calc.bin_selection)
    sel  = ((obs_d_star_sq > low_lim) & (obs_d_star_sq < high_lim) &
            (norma_obs.mean_I_array > 0))
    sel = sel.select(sel_calc_obs)

    self.obs_d_star_sq  = obs_d_star_sq.select( sel )
    self.calc_d_star_sq = calc_d_star_sq.select( sel_obs_calc ).select(sel)
    self.mean_obs       = norma_obs.mean_I_array.select(sel)
    self.mean_calc      = norma_calc.mean_I_array.select(
                            sel_obs_calc).select(sel)
    self.var_obs        = norma_obs.var_I_array.select(sel)
    self.var_calc       = norma_calc.var_I_array.select(
      sel_obs_calc).select(sel)

    # make an interpolator object please
    self.interpol = scale_curves.curve_interpolator( self.mind, self.maxd,
      self.m)
    # do the interpolation
    tmp_obs_d_star_sq  , self.mean_obs,self.obs_a  , self.obs_b  = \
      self.interpol.interpolate(self.obs_d_star_sq,self.mean_obs)
    self.obs_d_star_sq , self.var_obs,self.obs_a   , self.obs_b  = \
      self.interpol.interpolate(self.obs_d_star_sq, self.var_obs)
    tmp_calc_d_star_sq , self.mean_calc,self.calc_a, self.calc_b = \
      self.interpol.interpolate(self.calc_d_star_sq,self.mean_calc)
    self.calc_d_star_sq, self.var_calc,self.calc_a , self.calc_b = \
      self.interpol.interpolate(self.calc_d_star_sq,self.var_calc)

    self.mean_ratio_engine = chebyshev_polynome( mean_coefs.size(),
      low_lim-1e-3, high_lim+1e-3,mean_coefs)
    self.std_ratio_engine = chebyshev_polynome( std_coefs.size(),
      low_lim-1e-3, high_lim+1e-3,std_coefs)

    self.x = flex.double([0,0])

    self.low_lim_for_scaling = 1.0/(4.0*4.0) #0.0625
    selection = (self.calc_d_star_sq > self.low_lim_for_scaling)
    if (selection.count(True) == 0):
      raise Sorry("No reflections within required resolution range after "+
        "filtering.")
    self.weight_array = selection.as_double() / (2.0 * self.var_obs)
    assert (not self.weight_array.all_eq(0.0))

    self.mean   = flex.double( [1.0/(flex.sum(self.mean_calc) /
                                flex.sum(self.mean_obs)), 0.0 ] )
    self.sigmas = flex.double( [0.5, 0.5] )

    s = 1.0/(flex.sum(self.weight_array*self.mean_calc)/
             flex.sum(self.weight_array*self.mean_obs))
    b = 0.0
    self.sart_simplex = [ flex.double([s,b]), flex.double([s+0.1,b+1.1]),
                          flex.double([s-0.1,b-1.1]) ]
    self.opti = simplex.simplex_opt( 2, self.sart_simplex, self)

    sol = self.opti.get_solution()
    self.scale   = abs(sol[0])
    self.b_value = sol[1]

    self.modify_weights()
    self.all_bad_z_scores = self.weight_array.all_eq(0.0)
    if (not self.all_bad_z_scores):
      s = 1.0/(flex.sum(self.weight_array*self.mean_calc) /
               flex.sum(self.weight_array*self.mean_obs))
      b = 0.0
      self.sart_simplex = [ flex.double([s,b]), flex.double([s+0.1,b+1.1]),
                            flex.double([s-0.1,b-1.1]) ]
      self.opti = simplex.simplex_opt( 2, self.sart_simplex, self)
    #self.mean_calc = self.mean_calc*self.scale*flex.exp(self.calc_d_star_sq*self.b_value)

  def summary(self):
    i_scaled = flex.exp( self.calc_d_star_sq*self.b_value ) * \
                self.mean_calc * self.scale
    sel = (self.mean_obs > 0).iselection()
    ratio  = flex.log(i_scaled.select(sel) / self.mean_obs.select(sel))
    ratio_ = flex.double(self.mean_obs.size(), 0)
    ratio_.set_selected(sel, ratio)
    curves = [
      self.calc_d_star_sq,
      -ratio_, # observed
      self.curve( self.calc_d_star_sq ), # expected
      self.get_z_scores(self.scale, self.b_value)
    ]
    return summary(
      all_curves=curves,
      level=self.level,
      all_bad_z_scores=self.all_bad_z_scores)

  def modify_weights(self,level=5):
    z_scores = self.get_z_scores(self.scale, self.b_value)
    sel  = flex.double(list(flex.bool(z_scores<level)))
    self.weight_array = self.weight_array*sel

  def get_z_scores(self, scale, b_value):
    i_scaled = flex.exp( self.calc_d_star_sq*b_value )*self.mean_calc*scale
    sel = ((self.mean_obs > 0) & (i_scaled > 0)) .iselection()
    ratio  = self.mean_obs.select(sel) / i_scaled.select(sel)
    mean = self.curve( self.calc_d_star_sq ).select(sel)
    assert ratio.all_gt(0) # FIXME need to filter first!
    ratio = flex.log(ratio)
    var = self.std(self.calc_d_star_sq).select(sel)
    d_star_sq = self.calc_d_star_sq.select(sel)
    assert var.all_ne(0)
    z = flex.abs(ratio-mean)/var
    z_ = flex.double(self.mean_obs.size(), -1)
    z_.set_selected(sel, z)
    return z_

  def target(self,vector):
    v=1.0
    scale = abs(vector[0])
    b_value = vector[1]
    if b_value > 200.0:
      b_value = 200.0
    if b_value < -200.0:
      b_value = -200.0
    i_scaled = flex.exp( self.calc_d_star_sq*b_value )*self.mean_calc*scale
    ratio  = i_scaled / self.mean_obs
    curve = self.curve( self.calc_d_star_sq )
    result = ratio - flex.exp(curve)
    if (flex.max(result) > math.sqrt(sys.float_info.max)):
      raise OverflowError("Result array exceeds floating-point limit.")
    result = result*result
    wmax = flex.max(self.weight_array)
    assert (wmax != 0)
    if (wmax > 1) and (flex.max(result) > sys.float_info.max / wmax):
      raise OverflowError("Weighted result array will exceed floating-point "+
        "limit: %e" % flex.max(result))
    result = result*self.weight_array
    result = flex.sum( result )
    return result

  def curve(self,d_star_sq):
    result =  self.mean_ratio_engine.f( d_star_sq )
    return result

  def std(self,d_star_sq):
    result = self.std_ratio_engine.f( d_star_sq )
    return result

  def show_summary(self, out):
    return self.summary().show(out=out)

class summary(mmtbx.scaling.xtriage_analysis):
  def __init__(self,
      all_curves,
      level=6.0,
      all_bad_z_scores=False):
    self.table = data_plots.table_data(
      title="Relative Wilson plot",
      column_labels=["Max. resolution", "log(I_exp/I_obs)", "Reference curve",
        "Z-score"],
      graph_names=["Relative Wilson plot"],
      graph_labels=[("High resolution", "")],
      graph_columns=[list(range(4))],
      x_is_inverse_d_min=True,
      data=[ list(array) for array in all_curves ])
    self.cutoff = level
    self.all_bad_z_scores = all_bad_z_scores

  def n_outliers(self):
    ss,rr,ii,zz = self.data_as_flex_arrays()
    flagged = zz > self.cutoff
    return flagged.count(True)

  def data_as_flex_arrays(self):
    return [ flex.double(column) for column in self.table.data ]

  def _show_impl(self, out):
    ss,rr,ii,zz = self.data_as_flex_arrays()
    flagged = zz > self.cutoff
    sel_ss = ss.select(flagged)
    sel_z = zz.select(flagged)
    sel_r = rr.select(flagged)
    sel_i = ii.select(flagged)
    out.show_sub_header("Relative Wilson plot")
    out.show_text("""\
The relative Wilson plot compares the mean intensity of the observed data with
the mean intensity computed from the model, as a function of resolution.  This
curve is expected to fall off at low resolution if no contribution for bulk
solvent is provided for the calculated intensities, because the presence of
bulk solvent reduces the observed intensities at low resolution by reducing
the contrast.  At high resolution, the curve should be a straight line with a
slope that reflects the difference in overall B-factor between the model and
the data.  Compared to the normal Wilson plot, the relative Wilson plot is
more linear because the influence of favored distances between atoms, caused
by bonding and secondary structure, is cancelled out.
""")
    out.show_plot(self.table)
    if (self.all_bad_z_scores):
      out.warn("""\
All resolution shells have Z-scores above %4.2f sigma.  This is indicative of
severe problems with the input data, including processing errors or ice rings.
We recommend checking the logs for data processing and inspecting the raw
images.\n""" % self.cutoff)
    else :
      out.show_text("""\
All relative wilson plot outliers above %4.2f sigma are reported.
""" % self.cutoff)
    out.newline()
    rows = []
    if len(sel_ss) > 0:
      for s,z,r,i in zip(sel_ss,sel_z,sel_r,sel_i):
        sss = math.sqrt(1.0/s)
        rows.append([ "%8.2f" % sss, "%9.3e" % r, "%9.3e" % i, "%5.2f" % z ])
      table = table_utils.simple_table(
        column_headers=["d-spacing", "Obs. Log[ratio]", "Expected Log[ratio]",
          "Z-score"],
        table_rows=rows)
      out.show_table(table)
    else:
      out.show("The Relative wilson plot doesn't indicate any serious errors.")
