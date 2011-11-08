from mmtbx.scaling import absolute_scaling
from iotbx import data_plots
from scitbx.array_family import flex
from scitbx.math import scale_curves
from scitbx import simplex
from scitbx.math import chebyshev_polynome
import sys,math

low_lim = 0.00142857142857
high_lim = 0.914857142857

mean_coefs = flex.double([-0.43300151026105321, 0.33427752644857256, -0.36429668743046412, 0.32015342663362861, -0.33187593112421665, 0.25808435074324798, -0.21097660027094089, 0.19245214093632027, -0.16848242519036311, 0.15291501268526811, -0.12910039774102444, 0.099417816273142834, -0.088720990827720211, 0.095701367393211653, -0.10629070210548312, 0.10981060473882812, -0.098166077509690655, 0.081269422162739482, -0.070505309537945482, 0.062534977493702376, -0.052309784988656668, 0.043125649942558686, -0.033948257000056978, 0.025057256942873911, -0.017846896900315476, 0.013015472828176888, -0.0090867510653793015, 0.005271961028977065, -0.00066906935937178775, -0.001734805197078687, 0.0019473391336692384, -0.00077699241439945844, -0.00012166448304738755, 0.0022859384260305684, -0.0041703347977807507, 0.0048616662471392688, -0.0049936320989681648, 0.0062946025455139906, -0.0054911210539810868, 0.0040222680071224075, -0.0033087237096459769, 0.0042634379859494151, -0.0037060560347156168, 0.0026770515762505058, -0.0020954947717182685, 0.0035512064084911323, -0.0028854530832642875, 0.002100343979509825, -0.0014536705634179688, 0.0024292695349115174])

std_coefs = flex.double([0.22609807236783072, -0.051083004382385722, 0.10050223345603099, -0.059797469000968342, 0.078452075293119358, -0.061376061912225756, 0.046001019180730129, -0.046818252753277688, 0.037535878728343949, -0.031883497361025873, 0.031132854775465228, -0.026248228833806508, 0.025229855893282804, -0.022987539515026915, 0.018603638709078982, -0.020688685663116515, 0.021490882895477355, -0.019155463126466928, 0.018694555361791723, -0.017919220545523508, 0.01688432732243788, -0.016177982330096936, 0.013523772618558827, -0.011497460798395623, 0.010090183879313928, -0.0077311745570573494, 0.0069765868372828593, -0.0085904927919333608, 0.0079389398144072369, -0.0063106616713442193, 0.0072030470015979342, -0.0082688707324504833, 0.0075456582719407002, -0.0078122483377159966, 0.007131121698384397, -0.004898714984268495, 0.0045473543279292298, -0.0055478491205527948, 0.0041818036356804219, -0.0032442174724577502, 0.0035282617908206485, -0.0026738719276938735, 0.0012942832126333331, -0.001864418991069073, 0.001979588643470585, -0.001413729012970848, 0.00074827896319899767, -0.00089235624086152622, 0.00061639083311362331, -0.0007922443411235876])



class relative_wilson(object):
  def __init__(self, miller_obs, miller_calc, min_d_star_sq=0.0, max_d_star_sq=2.0, n_points=2000):
    self.obs  = miller_obs.deep_copy()
    self.calc = miller_calc.deep_copy()
    self.mind = min_d_star_sq
    self.maxd = max_d_star_sq
    self.m    = n_points
    self.n    = 2

    self.calc = self.calc.f_as_f_sq()

    self.norma_obs  = absolute_scaling.kernel_normalisation( self.obs,auto_kernel=True,n_bins=45, n_term=17)
    self.norma_calc = absolute_scaling.kernel_normalisation( self.calc,auto_kernel=True,n_bins=45,n_term=17)

    self.obs_d_star_sq  = self.norma_obs.d_star_sq_array
    self.calc_d_star_sq = self.norma_calc.d_star_sq_array
    sel  = flex.bool(self.obs_d_star_sq > low_lim)& flex.bool(self.obs_d_star_sq<high_lim)

    self.obs_d_star_sq = self.obs_d_star_sq.select( sel )
    self.calc_d_star_sq = self.calc_d_star_sq.select( sel )
    self.mean_obs       = self.norma_obs.mean_I_array.select(sel)
    self.mean_calc      = self.norma_calc.mean_I_array.select(sel)
    self.var_obs        = self.norma_obs.var_I_array.select(sel)
    self.var_calc       = self.norma_calc.var_I_array.select(sel)

    # make an interpolator object please
    self.interpol = scale_curves.curve_interpolator( self.mind, self.maxd, self.m)
    # do the interpolation
    tmp_obs_d_star_sq  , self.mean_obs,self.obs_a  , self.obs_b  = self.interpol.interpolate(self.obs_d_star_sq,self.mean_obs)
    self.obs_d_star_sq , self.var_obs,self.obs_a   , self.obs_b  = self.interpol.interpolate(self.obs_d_star_sq, self.var_obs)
    tmp_calc_d_star_sq , self.mean_calc,self.calc_a, self.calc_b = self.interpol.interpolate(self.calc_d_star_sq,self.mean_calc)
    self.calc_d_star_sq, self.var_calc,self.calc_a , self.calc_b = self.interpol.interpolate(self.calc_d_star_sq,self.var_calc)

    self.mean_ratio_engine = chebyshev_polynome( mean_coefs.size(), low_lim-1e-3, high_lim+1e-3,mean_coefs)
    self.std_ratio_engine = chebyshev_polynome( std_coefs.size(), low_lim-1e-3, high_lim+1e-3,std_coefs)


    self.x = flex.double([0,0])

    self.low_lim_for_scaling = 1.0/(4.0*4.0) #0.025
    selection = (self.calc_d_star_sq > self.low_lim_for_scaling)
    if (selection.count(True) == 0) :
      raise RuntimeError("No reflections meeting selection criteria.")
    self.weight_array = selection.as_double() / (2.0 * self.var_obs)
    assert (not self.weight_array.all_eq(0.0))

    self.mean   = flex.double( [1.0/(flex.sum(self.mean_calc)/flex.sum(self.mean_obs)), 0.0 ] )
    self.sigmas = flex.double( [0.5, 0.5] )

    s = 1.0/(flex.sum(self.weight_array*self.mean_calc)/flex.sum(self.weight_array*self.mean_obs))
    b = 0.0
    self.sart_simplex = [ flex.double([s,b]), flex.double([s+0.1,b+1.1]), flex.double([s-0.1,b-1.1]) ]
    self.opti = simplex.simplex_opt( 2, self.sart_simplex, self)

    sol = self.opti.get_solution()
    self.scale   = abs(sol[0])
    self.b_value = sol[1]

    self.modify_weights()
    s = 1.0/(flex.sum(self.weight_array*self.mean_calc)/flex.sum(self.weight_array*self.mean_obs))
    b = 0.0
    self.sart_simplex = [ flex.double([s,b]), flex.double([s+0.1,b+1.1]), flex.double([s-0.1,b-1.1]) ]
    self.opti = simplex.simplex_opt( 2, self.sart_simplex, self)
    #self.mean_calc = self.mean_calc*self.scale*flex.exp(self.calc_d_star_sq*self.b_value)

    self.show_summary()


  def modify_weights(self,level=5):
    z_scores = self.get_z_scores(self.scale, self.b_value)
    sel  = flex.double(list(flex.bool(z_scores<level)))
    self.weight_array = self.weight_array*sel


  def get_z_scores(self, scale, b_value):
    i_scaled = flex.exp( self.calc_d_star_sq*b_value )*self.mean_calc*scale
    ratio  = i_scaled / self.mean_obs
    mean = self.curve( self.calc_d_star_sq )
    ratio = flex.log(ratio)
    var = self.std(self.calc_d_star_sq)
    z = flex.abs(ratio-mean)/var
    return z

  def get_data_plot (self) :
    table = data_plots.table_data(
      title="Relative Wilson plot",
      column_labels=["Max. resolution", "log(I_exp/I_obs)", "Z-score",
        "Reference curve"],
      graph_names=["Relative Wilson plot"],
      graph_labels=[("High resolution", "")],
      graph_columns=[list(range(4))],
      x_is_inverse_d_min=True)
    ss,rr,zz,ii = self.get_all_curves()
    for s,z,r,i in zip(ss,zz,rr,ii):
      table.add_row([s,z,r,i])
    return table

  def get_all_curves(self):
    i_scaled = flex.exp( self.calc_d_star_sq*self.b_value )*self.mean_calc*self.scale
    ratio  = i_scaled / self.mean_obs
    ratio = flex.log(ratio)
    return self.calc_d_star_sq, -ratio, flex.abs(ratio)/self.std(self.calc_d_star_sq), self.curve( self.calc_d_star_sq )


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
    result = result*result
    result = result*self.weight_array
    result = flex.sum( result )
    return result


  def curve(self,d_star_sq):
    result =  self.mean_ratio_engine.f( d_star_sq )
    return result

  def std(self,d_star_sq):
    result = self.std_ratio_engine.f( d_star_sq )
    return result

  def show(self,out=None):
    if out is None:
      out = sys.stdout
    ss,rr,zz,ii = self.get_all_curves()
    for s,r,z,i in zip( ss,rr,zz,ii ):
      print >> out, s,r,z,i

  def show_summary(self,level=6.0,out=None):
    if out is None:
      out = sys.stdout

    ss,rr,zz,ii = self.get_all_curves()
    flagged = flex.bool(zz>level)
    sel_ss = ss.select(flagged)
    sel_z = zz.select(flagged)
    sel_r = rr.select(flagged)
    sel_i = ii.select(flagged)
    print >> out
    print >> out, "  Relative Wilson plot Analayses  "
    print >> out, "     - All relative wilson plot outliers above %4.2f sigma are reported"%level
    if len(sel_ss) > 0:
      for s,z,r,i in zip(sel_ss,sel_z,sel_r,sel_i):
        sss = math.sqrt(1.0/s)
        print >> out, "d-spacing: %5.2f  Z-score: %5.2f    Observed and expected Log[ratio]: %5.2e  %5.2e"%(sss,z,r,i)

    else:
      print >> out, "The Relative wilson plot doesn't indicate any serious errors"
    print >> out
    #self.show()
