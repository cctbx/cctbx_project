from cctbx.array_family import flex
from mmtbx import scaling
from mmtbx.scaling import absolute_scaling
from cctbx import uctbx
from cctbx import adptbx
from cctbx import sgtbx
from cctbx.sgtbx import pointgroup_tools
from cctbx import maptbx
from cctbx import crystal
import cctbx.sgtbx.lattice_symmetry
import cctbx.sgtbx.cosets
import scitbx.math
from scitbx import matrix
from scitbx.math import chebyshev_lsq
from scitbx.math import chebyshev_polynome
from scitbx.math import chebyshev_lsq_fit
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry
from cStringIO import StringIO
import math
import sys
from iotbx import data_plots
from libtbx import table_utils


class obliquity(object):
  def __init__(self, reduced_cell, rot_mx, deg=True):
    orth = matrix.sqr(reduced_cell.orthogonalization_matrix())
    frac = matrix.sqr(reduced_cell.fractionalization_matrix())
    r_info = rot_mx.info()
    self.type = r_info.type()
    self.u = rot_mx.info().ev()
    self.h = rot_mx.transpose().inverse().info().ev()
    self.t = orth * matrix.col(self.u)
    self.tau = matrix.row(self.h) * frac
    self.delta = self.t.accute_angle(self.tau, deg=deg)

class twin_law_quality(object):
  def __init__(self,
               xs,
               twin_law):
    self.xs = xs

    self.xs_niggli = xs.change_basis(
      xs.change_of_basis_op_to_niggli_cell())

    self.twin_in_nig = xs.change_of_basis_op_to_niggli_cell().apply(
      twin_law)
    self.cb_op = sgtbx.change_of_basis_op(
      str_xyz = self.twin_in_nig.r().as_xyz() )

    self.niggli_cell = self.xs_niggli.unit_cell()
    self.new_niggli_cell = self.cb_op.apply( self.niggli_cell )

  def delta_santoro(self):
    # santoros measure for quality of a twin law
    # relative (default) or absolute is possible
    old = self.niggli_cell.metrical_matrix()
    new = self.new_niggli_cell.metrical_matrix()
    old = flex.double( old )
    new = flex.double( new )
    delta = flex.abs(old - new)
    result = 0
    top = 0
    bottom = 0
    for oi, ni in zip(old,new):
      tmp = math.fabs(oi-ni)
      top += tmp
      bottom += math.fabs(oi)
    assert (bottom>0)
    delta = 100.0*top/bottom
    return( delta )

  def delta_le_page(self):
    rot_mx_current = self.cb_op.c().r().new_denominator(1)
    obl = obliquity(self.niggli_cell,
                    rot_mx_current )
    type_string = str(obl.type)+"-fold"
    return obl.delta, type_string

  def delta_lebedev(self):
    return self.twin_in_nig.r().lebedev_2005_perturbation(
      reduced_cell=self.niggli_cell)


  def strain_tensor(self):
    # this gives a tensor describing the deformation of the unit cell needed
    # to obtain a perfect match. the sum of diagonal elements describes
    # the change in volume, off diagonal components measure associated shear.
    #
    x1 = matrix.col( (1,0,0) )
    x2 = matrix.col( (0,1,0) )
    x3 = matrix.col( (0,0,1) )

    f1 = matrix.sqr(self.niggli_cell.fractionalization_matrix() )*x1
    f2 = matrix.sqr(self.niggli_cell.fractionalization_matrix() )*x2
    f3 = matrix.sqr(self.niggli_cell.fractionalization_matrix() )*x3

    xn1 =  matrix.sqr(self.new_niggli_cell.orthogonalization_matrix() )*f1
    xn2 =  matrix.sqr(self.new_niggli_cell.orthogonalization_matrix() )*f2
    xn3 =  matrix.sqr(self.new_niggli_cell.orthogonalization_matrix() )*f3

    u1 = math.sqrt( xn1.dot(xn1) )-1.0
    u2 = math.sqrt( xn2.dot(xn2) )-1.0
    u3 = math.sqrt( xn3.dot(xn3) )-1.0
    e11 = u1/1.0
    e22 = u2/1.0
    e33 = u3/1.0
    e12 = (u1 + u2)/2.0
    e13 = (u1 + u3)/2.0
    e23 = (u2 + u3)/2.0
    result=matrix.sqr( [e11,e12,e13,
                        e12,e22,e23,
                        e13,e23,e33 ] )
    return result


## python routines copied from iotbx.iotbx.reflection_statistics.
## Should be moved but are (for now) in a conveniant place.
class twin_law(object):
  def __init__(self,
               op,
               pseudo_merohedral_flag,
               axis_type,
               delta_santoro,
               delta_le_page,
               delta_lebedev):
    self.operator =  op
    self.twin_type = pseudo_merohedral_flag
    self.delta_santoro = delta_santoro
    self.delta_le_page = delta_le_page
    self.delta_lebedev = delta_lebedev
    self.axis_type = axis_type

class twin_laws(object):
  def __init__(self,
               miller_array,
               lattice_symmetry_max_delta=3.0,
               out=None):

    self.input = miller_array.eliminate_sys_absent(integral_only=True,
                                                   log=out)
    self.change_of_basis_op_to_niggli_cell \
      = self.input.change_of_basis_op_to_niggli_cell()

    self.minimum_cell_symmetry = crystal.symmetry.change_basis(
      self.input,
      cb_op=self.change_of_basis_op_to_niggli_cell)

    self.lattice_group = sgtbx.lattice_symmetry.group(
      self.minimum_cell_symmetry.unit_cell(),
      max_delta=lattice_symmetry_max_delta)

    self.intensity_symmetry = \
      self.minimum_cell_symmetry.reflection_intensity_symmetry(
        anomalous_flag=self.input.anomalous_flag())

    self.euclid = self.intensity_symmetry.space_group_info().type()\
      .expand_addl_generators_of_euclidean_normalizer(flag_k2l=True,
                                                      flag_l2n=True )

    self.operators = []
    self.m=0
    self.pm=0

    cb_op = self.change_of_basis_op_to_niggli_cell.inverse()
    for partition in sgtbx.cosets.left_decomposition(
      g=self.lattice_group,
      h=self.intensity_symmetry.space_group()
          .build_derived_acentric_group()
          .make_tidy()).partitions[1:]:
      if (partition[0].r().determinant() > 0):
        is_pseudo_merohedral=False
        twin_type =str("  M")
        self.m+=1
        euclid_check = sgtbx.space_group( self.euclid )
        try:
          euclid_check.expand_smx( partition[0] )
        except KeyboardInterupt: raise
        except:
          is_pseudo_merohedral=True
          twin_type = str(" PM")
          self.pm+=1
          self.m-=1

        if ( euclid_check.order_z() != self.euclid.order_z() ):
          is_pseudo_merohedral=True
          if is_pseudo_merohedral:
            twin_type = str(" PM")
            self.pm+=1
            self.m-=1

        tlq = twin_law_quality( miller_array,
                                cb_op.apply(partition[0]) )

        tl = twin_law( cb_op.apply(partition[0]),
                       str(twin_type),
                       tlq.delta_le_page()[1],
                       tlq.delta_santoro(),
                       tlq.delta_le_page()[0],
                       tlq.delta_lebedev()
                      )

        self.operators.append(tl)




  def show(self, out=None):
    if out is None:
      out=sys.stdout

    comments="""\
M:  Merohedral twin law
PM: Pseudomerohedral twin law"""

    if len(self.operators)!=0 :
      print >> out
      print >> out, "The following twin laws have been found:"
      print >> out
      table_labels = ('Type', 'Axis', 'R metric (%)', 'delta (le Page)', 'delta (Lebedev)', 'Twin law')
      table_rows = []
      for twin_law in self.operators:
        table_rows.append(
          [twin_law.twin_type,
           twin_law.axis_type,
           str("%5.3f"%(twin_law.delta_santoro)),
           str("%5.3f"%(twin_law.delta_le_page)),
           str("%5.3f"%(twin_law.delta_lebedev)),
           str(twin_law.operator.r().as_hkl())] )

      print >> out, table_utils.format([table_labels]+table_rows,
                                       comments=comments,
                                       has_header=True,
                                       separate_rows=False,
                                       prefix='| ',
                                       postfix=' |')
      print >> out
      print >> out, "%3.0f merohedral twin operators found"%(self.m)
      print >> out, "%3.0f pseudo-merohedral twin operators found"%(self.pm)
      print >> out, "In total, %3.0f twin operator were found"%(len(self.operators))
      print >> out
      print >> out
      assert (self.m + self.pm)==len(self.operators)
    else:
      print >> out
      print >> out, "%3.0f merohedral twin operators found"%(self.m)
      print >> out, "%3.0f pseudo-merohedral twin operators found"%(self.pm)
      print >> out, "In total, %3.0f twin operator were found"%(len(self.operators))
      print >> out
      print >> out
      assert (self.m + self.pm)==len(self.operators)


class wilson_normalised_intensities(object):
  """ making centric and acentric cut """
  def __init__(self,
               miller_array,
               normalise=True,
               out=None,
               verbose=0):

    if out is None:
      out = sys.stdout

    assert not miller_array.space_group().is_centric()
    if not miller_array.is_xray_intensity_array():
      miller_array = miller_array.f_as_f_sq()

    work_array =  miller_array.deep_copy()
    if normalise:
      normalizer = absolute_scaling.kernel_normalisation(
        miller_array, auto_kernel=True)
      work_array = normalizer.normalised_miller.deep_copy()
      work_array = work_array.select(work_array.data()>0)
    else:
      work_array = miller_array.deep_copy().set_observation_type(miller_array)
      work_array = work_array.select(work_array.data()>0)

    self.acentric = work_array.select_acentric().as_intensity_array()
    self.centric = work_array.select_centric().as_intensity_array()

    if (self.acentric.indices().size()<=0):
      raise Sorry("No acentric reflections available. Check your input")

    if verbose > -10:
      print >> out, "Number of centrics  :", self.centric.data().size()
      print >> out, "Number of acentrics :", self.acentric.data().size()



class detect_pseudo_translations(object):
  def __init__(self,
               miller_array,
               low_limit=10.0,
               high_limit=5.0,
               max_sites=100,
               height_cut=0.0,
               distance_cut=15.0,
               p_value_cut=0.05,
               out=None,verbose=0):
    if out is None:
      out=sys.stdout

    if miller_array.is_xray_intensity_array():
      miller_array = miller_array.f_sq_as_f()
    work_array = miller_array.resolution_filter(low_limit,high_limit)
    work_array = work_array.select(work_array.data()>0).set_observation_type(
      miller_array)
    if work_array.indices().size()<20:
      print >> out
      print >> out," WARNING: "
      print >> out,"  There are only %2.0f reflections between %3.1f and %3.1f A."%(
        work_array.indices().size(), low_limit, high_limit)
      print >> out,"  This might not be enough to obtain a good estimate"
      print >> out,"  of the presence or absense of pseudo translational"
      print >> out,"  symmetry."
    if work_array.indices().size()==0:
      raise Sorry("No low resolution reflections")


    if work_array.anomalous_flag():
      work_array = work_array.average_bijvoet_mates().set_observation_type(
        miller_array)

    everything_okai = True

    if (work_array.indices().size()<0):
      print >> out, \
         "The number of reflection between %3.1f and %3.1f Angstrom" \
         %( low_limit,
            high_limit )
      print >> out, "is equal to %i" %(work_array.indices().size())
      print >> out, " ##  This is not enough to obtain a reasonable estimate of"
      print >> out, " ##  the presence of translational NCS"
      everything_okai = False

    if everything_okai:


      patterson_map = work_array.patterson_map(
        symmetry_flags=maptbx.use_space_group_symmetry).apply_sigma_scaling()

      peak_list = patterson_map.tags().peak_search(
        map=patterson_map.real_map(),
        parameters=maptbx.peak_search_parameters())

      max_height = peak_list.heights()[0]

      sym_equiv_origin = sgtbx.sym_equiv_sites(
        unit_cell=patterson_map.unit_cell(),
        space_group=patterson_map.space_group(),
        original_site=(0,0,0))

      self.suspected_peaks = []

      if max_sites > peak_list.sites().size():
        max_sites = peak_list.sites().size()


      for i_peak in range(max_sites):
        height = peak_list.heights()[i_peak]/max_height*100.0
        site = peak_list.sites()[i_peak]
        dist_info = sgtbx.min_sym_equiv_distance_info(sym_equiv_origin, site)
        if (dist_info.dist() >= distance_cut):
          if (height >= height_cut):
            p_value = self.p_value(height)
            self.suspected_peaks.append( [dist_info.sym_op()*site,
                                          height, p_value,
                                          dist_info.dist()] )
      if len(self.suspected_peaks)==0:

        print >> out
        print >> out, "No patterson vectors with a length larger then"
        print >> out, "%5.2f found. removing distance constraint"%(distance_cut)
        print >> out
        distance_cut = 1e-3
        for i_peak in range(max_sites):
          height = peak_list.heights()[i_peak]/max_height*100.0
          site = peak_list.sites()[i_peak]
          dist_info = sgtbx.min_sym_equiv_distance_info(sym_equiv_origin, site)
          if (dist_info.dist() >= distance_cut):
            if (height >= height_cut):
              p_value = self.p_value(height)
              self.suspected_peaks.append( [dist_info.sym_op()*site,
                                            height, p_value,
                                            dist_info.dist()] )



      self.p_value_cut = p_value_cut
      self.mod_h = 2
      self.mod_k = 2
      self.mod_l = 2
      if everything_okai:
        self.high_peak = self.suspected_peaks[0][1]
        self.high_peak_distance = self.suspected_peaks[0][3]
        self.high_peak_xyz = self.suspected_peaks[0][0]
        self.high_p_value = self.suspected_peaks[0][2]
        if( self.high_p_value <= self.p_value_cut):
          self.guesstimate_mod_hkl()
        if verbose > 0:
          self.show(out)


  def guesstimate_mod_hkl(self):
    tmp_mod_h = 1.0/(self.high_peak_xyz[0]+1.0e-6)
    tmp_mod_k = 1.0/(self.high_peak_xyz[1]+1.0e-6)
    tmp_mod_l = 1.0/(self.high_peak_xyz[2]+1.0e-6)
    tmp_mod_h = int(math.fabs(tmp_mod_h)+0.5)
    if (tmp_mod_h>=8):
      tmp_mod_h = 2
    tmp_mod_k = int(math.fabs(tmp_mod_k)+0.5)
    if (tmp_mod_k>=8):
      tmp_mod_k = 2
    tmp_mod_l = int(math.fabs(tmp_mod_l)+0.5)
    if (tmp_mod_l>=8):
      tmp_mod_l = 2
    self.mod_h = tmp_mod_h
    self.mod_k = tmp_mod_k
    self.mod_l = tmp_mod_l

  def p_value(self, peak_height):
    x= peak_height/100.0
    result=None
    if x<1.0:
      x = x/(1.0-x)
      a = 0.06789
      b = 3.5628
      result = 1.0 - math.exp(- ((x/a)**(-b)) )
    else:
      result=0.0
    return result

  def show(self,out=None):
    if out is None:
      out = sys.stdout
    print >> out
    print >> out," Largest patterson peak with length larger then 15 Angstrom "
    print >> out
    print >> out," Frac. coord.        :%8.3f %8.3f %8.3f" %(self.high_peak_xyz)
    print >> out," Distance to origin  :%8.3f" %(self.high_peak_distance)
    print >> out," Height (origin=100) :%8.3f" %(self.high_peak)
    print >> out," p_value(height)     :%12.3e" %(self.high_p_value)
    print >> out
    print >> out,"   The reported p_value has the following meaning: "
    print >> out,"     The probability that a peak of the specified height "
    print >> out,"     or larger is found in a Patterson function of a "
    print >> out,"     macro molecule that does not have any translational"
    print >> out,"     pseudo symmetry is equal to %10.3e "%(self.high_p_value)
    print >> out,"     p_values smaller then 0.05 might indicate "
    print >> out,"     weak translation pseudo symmetry, or the self vector of "
    print >> out,"     a large anomalous scatterer such as Hg, whereas values "
    print >> out,"     smaller then 1e-3 are a very strong indication for "
    print >> out,"     the presence of translational pseudo symmetry."
    print >> out


    if self.high_p_value <= self.p_value_cut:

      print >> out
      print >> out, "The full list of patterson peaks is: "
      print >> out
      print >> out, "  x      y      z            height   p-value(height)"
      for ii in range(len(self.suspected_peaks)):
        print >> out, "(%6.3f,%6.3f,%6.3f ) :"%(
          self.suspected_peaks[ii][0]),

        print >> out,"%8.3f   (%9.3e)"%(
          self.suspected_peaks[ii][1],
          self.suspected_peaks[ii][2])

        if self.suspected_peaks[ii][2] > self.p_value_cut:
          break


class wilson_moments(object):
  def __init__(self,
               acentric_z,
               centric_z,
               out=None,
               verbose=0):
    if out is None:
      out=sys.stdout

    self.centric_i_ratio = None
    self.centric_f_ratio = None
    self.centric_e_sq_minus_one = None

    self.centric_i_ratio_library= [3.0,2.0]
    self.centric_f_ratio_library = [0.637,0.785]
    self.centric_e_sq_minus_one_library = [0.968,0.736]

    self.acentric_i_ratio = None
    self.acentric_f_ratio = None
    self.acentric_abs_e_sq_minus_one = None

    self.acentric_i_ratio_library= [2.0,1.5]
    self.acentric_f_ratio_library = [0.785,0.885]
    self.acentric_e_sq_minus_one_library = [0.736, 0.541]

    self.compute_ratios(
      acentric_z.data()/acentric_z.epsilons().data().as_double(),
      centric_z.data()/centric_z.epsilons().data().as_double())

    self.centric_present = True
    if centric_z.data().size()==0:
      self.centric_present=False

    if verbose>0:
      self.show(out)

  def compute_ratios(self, ac, c):

    if (ac.size()>0):
      mean_z = flex.mean( ac )
      mean_z_sq = flex.mean( ac*ac )
      mean_e = flex.mean( flex.sqrt(ac) )

      self.acentric_i_ratio = mean_z_sq / (mean_z*mean_z)
      self.acentric_f_ratio = mean_e*mean_e/mean_z
      self.acentric_abs_e_sq_minus_one = flex.mean( flex.abs(ac - 1.0) )

    if (c.size()>0):
      mean_z = flex.mean( c )
      mean_z_sq = flex.mean( c*c )
      mean_e = flex.mean( flex.sqrt(c) )

      self.centric_i_ratio = mean_z_sq / (mean_z*mean_z)
      self.centric_f_ratio = mean_e*mean_e/mean_z
      self.centric_abs_e_sq_minus_one = flex.mean( flex.abs(c - 1.0) )

  def show(self,out=None):
    if out is None:
      out = sys.stdout
    print >> out
    print >> out
    print >> out, "Wilson ratio and moments "
    print >> out
    print >> out, "Acentric reflections "
    print >> out, "   <I^2>/<I>^2    :%4.3f   (untwinned: %4.3f; perfect twin %4.3f)"\
          %(self.acentric_i_ratio,
            self.acentric_i_ratio_library[0],
            self.acentric_i_ratio_library[1])
    print >> out, "   <F>^2/<F^2>    :%4.3f   (untwinned: %4.3f; perfect twin %4.3f)"\
          %(self.acentric_f_ratio,
            self.acentric_f_ratio_library[0],
            self.acentric_f_ratio_library[1])
    print >> out, "   <|E^2 - 1|>    :%4.3f   (untwinned: %4.3f; perfect twin %4.3f)"\
          %(self.acentric_abs_e_sq_minus_one,
            self.acentric_e_sq_minus_one_library[0],
            self.acentric_e_sq_minus_one_library[1])
    print >> out
    print >> out
    if self.centric_present:
      print >> out, "Centric reflections "
      print >> out, "   <I^2>/<I>^2    :%4.3f   (untwinned: %4.3f; perfect twin %4.3f)"\
            %(self.centric_i_ratio,
              self.centric_i_ratio_library[0],
              self.centric_i_ratio_library[1])
      print >> out, "   <F>^2/<F^2>    :%4.3f   (untwinned: %4.3f; perfect twin %4.3f)"\
            %(self.centric_f_ratio,
              self.centric_f_ratio_library[0],
              self.centric_f_ratio_library[1])
      print >> out, "   <|E^2 - 1|>    :%4.3f   (untwinned: %4.3f; perfect twin %4.3f)"\
            %(self.centric_abs_e_sq_minus_one,
              self.centric_e_sq_minus_one_library[0],
              self.centric_e_sq_minus_one_library[1] )
      print >> out
      print >> out




class n_z_test(object):
  def __init__(self,
               normalised_acentric,
               normalised_centric,
               out=None,verbose=0):
    if out is None:
      out = sys.stdout

    centric_available = True
    acentric_available = True
    if normalised_centric.data().size() == 0:
      centric_available = False
    if normalised_acentric.data().size() == 0:
      acentric_available = False


    n_z = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    ac_theory = flex.double([0.0000, 0.0952, 0.1813, 0.2592, 0.3297, 0.3935,
                             0.4512, 0.5034, 0.5507, 0.5934, 0.6321])

    c_theory =  flex.double([0.0000, 0.2481, 0.3453, 0.4187, 0.4738, 0.5205,
                             0.5614, 0.5972, 0.6289, 0.6572, 0.6833])
    ac_obs = flex.double(11,0)
    c_obs = flex.double(11,0)


    for ii in range(10):
      ac_obs[ii+1] = (
        (normalised_acentric.data() < (ii+1.0)/10.0) ).count(True)
      if (centric_available):
        c_obs[ii+1] = (
          (normalised_centric.data() < (ii+1.0)/10.0) ).count(True)
    if acentric_available:
      ac_obs = ac_obs/float(normalised_acentric.data().size())
    if centric_available:
      c_obs = c_obs/float(normalised_centric.data().size())

    max_deviation_ac = flex.max( flex.abs(ac_obs-ac_theory)  )
    max_deviation_c = flex.max( flex.abs(c_obs-c_theory) )


    n_z_less_then_one_ac = (
      normalised_acentric.data() < 1.0  ).count(True)
    n_z_less_then_one_c =(
      normalised_centric.data() < 1.0  ).count(True)

    d_kolmogorov_smirnov_ac = max_deviation_ac/ac_theory[10]*math.sqrt(
      n_z_less_then_one_ac )
    d_kolmogorov_smirnov_c = max_deviation_c/c_theory[10]*math.sqrt(
      n_z_less_then_one_c  )

    z = flex.double(range(11))/10.0

    self.z = z
    self.ac_obs = ac_obs
    self.ac_untwinned = ac_theory
    self.frac_ac_lt_1 = 0
    if normalised_acentric.data().size()>0:
      self.frac_ac_lt_1 = n_z_less_then_one_ac/normalised_acentric.data().size()
    self.max_diff_ac = max_deviation_ac
    self.kolmogorov_smirnoff_ac = max_deviation_ac/ac_theory[10]*math.sqrt(
      n_z_less_then_one_ac )

    self.c_obs = c_obs
    self.c_untwinned = c_theory

    self.frac_c_lt_1 = 0
    if normalised_centric.data().size()>0:
      self.frac_c_lt_1 = n_z_less_then_one_c/normalised_centric.data().size()
    self.max_diff_c = max_deviation_c
    self.kolmogorov_smirnoff_c = max_deviation_c/c_theory[10]*math.sqrt(
      n_z_less_then_one_c )

    self.mean_diff_ac = flex.sum(self.ac_obs - self.ac_untwinned)/11.0
    self.mean_diff_c = flex.sum(self.c_obs - self.c_untwinned)/11.0

    if verbose > 0:
      self.show(out)



  def show(self,out=None):
    if out is None:
      out = sys.stdout
    print >> out
    print >> out,"NZ test (0<=z<1) to detect twinning and possible translational NCS"
    print >> out
    print >> out
    print >> out,"-----------------------------------------------"
    print >> out,"|  Z  | Nac_obs | Nac_theo | Nc_obs | Nc_theo |"
    print >> out,"-----------------------------------------------"
    for ii in range(11):
      print >> out,"|%4.1f | %7.3f | %8.3f | %6.3f | %7.3f |" \
            %(ii/10.0,
              self.ac_obs[ii],
              self.ac_untwinned[ii],
              self.c_obs[ii],
              self.c_untwinned[ii])

    sign_ac = '+'
    if self.mean_diff_ac < 0:
      sign_ac = '-'

    sign_c = '+'
    if self.mean_diff_c < 0:
      sign_c = '-'

    print >> out,"-----------------------------------------------"
    print >> out,"| Maximum deviation acentric      :  %4.3f    |" \
          %(self.max_diff_ac)
    print >> out,"| Maximum deviation centric       :  %4.3f    |" \
          %(self.max_diff_c)
    print >> out,"|                                             |"
    print >> out,"| <NZ(obs)-NZ(twinned)>_acentric  : %1s%4.3f    |" \
          %(sign_ac,math.fabs(self.mean_diff_ac))
    print >> out,"| <NZ(obs)-NZ(twinned)>_centric   : %1s%4.3f    |" \
          %(sign_c,math.fabs(self.mean_diff_c))
    print >> out,"-----------------------------------------------"




class britton_test(object):
  def __init__(self,
               twin_law,
               miller_array,
               cc_cut_off=0.995,
               out=None,
               verbose=0):
    if out is None:
      out = sys.stdout


    result = [0.5,1.0,0,0]
    miller_array = miller_array.select(
    miller_array.data()>0).set_observation_type(miller_array)

    if not miller_array.is_xray_intensity_array():
      miller_array = miller_array.f_as_f_sq()

    britton_plot_array = []
    britton_plot_alpha = flex.double(range(50))/101.0
    detwin_object = scaling.detwin(miller_array.indices(),
                                   miller_array.data(),
                                   miller_array.sigmas(),
                                   miller_array.space_group(),
                                   miller_array.anomalous_flag(),
                                   twin_law)

    for ii in range(50):
      alpha = (ii)/101.0
      negative_fraction = detwin_object.detwin_with_alpha(alpha)
      britton_plot_array.append(negative_fraction)

    britton_plot_array = flex.double(britton_plot_array)
    britton_plot_array = britton_plot_array - britton_plot_array[0]

    if flex.min( britton_plot_array )==flex.max( britton_plot_array ):
      not_done=False
      estimated_alpha = 0.5
    else:
      estimated_alpha = 0.0
      not_done=True

    while not_done:
      for ii in range(48):
        alpha = ii/101.0
        britton_range = (flex.double(range(ii,50)))/101.0
        britton_obs = britton_plot_array[ii:50]
        result = self.get_alpha(britton_range,britton_obs)
        if result[1]>=cc_cut_off:
          estimated_alpha = result[0]
          not_done=False
          break
      cc_cut_off-=0.005

    ## reset the cc_cut_off one step back
    cc_cut_off+=0.005

    self.alpha_cut = ii/101.0


    britton_plot_fit = flex.double(50,0)
    for ii in range(50):
      alpha = ii/101.0
      if (alpha<estimated_alpha):
        britton_plot_fit[ii]=0.0
      else:
        britton_plot_fit[ii]= result[2] + alpha*result[3]


    self.estimated_alpha = estimated_alpha
    self.correlation = result[1]
    self.britton_alpha =  britton_plot_alpha
    self.britton_obs = britton_plot_array
    self.britton_fit = britton_plot_fit

    if verbose > 0:
      self.show(out)

  def get_alpha(self, x, y):
    assert x.size() == y.size()
    mean_x = flex.mean(x)
    mean_y = flex.mean(y)
    var_x = flex.mean(x*x)-mean_x*mean_x
    var_y = flex.mean(y*y)-mean_y*mean_y
    covar_xy = flex.mean(x*y)-mean_x*mean_y

    N = float(x.size())
    m = flex.sum(x*x)- N*mean_x*mean_x
    b = covar_xy/(var_x+1.0e-6)
    a = mean_y - b*mean_x
    correlation = covar_xy/(math.sqrt(var_x*var_y)+1.0e-6)
    return [-a/(b+1.0e-6) ,  correlation, a, b]


  def show(self, out=None):
    if out is None:
      out = sys.stdout
    print >> out
    print >> out
    print >> out, "Britton analyses"
    print >> out
    print >> out,"  Extrapolation performed on  %3.2f < alpha < 0.495 "\
          %(self.alpha_cut)
    print >> out,"  Estimated twin fraction: %4.3f"%(self.estimated_alpha)
    print >> out,"  Correlation: %5.4f"%(self.correlation)


class h_test(object):
  def __init__(self,
               twin_law,
               miller_array,
               fraction=0.50,
               out=None, verbose=0):
    if out is None:
      out = sys.stdout

    self.fraction = fraction
    miller_array = miller_array.select(
      miller_array.data()>0).set_observation_type(miller_array)

    if not miller_array.is_xray_intensity_array():
      miller_array = miller_array.f_as_f_sq()
    if miller_array.is_real_array():

      acentric_data =  miller_array.select_acentric().set_observation_type(
        miller_array)
      h_test_object  = scaling.h_test(acentric_data.indices(),
                                      acentric_data.data(),
                                      acentric_data.sigmas(),
                                      acentric_data.space_group(),
                                      acentric_data.anomalous_flag(),
                                      twin_law,
                                      fraction)

      self.mean_h = h_test_object.mean_h()
      self.mean_h2 = h_test_object.mean_h2()
      self.estimated_alpha = h_test_object.alpha()
      self.alpha_from_mean_h = (self.mean_h*2.0-1.0)/-2.0
      self.h_array = h_test_object.h_array()
      self.h_values = h_test_object.h_values()
      self.cumul_obs = h_test_object.h_cumul_obs()
      self.cumul_fit = h_test_object.h_cumul_fit()
      if verbose > 0:
        self.show(out)

  def show(self,out=None):
    if out is None:
      out = sys.stdout

    print >> out
    print >> out
    print >> out,"Results of the H-test on a-centric data: "
    print >> out
    print >> out," (Only %3.1f%% of the strongest twin pairs were used)"\
          %(self.fraction*100.0)
    print >> out
    print >> out,"mean |H| : %4.3f" %(self.mean_h) ,\
          "  (0.50: untwinned; 0.0: 50% twinned)"
    print >> out,"mean H^2 : %4.3f" %(self.mean_h2),\
          "  (0.33: untwinned; 0.0: 50% twinned)"
    print >> out,"Estimation of twin fraction via mean |H|: %4.3f" \
          %(self.alpha_from_mean_h)
    print >> out,"Estimation of twin fraction via cum. dist. of H: %4.3f" \
          %( self.estimated_alpha )
    print >> out


class l_test(object):
  def __init__(self, miller_array,
               parity_h=2.0,
               parity_k=2.0,
               parity_l=2.0,
               out=None,verbose=0):
    if out is None:
      out=sys.stdout

    acentric_data = miller_array.select_acentric().set_observation_type(
      miller_array)
    if not miller_array.is_xray_intensity_array():
      acentric_data = acentric_data.f_as_f_sq()
    self.parity_h = parity_h
    self.parity_k = parity_k
    self.parity_l = parity_l

    l_stats = scaling.l_test( acentric_data.indices(),
                              acentric_data.data()/\
                               acentric_data.epsilons().data().as_double(),
                              acentric_data.space_group(),
                              acentric_data.anomalous_flag(),
                              parity_h,
                              parity_k,
                              parity_l,
                              8);

    self.mean_l = l_stats.mean_l()
    self.mean_l2 = l_stats.mean_l2()

    self.l_cumul = l_stats.cumul()
    self.l_values = flex.double(range(self.l_cumul.size()))/float(
      self.l_cumul.size())
    self.l_cumul_untwinned = self.l_values
    self.l_cumul_perfect_twin = self.l_values*(
      3.0-self.l_values*self.l_values)/2.0

    self.ml_alpha = l_stats.ml_alpha()
    if verbose > 0:
      self.show(out)

  def show(self,out=None):
    if out is None:
      out=sys.stdout
    print >> out
    print >> out
    print >> out," L test for acentric data"
    print >> out
    print >> out, " using difference vectors (dh,dk,dl) of the form: "
    print >> out, "(%ihp,%ikp,%ilp)"%(self.parity_h,self.parity_k,self.parity_l)
    print >> out, "  where hp, kp, and lp are random signed integers such that "
    print >> out, "  2 <= |dh| + |dk| + |dl| <= 8 "
    print >> out
    print >> out, "  Mean |L|   :%4.3f  (untwinned: 0.500; perfect twin: 0.375)"\
          %(self.mean_l)
    print >> out, "  Mean  L^2  :%4.3f  (untwinned: 0.333; perfect twin: 0.200)"\
          %(self.mean_l2)
    print >> out
    print >> out, "  The distribution of |L| values indicates a twin fraction of"
    print >> out, "  %3.2f. Note that this estimate is not as reliable as obtained"\
          %(self.ml_alpha)
    print >> out,"  via a Britton plot or H-test if twin laws are available. "
    print >> out
    print >> out




class twin_law_dependend_twin_tests(object):
  """Twin law dependent test results"""
  def __init__(self,
               twin_law,
               miller_array,
               out=None,
               verbose=0,
               miller_calc=None):

    acentric_data = miller_array.select_acentric().set_observation_type(
      miller_array)

    self.twin_law = twin_law

    self.h_test = h_test(twin_law.operator.as_double_array()[0:9],
                         miller_array = acentric_data,
                         out=out,
                         verbose=verbose)

    self.britton_test = britton_test(twin_law.operator.as_double_array()[0:9],
                                     acentric_data,
                                     out=out,
                                     verbose=verbose)

    self.r_values = r_values(miller_array,
                             twin_law.operator.as_double_array()[0:9],
                             miller_calc,
                             out)




class twin_results_interpretation(object):
  def __init__(self,
               nz_test,
               wilson_ratios,
               l_test,
               translational_pseudo_symmetry=None,
               twin_law_related_test=None,
               symmetry_issues=None,

               maha_l_cut=3.5,
               patterson_p_cut=0.01,

               out=None):

    self.maha_l_cut = maha_l_cut
    self.patterson_p_cut = patterson_p_cut


    self.twin_results = twin_results_summary()
    # patterson analyses
    if translational_pseudo_symmetry is not None:
      self.twin_results.patterson_height = translational_pseudo_symmetry.high_peak
      self.twin_results.patterson_p_value = translational_pseudo_symmetry.high_p_value

    # wilson statistics moments, etc
    self.twin_results.i_ratio = wilson_ratios.acentric_i_ratio
    self.twin_results.f_ratio = wilson_ratios.acentric_f_ratio
    self.twin_results.e_sq_minus_1 = wilson_ratios.acentric_abs_e_sq_minus_one

    # l test
    self.twin_results.l_mean = l_test.mean_l
    self.twin_results.l_sq_mean = l_test.mean_l2
    self.compute_maha_l()

    # twin dependent tests
    self.twin_results.n_twin_laws = len(twin_law_related_test)
    for twin_item in twin_law_related_test:
      self.twin_results.twin_laws.append(
        twin_item.twin_law.operator.r().as_hkl() )
      print twin_item.twin_law.twin_type[0]
      self.twin_results.twin_law_type.append(
        str(twin_item.twin_law.twin_type) )

      self.twin_results.r_obs.append(
        twin_item.r_values.r_abs_obs)
      self.twin_results.r_calc.append(
        twin_item.r_values.r_abs_calc)

      self.twin_results.britton_alpha.append(
        twin_item.britton_test.estimated_alpha)
      self.twin_results.h_alpha.append(
        twin_item.h_test.estimated_alpha)

    if self.twin_results.n_twin_laws > 0 :
      self.twin_results.most_worrysome_twin_law = flex.max_index(
        flex.double(self.twin_results.britton_alpha) )
      self.twin_results.suspected_point_group = symmetry_issues.pg_choice
      self.twin_results.possible_sgs = symmetry_issues.sg_possibilities
      self.twin_results.input_point_group = symmetry_issues.pg_low_prim_set_name

    else:
      self.twin_results.input_point_group = None
      self.twin_results.most_worrysome_twin_law = None
      self.twin_results.suspected_point_group = None
      self.twin_results.possible_sgs = None


    # These items will hold the 'verdicts'
    self.patterson_verdict = StringIO()
    self.patterson_short = None

    self.twinning_verdict = StringIO()
    self.twinning_short = None

    self.space_group_issues_verdict = StringIO()
    self.space_group_issues_short = None

    if translational_pseudo_symmetry is not None:
      self.analyse_pseudo_translational_symmetry()
    self.analyse_intensity_stats()


    self.twin_results.verdict = self.patterson_verdict.getvalue() + self.twinning_verdict.getvalue()
    self.twin_results.show(out=out)


  def compute_maha_l(self):
    maha_l = 117820.0
    maha_l2 = 106570
    maha_ll2= -212319
    maha_mean_l = 0.487758242
    mama_mean_l2 = 0.322836996
    tmp_l = self.twin_results.l_mean - maha_mean_l
    tmp_l2 = self.twin_results.l_sq_mean - mama_mean_l2
    maha_distance_l = tmp_l*tmp_l*maha_l +\
                      tmp_l2*tmp_l2*maha_l2 +\
                      tmp_l*tmp_l2*maha_ll2
    maha_distance_l = math.sqrt(maha_distance_l)
    self.twin_results.maha_l = maha_distance_l


  def analyse_pseudo_translational_symmetry(self):
    if self.twin_results.patterson_p_value <= self.patterson_p_cut:
      print >> self.patterson_verdict,\
      "The analyses of the Patterson function reveals a significant off-origin"
      print >> self.patterson_verdict,\
      "peak that is %3.2f %s of the origin peak, indicating pseudo translational symmetry."%(
        self.twin_results.patterson_height,"%")
      print >> self.patterson_verdict,\
      "The chance of finding a peak of this or larger height by random in a "
      print >> self.patterson_verdict,\
      "structure without pseudo translational symmetry is equal to the %5.4e."%(self.twin_results.patterson_p_value)
      if self.twin_results.i_ratio > 2:
        print >> self.patterson_verdict,\
        "The detected tranlational NCS is most likely also responsible for the elevated intensity ratio."
      print >> self.patterson_verdict,\
        "See the relevant section of the logfile for more details."

    else:
      print >> self.patterson_verdict,\
       "The largest off-origin peak in the Patterson function is %3.2f%s of the "%(self.twin_results.patterson_height,"%")
      print >> self.patterson_verdict,\
            "height of the origin peak. No significant pseudotranslation is detected."
      print >> self.patterson_verdict

  def analyse_intensity_stats(self):
    if self.twin_results.maha_l >= self.maha_l_cut:
      if self.twin_results.l_mean < 0.5 :
        print >> self.twinning_verdict, \
          "The results of the L-test indicate that the intensity statistics"
        print >> self.twinning_verdict, \
          "are significantly different then is expected from good to reasonable,"
        print >> self.twinning_verdict, \
          "untwinned data."
        if self.twin_results.n_twin_laws > 0:
          print >> self.twinning_verdict, \
            "As there are twin laws possible given the crystal symmetry, twinning could"
          print >> self.twinning_verdict, \
            "be the reason for the departure of the intensity statistics from normality."
          print >> self.twinning_verdict, \
            "It might be worthwhile carrying refinement with a twin specific target function."
          self.twinning_short=True
          if not self.twin_results.input_point_group == self.twin_results.suspected_point_group:
            print >> self.twinning_verdict
            print >> self.twinning_verdict,\
              "Note that the symmetry of the intensities suggest that the assumed space group"
            print >> self.twinning_verdict, \
              "is too low. As twinning is however suspected, it is not immediuatly clear if this"
            print >> self.twinning_verdict, \
              "is the case. Carefull reprocessing and (twin)refinement for all cases might resolve"
            print >> self.twinning_verdict, \
              "this question."

        if self.twin_results.n_twin_laws == 0:
          print >> self.twinning_verdict, \
            "As there are no twin laws possible given the crystal symmetry, there could be"
          print >> self.twinning_verdict, \
            "a number of reasons for the departure of the intensity statistics from normality."
          print >> self.twinning_verdict, \
            "Overmerging pseudo-symmetric or twinned data, intenisty to amplitude conversion problems"
          print >> self.twinning_verdict, \
            " as well as bad data quality might be possible reasons."
          print >> self.twinning_verdict, \
            "It could be worthwhile considering reprocessing the data."
          self.twinning_short=None

      else:
        self.twinning_short=None
        print >> self.twinning_verdict, \
           "The results of the L-test indicate that the intensity statistics"
        print >> self.twinning_verdict, \
           "Show more centric character then is expected for acentric data."
        if self.twin_results.patterson_p_value <= self.patterson_p_cut:
          print >> self.twinning_verdict, \
            "This behavoir might be explained by the presence of the detected pseudo translation."
          self.twinning_short=False

    else:
      self.twinning_short=False
      print >> self.twinning_verdict, \
        "The results of the L-test indicate that the intensity statistics"
      print >> self.twinning_verdict, \
        "behave as expected. No twinning is suspected."
      if not (self.twin_results.input_point_group == self.twin_results.suspected_point_group):
        print >> self.twinning_verdict, \
          "The symmetry of the lattice and intensity however suggests that the"
        print >> self.twinning_verdict, \
              "input space group is too low. See the relevant sections of the log"
        print >> self.twinning_verdict, \
          "file for more details on your choice of space groups."

      if self.twin_results.n_twin_laws > 0:
        if (self.twin_results.input_point_group == self.twin_results.suspected_point_group):
          print >> self.twinning_verdict, \
            "Even though no twinning is suspected, it might be worthwhile carrying out "
          print >> self.twinning_verdict, \
            "a refinement using a dedicated twin target anyway, as twinned structures with"
          print >> self.twinning_verdict, \
            "low twin fractions are difficult to distinguish from non-twinned structures."
          print >> self.twinning_verdict
          if self.twin_results.most_worrysome_twin_law != None:
            if self.twin_results.britton_alpha[ self.twin_results.most_worrysome_twin_law ]> 0.05:
              print >> self.twinning_verdict,\
                    "The correlation between the intensities related by the twin law %s with an"%(
                self.twin_results.twin_laws[ self.twin_results.most_worrysome_twin_law ])
              print >> self.twinning_verdict,\
                    "estimated twin fraction of %3.2f %s "%(
                self.twin_results.britton_alpha[ self.twin_results.most_worrysome_twin_law],
                "%")
              print >> self.twinning_verdict,\
                    "is most likely due to an NCS axis parallel to the twin axis. This can be verified by"
              print >> self.twinning_verdict,\
                    "supplying calculated data as well."
        else:
          print >> self.twinning_verdict,\
            "As the symmetry is suspected to be incorrect, it is advicable to reconsider data processing."
          print >> self.twinning_verdict,\
            " "










class twin_results_summary(object):
  def __init__(self):
    self.patterson_p_value=None
    self.patterson_height=None

    self.i_ratio=None
    self.f_ratio=None
    self.e_sq_minus_1=None
    self.l_mean=None
    self.l_sq_mean=None
    self.maha_l=None

    self.n_twin_laws=None
    self.twin_laws=[]
    self.twin_law_type=[]
    self.r_obs=[]
    self.r_calc=[]
    self.britton_alpha=[]
    self.h_alpha=[]

    self.most_worrysome_twin_law = None

    self.input_point_group = None
    self.suspected_point_group = None
    self.possible_sgs = []

    self.table = None

    self.verdict=None
    self.in_short=None

  def make_sym_op_table(self):
    if self.r_calc[0]==None:
      legend = ('Operator',
                'type',
                'R obs.',
                'Britton alpha',
                'H alpha')

      table_data = []
      for item in range( len(self.twin_laws) ):
        tmp = [ self.twin_laws[item],
                self.twin_law_type[item],
                str("%4.3f"%(self.r_obs[item])),
                str("%4.3f"%(self.britton_alpha[item])),
                str("%4.3f"%(self.h_alpha[item])) ]
        table_data.append( tmp )


    else:
      legend = ('Operator',
                 'type',
                 'R_abs obs.',
                 'R_abs calc.',
                 'Britton alpha',
                 'H alpha')

      table_data = []
      for item in range( len(self.twin_laws) ):
        tmp = [ self.twin_laws[item],
                self.twin_law_type[item],
                str("%4.3f"%(self.r_obs[item])),
                str("%4.3f"%(self.r_calc[item])),
                str("%4.3f"%(self.britton_alpha[item])),
                str("%4.3f"%(self.h_alpha[item])) ]
        table_data.append( tmp )




    self.table = table_utils.format( [legend]+table_data,
                                comments=None,
                                has_header=True,
                                separate_rows=False,
                                prefix='| ',
                                postfix=' |')




  def show(self,out=None):
    if out is None:
      out = sys.stdout
    print >> out
    print >> out
    print >> out
    print >> out
    print >> out, "-------------------------------------------------------------------------------"
    print >> out, "Twinning and intensity statistics summary (acentric data):"
    print >> out
    print >> out, "Statistics independent of twin laws"
    print >> out, "  - <I^2>/<I>^2 : %5.3f"%(self.i_ratio)
    print >> out, "  - <F>^2/<F^2> : %5.3f"%(self.f_ratio)
    print >> out, "  - <|E^2-1|>   : %5.3f"%(self.e_sq_minus_1)
    print >> out, "  - <|L|>, <L^2>: %5.3f, %4.3f"%(self.l_mean,self.l_sq_mean)
    print >> out, "       Multivariate Z score L-test: %5.3f "%( self.maha_l )
    print >> out, "       The multivariate Z score is a quality measure of the given"
    print >> out, "       spread in intensities. Good to reasonable data is expected"
    print >> out, "       to have a Z score lower than 3.5. "
    print >> out, "       Large values can indicate twinning, but small values do not"
    print >> out, "       neccesarily exclude it. "
    print >> out
    print >> out
    if len(self.twin_laws)>0:
      print >> out, "Statistics depending on twin laws"
      self.make_sym_op_table()
      print >> out, self.table
    else:
      print >> out, "No (pseudo)merohedral twin laws were found."
      print >> out


    if self.patterson_height is not None:
      print >> out
      print >> out, "Patterson analyses"
      print >> out, "  - Largest peak height   : %5.3f"%(self.patterson_height)
      print >> out, "   (correpsonding p value : %8.3e)"%(self.patterson_p_value)
      print >> out
    print >> out
    print >> out, self.verdict
    print >> out, "-------------------------------------------------------------------------------"



class symmetry_issues(object):
  def __init__(self,
               miller_array,
               max_delta=3.0,
               r_cut=0.05,
               scoring_function=None,
               out=None):

    self.out = out
    if self.out == None:
      self.out = sys.stdout

    if scoring_function == None:
      self.scoring_function=[ 0.08, 75.0, 0.08, 75.0 ]

    self.miller_array = miller_array
    self.xs_input = crystal.symmetry(miller_array.unit_cell(),
                                     space_group=miller_array.space_group() )
    self.explore_sg = pointgroup_tools.space_group_graph_from_cell_and_sg(
      self.xs_input.unit_cell(),
      self.xs_input.space_group(),
      max_delta=max_delta)

    self.pg_of_input_sg = self.xs_input.space_group().build_derived_group(
      False,False)

    self.pg_input_name = str(sgtbx.space_group_info(group=self.pg_of_input_sg))

    self.pg_low_prim_set = self.explore_sg.pg_low_prim_set
    self.pg_low_prim_set_name = str(sgtbx.space_group_info(group=self.explore_sg.pg_low_prim_set))

    self.pg_lattice_name = str( sgtbx.space_group_info( group = self.explore_sg.pg_high ) )

    self.miller_niggli =self.miller_array.change_basis(
      self.xs_input.change_of_basis_op_to_niggli_cell()
      )

    self.ops_and_r_pairs = {}
    self.pg_r_used_table = {}
    self.pg_max_r_used_table = {}
    self.pg_r_unused_table = {}

    self.pg_r_unused_split = {}
    self.pg_r_used_split = {}

    self.pg_min_r_unused_table = {}
    self.pg_scores = {}
    self.pg_choice = None
    self.sg_possibilities = []

    self.legend = None
    self.table_data = []
    self.table = None

    self.make_r_table()
    self.make_pg_r_table()
    self.score_all()
    self.wind_up()
    self.show(self.out)


  def make_r_table(self):
    tmp_buffer=StringIO()
    # please find all missing sym ops
    start = str(sgtbx.space_group_info(group=self.pg_low_prim_set))
    tmp_key = self.explore_sg.pg_graph.graph.edge_objects[ start ].keys()[0]
    tmp_edge = self.explore_sg.pg_graph.graph.edge_objects[ start ][tmp_key]

    for symop in tmp_edge.return_used():
      # get the r value please for this symop
      r_value = r_values( self.miller_niggli,
                          symop.as_double_array()[0:9],
                          out=tmp_buffer)
      top = r_value.r_abs_top_obs
      bottom = r_value.r_abs_bottom_obs
      self.ops_and_r_pairs.update( {symop.r().as_hkl():
                                    (top, bottom)} )

    for symop in tmp_edge.return_unused():
      r_value = r_values( self.miller_niggli,
                          symop.as_double_array()[0:9],
                          out=tmp_buffer)
      top = r_value.r_abs_top_obs
      bottom = r_value.r_abs_bottom_obs
      self.ops_and_r_pairs.update( {symop.r().as_hkl():
                                    (top, bottom)} )

  def score_all(self):
    # loop over all pg's
    for pg in self.pg_max_r_used_table:
      score_used = 0
      score_unused = 0
      tmp_const=1e-250
      if len( self.pg_r_used_split[pg] ) > 0 :
        min_used = flex.min( flex.double(self.pg_r_used_split[pg]) )
        max_used = flex.max( flex.double(self.pg_r_used_split[pg]) )
        score_used += math.log( max( 0.5*( 1.0-math.tanh((
          min_used -
          self.scoring_function[0])*self.scoring_function[1] )),tmp_const))
        score_used += math.log( max(0.5*( 1.0-math.tanh((
          max_used -
          self.scoring_function[0])*self.scoring_function[1] )),tmp_const))



      if len( self.pg_r_unused_split[pg] ) > 0 :
        min_unused = flex.min( flex.double(self.pg_r_unused_split[pg]) )
        max_unused = flex.max( flex.double(self.pg_r_unused_split[pg]) )

        score_unused += math.log( max(0.5*( 1.0-math.tanh((
          -max_unused +
          self.scoring_function[0])*self.scoring_function[1] )),tmp_const))
        score_unused += math.log( max(0.5*( 1.0-math.tanh((
          -min_unused +
          self.scoring_function[0])*self.scoring_function[1] )),tmp_const))

      final_score = -score_used -score_unused
      self.pg_scores.update( {pg:final_score} )



  def make_pg_r_table(self):
    start = str(sgtbx.space_group_info(group=self.pg_low_prim_set))
    end = str(sgtbx.space_group_info(group=self.explore_sg.pg_high))

    for pg in self.explore_sg.pg_graph.graph.node_objects:
      r, max_r, min_r, all_r = self.get_r_value_total(start,pg)
      self.pg_r_used_table.update( {pg:r} )
      self.pg_max_r_used_table.update( {pg:max_r} )
      self.pg_r_used_split.update( {pg:all_r} )

      r, max_r, min_r, all_r = self.get_r_value_total(pg,end)
      self.pg_r_unused_table.update( {pg:r} )
      self.pg_min_r_unused_table.update( {pg:min_r} )
      self.pg_r_unused_split.update( {pg:all_r} )


  def get_r_value_total(self,
                        start_pg,
                        end_pg):
    top = 0.0
    bottom = 0.0
    result=None
    max_r=-100.0
    min_r= 10000.0
    # please find the shortest path from start to pg
    path = self.explore_sg.pg_graph.graph.find_shortest_path(
      start_pg,end_pg)
    all_r = []
    if len(path)>1:
      for ii in xrange(len(path)-1):
        start_point = path[ii]
        end_point = path[ii+1]
        # get the ops for this transformation please
        tmp_edge = self.explore_sg.pg_graph.graph.edge_objects[start_point][end_point]
        tmp_ops_used = tmp_edge.symops_used
        tmp_ops_unused = tmp_edge.symops_unused
        # for each used symop find the entry in the table
        # in the same row, we have to find an entry that is in the r table

        for this_s in tmp_ops_used:
          for trial_coset in self.explore_sg.coset_table:
            found_it = False
            for trial_s in trial_coset:
              if trial_s == str(this_s.r().as_hkl()):
                found_it = True
            if found_it:
              for trial_s in trial_coset:
                if self.ops_and_r_pairs.has_key( trial_s ):
                  # key is there, update please
                  top+=self.ops_and_r_pairs[ trial_s ][0]
                  bottom+=self.ops_and_r_pairs[ trial_s ][1]
                  current_r = 0
                  if self.ops_and_r_pairs[ trial_s ][1]>0:
                    current_r = self.ops_and_r_pairs[ trial_s ][0]/self.ops_and_r_pairs[ trial_s ][1]
                    all_r.append(current_r)
                  if (current_r> max_r):
                    max_r = current_r
                  if current_r <= min_r:
                    min_r = current_r

      if bottom==0:
        result = 0
      else:
        result = top/bottom
    if max_r<0:
      max_r = None
    if min_r > 100.0:
      min_r = None
    return result, max_r, min_r, all_r

  def string_it(self, xin, format):
    if xin==None:
      xin='None'
    else:
      xin=format%(xin)
    return xin

  def wind_up(self):
    self.legend = ('Point group',
                   'mean R_used',
                   'max R_used',
                   'mean R_unused',
                   'min R_unused',
                   'choice')

    self.table_data = []
    min_score = 2e+9
    for pg in self.pg_scores:

      tmp = [pg,
             self.string_it(self.pg_r_used_table[pg], "%4.3f"),
             self.string_it(self.pg_max_r_used_table[pg],"%4.3f"),
             self.string_it(self.pg_r_unused_table[pg],"%4.3f"),
             self.string_it(self.pg_min_r_unused_table[pg],"%4.3f"),
             "    " ]
      self.table_data.append( tmp )
      if self.pg_scores[ pg ] < min_score:
        min_score = self.pg_scores[ pg ]
        self.pg_choice = pg

    for row in self.table_data:
      if self.pg_choice == row[0]:
        row[ 5 ] = "<---"

    self.table = table_utils.format([self.legend]+self.table_data,
                                    comments=None,
                                    has_header=True,
                                    separate_rows=False,
                                    prefix='| ',
                                    postfix=' |')
    # store the possible spacegroups and the change of basis ops
    # for the most likely pg
    for xs in self.explore_sg.pg_graph.graph.node_objects[
      self.pg_choice ].allowed_xtal_syms:
      self.sg_possibilities.append( xs )



  def show(self, out=None):
    if out == None:
      out = sys.stdout
    print >> out
    print >> out
    print >> out, "Exploring higher metric symmetry"
    print >> out
    print >> out, "Point group of data as discted by the space group is", self.pg_input_name
    print >> out, "  the point group in the niggli setting is", sgtbx.space_group_info( group=self.pg_low_prim_set )
    print >> out, "The point group of the lattice is", self.pg_lattice_name
    print >> out, "A summary of R values for various possible point groups follow."
    print >> out
    print >> out, self.table
    print >> out
    print >> out, "R_used: mean and maximum R value for symmetry operators *used* in this point group"
    print >> out, "R_unused: mean and minimum R value for symmetry operators *not used* in this point group"
    print >> out, " The likely point group of the data is: ",self.pg_choice
    print >> out
    print >> out, "Possible space groups in this point groups are:"
    for sg in self.sg_possibilities:
      sg[0].show_summary(f=out, prefix= "   ")
      print >> out
    print >> out, "Note that this analyses does not take into account the effects of twinning."
    print >> out, "If the data is (allmost) perfectly twinned, the symmetry will appear to be"
    print >> out, "higher than it actually is."
    print >> out




class r_values(object):
  def __init__(self,
               miller_obs,
               twin_law,
               miller_calc=None,
               out=None,
               n_reflections=400):

    self.obs = miller_obs.deep_copy()
    self.calc=None

    self.n_reflections=n_reflections

    if miller_calc is not None:
      self.calc = miller_calc.deep_copy()
      self.obs, self.calc = self.obs.common_sets( self.calc )
    self.twin_law = twin_law
    self.out=out
    self.rvsr_interpretation=None
    self.table=None
    if out is None:
      self.out=sys.stdout

    self.r_abs_obs=None
    self.r_abs_calc=None
    self.r_sq_obs=None
    self.r_sq_calc=None

    self.r_abs_top_obs=None
    self.r_abs_bottom_obs=None
    self.r_sq_top_obs = None
    self.r_sq_bottom_obs = None

    self.d_star_sq = []
    self.r_sq_obs_reso = []
    self.r_sq_calc_reso = []

    self.guess=None
    #self.resolution_dependent_r_values()
    self.r_vs_r( self.obs,self.calc)
    self.r_vs_r_classification()

    self.show()

  def resolution_dependent_r_values(self):
    # bin the data please in the specified number of bins
    self.obs.setup_binner(reflections_per_bin=self.n_reflections)

    # now we have to loop over all bins please
    for bin_number in self.obs.binner().range_used():
      selection =  self.obs.binner().selection( bin_number ).iselection()
      tmp_obs = self.obs.select( selection )
      tmp_calc=None
      if self.calc is not None:
        tmp_calc = self.calc.select( selection )
      else:
        self.r_sq_calc_reso = None

      self.r_vs_r(tmp_obs, tmp_calc)


      up = self.obs.binner().bin_d_range(bin_number)[0]
      down = self.obs.binner().bin_d_range(bin_number)[1]

      self.d_star_sq.append( 0.5/(up*up) + 0.5/(down*down) )
      self.r_sq_obs_reso.append( self.r_sq_obs )
      if self.calc is not None:
        self.r_sq_calc_reso.append(self.r_sq_calc )



  def r_vs_r(self, input_obs, input_calc):
    # in order to avoid certain issues, take common sets please
    obs = None
    calc = None
    if input_calc is not None:
      obs, calc = input_obs.common_sets( input_calc )
    else:
      obs = input_obs.deep_copy()
    # make sure we have intensities
    if not obs.is_xray_intensity_array():
      obs = obs.f_as_f_sq()
    if calc is not None:
      if not calc.is_xray_intensity_array():
        calc = calc.f_as_f_sq()
    obs_obj = scaling.twin_r( obs.indices(),
                              obs.data(),
                              obs.space_group(),
                              obs.anomalous_flag(),
                              self.twin_law )

    self.r_abs_obs = obs_obj.r_abs_value()
    self.r_sq_obs = obs_obj.r_sq_value()

    self.r_sq_top_obs, self.r_sq_bottom_obs = obs_obj.r_sq_pair()
    self.r_abs_top_obs, self.r_abs_bottom_obs = obs_obj.r_abs_pair()


    if calc is not None:
      calc_obj =  scaling.twin_r( calc.indices(),
                                  calc.data(),
                                  calc.space_group(),
                                  calc.anomalous_flag(),
                                  self.twin_law )
      self.r_abs_calc = calc_obj.r_abs_value()
      self.r_sq_calc = calc_obj.r_sq_value()




  def r_vs_r_classification(self):
    self.rvsr_interpretation = [
      [ "0.1", "0.1", "Misspecified (too low) crystal symmetry \n ",None,""] ,
      [ "0.1", "0.3", "Twin with NCS parallel to twin operator \n ",None,""],
      [ "0.1", "0.5", "Close to perfect twinning \n ",None,""],
      [ "0.3", "0.3", "Data is not twinned, but NCS is parallel \n   to putative twin operator",None,""],
      [ "0.3", "0.5", "Partially twinned data \n ",None,""],
      [ "0.5", "0.5", "No twinning or parallel NCS",None,""]
      ]
    guess = None
    min=1.0
    count=0

    if self.r_abs_calc is not None:
      for test_class in  self.rvsr_interpretation:
        tmp_obs = (float(test_class[0])-self.r_abs_obs)**2.0
        tmp_calc = (float(test_class[1])-self.r_abs_calc)**2.0

        d =  math.sqrt( tmp_obs+tmp_calc )
        test_class[3]="%4.3f"%(d)
        if d < min:
          min = d
          guess = count
        count+=1
    self.guess=guess
    if guess is not None:
      self.rvsr_interpretation[ guess ][4]="<---"

  def show(self):
    print >> self.out
    print >> self.out, "R vs R statistic:"
    print >> self.out, "  R_abs_twin = <|I1-I2|>/<|I1+I2|>"
    print >> self.out, "  Lebedev, Vagin, Murshudov. Acta Cryst. (2006). D62, 83-95"
    print >> self.out
    print >> self.out, "   R_abs_twin observed data   : %4.3f"%(self.r_abs_obs)
    if self.r_abs_calc is not None:
      print >> self.out , "   R_abs_twin calculated data : %4.3f"%(self.r_abs_calc)

    print >> self.out
    print >> self.out, "  R_sq_twin = <(I1-I2)^2>/<(I1+I2)^2>"
    print >> self.out , "   R_sq_twin observed data    : %4.3f"%(self.r_sq_obs)
    if self.r_abs_calc is not None:
      print >> self.out , "   R_sq_twin calculated data  : %4.3f"%(self.r_sq_calc)
      print >> self.out
      print >> self.out , "  The following table can be used in a simple classification scheme"
      print >> self.out,  "  to determine the presence of twinning and or presence of pseudo symmetry."
      print >> self.out
      ## make a neat table please
      legend=('R_abs_obs','R_abs_calc', 'Class', 'distance', 'choice')
      self.table = table_utils.format( [legend]+self.rvsr_interpretation,
                                       comments=None,
                                       has_header=True,
                                       separate_rows=False,
                                       prefix='| ',
                                       postfix=' |')
      print >> self.out, self.table
      print >> self.out

    else:
      print >> self.out , "  No calculated data available."
      print >> self.out , "  R_twin for calculated data not determined."
    print >> self.out





class twin_analyses(object):
  """ Perform various twin related tests"""
  def __init__(self,
               miller_array,
               d_star_sq_low_limit=None,
               d_star_sq_high_limit=None,
               d_hkl_for_l_test=None,
               normalise=True, ## If normalised is true, normalisation is done
               out=None,
               out_plots = None,
               verbose = 1,
               miller_calc=None):

    ## If resolution limits are not specified
    ## use full resolution limit
    if  d_star_sq_high_limit is None:
      d_star_sq_high_limit = flex.min(miller_array.d_spacings().data())
      d_star_sq_high_limit = d_star_sq_high_limit**2.0
      d_star_sq_high_limit = 1.0/d_star_sq_high_limit
    if  d_star_sq_low_limit is None:
      d_star_sq_low_limit = flex.max(miller_array.d_spacings().data())
      d_star_sq_low_limit = d_star_sq_low_limit**2.0
      d_star_sq_low_limit = 1.0/d_star_sq_low_limit
    if d_hkl_for_l_test is None:
      d_hkl_for_l_test=[2.0,2.0,2.0]

    if out is None:
      out = sys.stdout

    ## sanity check on miller array
    if miller_array.observation_type() is None:
      raise RuntimeError("Observation type unknown")
    if miller_array.is_real_array():
      if miller_array.is_xray_intensity_array():
        miller_array = miller_array.f_sq_as_f()
    else:
      raise RuntimeError("Observations should be a real array.")

    print >> out, "Using data between %4.2f to %4.2f Angstrom."\
          %(math.sqrt(1./d_star_sq_low_limit),
            math.sqrt(1./d_star_sq_high_limit))
    print >> out

    miller_array = miller_array.resolution_filter(
      math.sqrt(1.0/d_star_sq_low_limit),
      math.sqrt(1.0/d_star_sq_high_limit))
    ## Determine possible twin laws
    print >> out, "Determining possible twin laws."
    possible_twin_laws = twin_laws(miller_array)
    possible_twin_laws.show(out=out)
    ##-----------------------------

    self.normalised_intensities = wilson_normalised_intensities(
      miller_array, normalise=normalise, out=out, verbose=verbose)

    ## Try to locat e pseudo translational symm.
    ## If no refls are available at low reso,
    ## an exception is thrown and caught here not to disturb things too much
    self.translation_pseudo_symmetry = None
    try:
      self.translation_pseudo_symmetry = detect_pseudo_translations(
        miller_array,
        out=out, verbose=verbose)
    except Sorry: pass

    centric_cut = self.normalised_intensities.centric

    acentric_cut = self.normalised_intensities.acentric

    self.wilson_moments = wilson_moments(
      acentric_cut,
      centric_cut, out=out, verbose=verbose)

    self.nz_test = n_z_test(
      acentric_cut,
      centric_cut,
      out=out, verbose=verbose)

    self.l_test=None
    if self.translation_pseudo_symmetry is not None:
      self.l_test = l_test(
        acentric_cut,
        self.translation_pseudo_symmetry.mod_h,
        self.translation_pseudo_symmetry.mod_k,
        self.translation_pseudo_symmetry.mod_l,
        out=out, verbose=verbose)
    else:
      self.l_test = l_test(
        acentric_cut,
        2,2,2,
        out=out, verbose=verbose)
    ##--------------------------

    if out_plots is not None:


      ## NZ test
      nz_test_plot  = data_plots.plot_data(
        plot_title = 'NZ test',
        x_label = 'z',
        y_label = 'P(Z>=z)',
        x_data = self.nz_test.z,
        y_data = self.nz_test.ac_obs,
        y_legend = 'Acentric observed',
        comments = 'NZ test, acentric and centric data')
      nz_test_plot.add_data(
        y_data = self.nz_test.ac_untwinned,
        y_legend = 'Acentric untwinned')
      nz_test_plot.add_data(
        y_data = self.nz_test.c_obs,
        y_legend = 'Centric observed')
      nz_test_plot.add_data(
        y_data = self.nz_test.c_untwinned,
        y_legend = 'Centric untwinned')
      data_plots.plot_data_loggraph(nz_test_plot,out_plots)
      ## L test
      l_test_plot  = data_plots.plot_data(
        plot_title = 'L test,acentric data',
        x_label = '|l|',
        y_label = 'P(L>=l)',
        x_data = self.l_test.l_values,
        y_data = self.l_test.l_cumul,
        y_legend = 'Observed',
        comments = 'L test, acentric data')
      l_test_plot.add_data(self.l_test.l_cumul_untwinned,
                           'Acentric theory')

      l_test_plot.add_data(self.l_test.l_cumul_perfect_twin,
                           'Acentric theory, perfect twin')
      data_plots.plot_data_loggraph(l_test_plot,out_plots)
      ##------------------------

    ##--------------------------

    self.n_twin_laws = len(possible_twin_laws.operators)

    self.twin_law_dependent_analyses = []



    for ii in range(self.n_twin_laws):
      print >> out
      print >> out,"---------------------------------------------"
      print >> out," Analysing possible twin law : ", \
            possible_twin_laws.operators[ii].operator.r().as_hkl()
      print >> out,"---------------------------------------------"

      tmp_twin_law_stuff = twin_law_dependend_twin_tests(
        possible_twin_laws.operators[ii],
        miller_array,
        out=out,
        verbose=verbose,
        miller_calc=miller_calc)

      self.twin_law_dependent_analyses.append( tmp_twin_law_stuff )

      ## Plotting section
      ##    Britton plot
      britton_plot = data_plots.plot_data(
        plot_title = 'Britton plot for twin law '\
        + possible_twin_laws.operators[ii].operator.r().as_hkl(),
        x_label = 'alpha',
        y_label = 'percentage negatives',
        x_data = tmp_twin_law_stuff.britton_test.britton_alpha,
        y_data = tmp_twin_law_stuff.britton_test.britton_obs,
        y_legend = 'percentage negatives',
        comments = 'percentage negatives')
      britton_plot.add_data(tmp_twin_law_stuff.britton_test.britton_fit,
                            'fit')
      if out_plots is not None:
        data_plots.plot_data_loggraph(britton_plot,out_plots)
      ##    H test
      h_plot = data_plots.plot_data(
        plot_title = 'H test for possible twin law '\
        +possible_twin_laws.operators[ii].operator.r().as_hkl(),
        x_label = 'H',
        y_label = 'S(H)',
        x_data = tmp_twin_law_stuff.h_test.h_array,
        y_data = tmp_twin_law_stuff.h_test.cumul_obs,
        y_legend = 'Observed S(H)',
        comments = 'H test for Acentric data')
      h_plot.add_data(tmp_twin_law_stuff.h_test.cumul_fit, 'Fitted S(H)')
      if out_plots is not None:
        data_plots.plot_data_loggraph(h_plot,out_plots)


      # now we can check for space group related issues
    self.check_sg = None
    if self.n_twin_laws > 0:
      self.check_sg = symmetry_issues(
        miller_array,
        3.0,
        out=out)

    ##--------------------------
    self.twin_summary = twin_results_interpretation(
      self.nz_test,
      self.wilson_moments,
      self.l_test,
      self.translation_pseudo_symmetry,
      self.twin_law_dependent_analyses,
      self.check_sg,
      out=out)




def twin_analyses_brief(miller_array,
                        cut_off=4.0,
                        out = None,
                        verbose=0):
  """
  A very brief twin analyses and tries to answer the question whether or
  not the data is twinned.
  possible outputs and the meaning:
  - False: data is not twinned
  - True : data does not behave as expected. One possible explanantion
           is twinning
  - None : data does not behave as expected, and might or might not be
           due to twinning.
           Also gives none when something messes up.
  """

  out_tmp = StringIO()
  out_tmp_plot = StringIO()
  twin_results = None
  twinned=None
  try:
    twin_results = twin_analyses(miller_array,
                                 d_star_sq_low_limit=1.0/100.0,
                                 d_star_sq_high_limit=1.0/(0.001**2.0),
                                 out = out_tmp,
                                 out_plots = out_tmp_plot,
                                 verbose=verbose)
  except Sorry, RuntimeError: pass


  if out is None:
    out = sys.stdout
  if twin_results is not None:
    if verbose>0:
      twin_results.twin_summary.twin_results.show(out=out)
    if (twin_results.twin_summary.twin_results.maha_l>cut_off):
      if twin_results.twin_summary.twin_results.l_mean <= 0.48:
        twinned = True
    if (twin_results.twin_summary.twin_results.maha_l<=cut_off):
        twinned = False

  return(twinned)
