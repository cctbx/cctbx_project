from __future__ import division
from cctbx.array_family import flex
from mmtbx import scaling
from mmtbx.scaling import absolute_scaling
from cctbx import uctbx
from cctbx import adptbx
from cctbx import sgtbx
from cctbx.sgtbx import pointgroup_tools
from cctbx import maptbx
from cctbx import crystal
import cctbx.xray
import cctbx.sgtbx.lattice_symmetry
import cctbx.sgtbx.cosets
import scitbx.math
from scitbx import matrix
from scitbx.math import chebyshev_lsq
from scitbx.math import chebyshev_polynome
from scitbx.math import chebyshev_lsq_fit
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry
from libtbx.str_utils import StringIO # XXX: pickle support
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
      symbol = self.twin_in_nig.r().as_xyz() )

    self.niggli_cell = self.xs_niggli.unit_cell()
    self.new_niggli_cell = self.xs_niggli.change_basis( self.cb_op ).unit_cell()

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


## python routines copied from iotbx.reflection_statistics.
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
        except KeyboardInterrupt: raise
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

      self.twin_law_table = [table_labels] + table_rows
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
      print >> out, " The presence of twin laws indicates the following: "
      print >> out, "        The symmetry of the lattice (unit cell) is higher (has more elements) "
      print >> out, "        than the point group of the assigned space group."
      print >> out
      print >> out, " There are four likely scenarios associated with the presence of twin laws: "
      print >> out, "    i.  The assigned space group is incorrect (too low)."
      print >> out, "   ii.  The assigned space group is correct and the data *is not* twinned."
      print >> out, "  iii.  The assigned space group is correct and the data *is* twinned."
      print >> out, "   iv.  The assigned space group is not correct (too low) and at the same time, the data *is* twinned."
      print >> out
      print >> out, " Xtriage tries to distinguish between these cases by inspecting the intensity statistics."
      print >> out, " It never hurts to carefully inspect statistics yourself and make sure that the automated "
      print >> out, " interpretation is correct."
      print >> out

      assert (self.m + self.pm)==len(self.operators)



      coset_table = None
      coset_table = sgtbx.cosets.left_decomposition(self.lattice_group, self.intensity_symmetry.space_group().build_derived_acentric_group()\
            .make_tidy())
      it_works = False

      try:
        tmp = self.lattice_group.change_basis( self.change_of_basis_op_to_niggli_cell.inverse() )
        it_works = True
      except: pass

      if it_works:
        print >> out
        print >> out
        print >> out, "Details of automated twin law derivation"
        print >> out, "----------------------------------------"
        print >> out, "Below, the results of the coset decomposition are given. "
        print >> out, "Each coset represents a single twin law, and all symmetry equivalent twin laws are given."
        print >> out, "For each coset, the operator in (x,y,z) and (h,k,l) notation are given. "
        print >> out, "The direction of the axis (in fractional coordinates), the type and possible offsets are given as well."
        print >> out, "Furthermore, the result of combining a certain coset with the input space group is listed.  "
        print >> out, "This table can be usefull when comparing twin laws generated by xtriage with those listed in lookup tables"
        print >> out, "In the table subgroup H denotes the *presumed intensity symmetry*. Group G is the symmetry of the lattice."
        print >> out

        coset_table.show(out=out, cb_op=self.change_of_basis_op_to_niggli_cell.inverse())

        print >> out
        print >> out, "Note that if group H is centered (C,P,I,F), elements corresponding to centering operators are omitted."
        print >> out, "(This is because internally the calculations are done with the symmetry of the reduced cell)"
        print >> out
        print >> out
        print >> out
        print >> out
        print >> out
        print >> out
      else:
        print >> out
        print >> out, "The present crystal symmetry does not allow to have the its lattice symmetry expressed in the setting desired."
        print >> out, "Becasue of this, a full coset table cannot be produced. Working with the data in the reduced cell will solve this."
        print >> out

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

    self.all = work_array
    self.acentric = work_array.select_acentric().as_intensity_array()
    self.centric = work_array.select_centric().as_intensity_array()

    if (self.acentric.indices().size()<=0):
      raise Sorry("No acentric reflections available. Check your input")

    if verbose > -10:
      print >> out, "Splitting data in centrics and acentrics"
      print >> out, "  Number of centrics  :", self.centric.data().size()
      print >> out, "  Number of acentrics :", self.acentric.data().size()



class detect_pseudo_translations(object):
  def __init__(self,
               miller_array,
               low_limit=10.0,
               high_limit=5.0,
               max_sites=100,
               height_cut=0.0,
               distance_cut=15.0,
               p_value_cut=0.05,
               completeness_cut=0.75,
               cut_radius=3.5,
               min_cubicle_edge=5.0,
               out=None,verbose=0):
    if out is None:
      out=sys.stdout

    if miller_array.is_xray_intensity_array():
      miller_array = miller_array.f_sq_as_f()
    work_array = miller_array.resolution_filter(low_limit,high_limit)
    work_array = work_array.select(work_array.data()>0).set_observation_type(
      miller_array)

    self.space_group = miller_array.space_group()
    self.unit_cell = miller_array.unit_cell()

    if work_array.completeness()<completeness_cut:
      print >> out
      print >> out," WARNING (twin_analysis):"
      print >> out,"  The completeness is only %3.2f between %3.1f and %3.1f A."%(
        work_array.completeness(), low_limit, high_limit)
      print >> out,"  This might not be enough to obtain a good estimate"
      print >> out,"  of the presence or absence of pseudo translational"
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


    self.xs = crystal.symmetry(unit_cell=miller_array.unit_cell(),
                               space_group=miller_array.space_group())

    if everything_okai:


      patterson_map = work_array.patterson_map(
        symmetry_flags=maptbx.use_space_group_symmetry).apply_sigma_scaling()

      peak_search_parameters = maptbx.peak_search_parameters(
        peak_search_level=1,
        interpolate=True,
        min_distance_sym_equiv=1e-4,
        general_positions_only=False,
        effective_resolution=None,
        min_cross_distance=cut_radius,
        min_cubicle_edge=min_cubicle_edge)

      cluster_analysis = patterson_map.peak_search(
        parameters=peak_search_parameters)

      peak_list = cluster_analysis.all(max_clusters=1000)

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
        print >> out, "No Patterson vectors with a length larger than"
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

  def suggest_new_space_groups(self,t_den=144,out=None):
    if out is None:
      out=sys.stdout

    symops = []
    sgs = []
    with_operator = []
    new_cell = []

    for peak in self.suspected_peaks:
      if peak[2] < self.p_value_cut:
        xyz = peak[0]
        dx = self.closest_rational( xyz[0] )
        dy = self.closest_rational( xyz[1] )
        dz = self.closest_rational( xyz[2] )
        additional_symop = None
        if ([dx,dy,dz]).count(None)==0:
          additional_symop = "x%s, y%s, z%s"%(
            dx, dy, dz )
          if additional_symop is not None:
            symops.append( additional_symop )


    tmp_space_group = sgtbx.space_group_info( group=self.space_group )
    try:
      tmp_space_group = sgtbx.space_group_info( str(tmp_space_group), space_group_t_den=t_den)
    except KeyboardInterrupt: raise
    except : pass

    for so in symops:
      sg_str = None
      try:
        new_tmp_sg = tmp_space_group.group().make_tidy()
        smx = sgtbx.rt_mx( so )
        smx = smx.new_denominators( new_tmp_sg.r_den(), new_tmp_sg.t_den() )
        new_tmp_sg.expand_smx( smx )
        sg_str = str( sgtbx.space_group_info( group = new_tmp_sg,  space_group_t_den=t_den  ) )
        to_ref_set = sgtbx.space_group_info(
          group = new_tmp_sg,  space_group_t_den=t_den  ).change_of_basis_op_to_reference_setting()
      except: pass
      if sg_str not in sgs:
        if sg_str is not None:
          sgs.append( sg_str )
          with_operator.append( so )
          new_cell.append( self.unit_cell.change_basis( to_ref_set.c_inv().r().as_double() ) )
        else:
          sgs.append( None )
          new_cell.append( self.unit_cell )
    number_of_new_sgs = len(symops)-sgs.count(None)
    if number_of_new_sgs>0:
      print >> out
      print >> out, " If the observed pseudo translationals are crystallographic"
      print >> out, " the following spacegroups and unit cells are possible: "
      print >> out
      print >> out, " %s                %s         %s  "%("space group", "operator", "unit cell of reference setting")

      for sg,op,uc in zip(sgs,with_operator,new_cell):
        if sg is not None:
          print >> out, " %20s  %20s  (%5.2f, %5.2f, %5.2f,  %5.2f, %5.2f, %5.2f)"%(
            sg,
            op,
            uc.parameters()[0],
            uc.parameters()[1],
            uc.parameters()[2],
            uc.parameters()[3],
            uc.parameters()[4],
            uc.parameters()[5] )
        print >> out


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

  def closest_rational(self, fraction, eps=0.02,return_text=True):
    tmp_fraction = abs(fraction)
    sign = "+"
    if fraction < 0:
      sign = "-"

    num = None
    den = None
    den_list = [2,3,4,5,6]
    fraction_trials = []
    min_del=10.0
    best_frac=None
    for trial_den in den_list:
      start_num=0
      if trial_den > 2:
        start_num = 1
      for trial_num in xrange(start_num,trial_den):
        tmp = (sign, trial_num, trial_den)
        tmp_frac = trial_num/trial_den
        delta = abs(tmp_frac - tmp_fraction )
        if delta < min_del:
          best_frac = tmp
          min_del = float( delta )
    result = None
    if min_del < eps:
      if return_text:
        if best_frac[1]==0:
          result = ""
        else:
          result = "%s%s/%s"%(best_frac[0],best_frac[1],best_frac[2])
      else:
        result = best_frac
    return result


  def show(self,out=None):
    if out is None:
      out = sys.stdout
    print >> out
    print >> out," Patterson analyses"
    print >> out, "------------------"
    print >> out
    print >> out," Largest Patterson peak with length larger than 15 Angstrom "
    print >> out
    self.peak_info = """\
 Frac. coord.        : %8.3f %8.3f %8.3f
""" % (self.high_peak_xyz)
    self.peak_info += """\
 Distance to origin  : %8.3f
 Height (origin=100) : %8.3f
 p_value(height)     : %12.3e
""" % (self.high_peak_distance,self.high_peak, self.high_p_value)
    print >> out, self.peak_info
    print >> out
    self.peak_meaning = """\
   The reported p_value has the following meaning:
     The probability that a peak of the specified height
     or larger is found in a Patterson function of a
     macro molecule that does not have any translational
     pseudo symmetry is equal to %10.3e.
     p_values smaller than 0.05 might indicate
     weak translational pseudo symmetry, or the self vector of
     a large anomalous scatterer such as Hg, whereas values
     smaller than 1e-3 are a very strong indication for
     the presence of translational pseudo symmetry.
""" % self.high_p_value
    print >> out, self.peak_meaning
    print >> out


    if self.high_p_value <= self.p_value_cut:

      print >> out
      print >> out, "The full list of Patterson peaks is: "
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
      self.suggest_new_space_groups(out=out)

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

    self.sign_ac = '+'
    if self.mean_diff_ac < 0:
      self.sign_ac = '-'

    self.sign_c = '+'
    if self.mean_diff_c < 0:
      self.sign_c = '-'
    print >> out,"-----------------------------------------------"
    print >> out,"| Maximum deviation acentric      :  %4.3f    |" \
          %(self.max_diff_ac)
    print >> out,"| Maximum deviation centric       :  %4.3f    |" \
          %(self.max_diff_c)
    print >> out,"|                                             |"
    print >> out,"| <NZ(obs)-NZ(twinned)>_acentric  : %1s%4.3f    |" \
          %(self.sign_ac,math.fabs(self.mean_diff_ac))
    print >> out,"| <NZ(obs)-NZ(twinned)>_centric   : %1s%4.3f    |" \
          %(self.sign_c,math.fabs(self.mean_diff_c))
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

    self.max_iter=1000
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
    icount=0
    while not_done:
      icount+=1
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
      if icount > self.max_iter: # a nasty fail safe
        break
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
    correlation = covar_xy/(math.sqrt(abs(var_x*var_y))+1.0e-6)
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
    print >> out,"Results of the H-test on acentric data: "
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

class ml_murray_rust_with_ncs(object):
  def __init__(self,
               miller_array,
               twin_law,
               out,
               n_bins=10,
               calc_data=None,
               start_alpha=None):
    if out == None:
      out = sys.stdout

    self.twin_law = twin_law

    self.twin_cap = 0.45
    self.d_cap = 0.95

    assert miller_array.is_xray_intensity_array()
    assert miller_array.sigmas() is not None

    tmp_miller_array = miller_array.deep_copy().set_observation_type( miller_array )
    # make a binner
    self.binner = tmp_miller_array.setup_binner( n_bins=n_bins )
    self.out = out
    self.f = None
    self.cycle = 0
    self.ml_object = scaling.ml_twin_with_ncs( tmp_miller_array.data(),
                                               tmp_miller_array.sigmas(),
                                               tmp_miller_array.indices(),
                                               self.binner.bin_indices(),
                                               tmp_miller_array.space_group(),
                                               tmp_miller_array.anomalous_flag(),
                                               twin_law,
                                               tmp_miller_array.unit_cell(),
                                               4 )
    self.n = 1+n_bins
    self.x = flex.double([-3])
    if start_alpha is not None:
      if start_alpha > 0.45:
         start_alpha = 0.40
      if start_alpha <=0 :
         start_alpha = 0.05
      tmp = -math.log(self.twin_cap/start_alpha-1.0)
      self.x = flex.double([tmp])

    for ii in xrange(n_bins):
      self.x.append( -1.0 - ii/20.0 )

    self.message()

    term_parameters = scitbx.lbfgs.termination_parameters( max_iterations = 1000 )
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self, termination_params=term_parameters )


    scitbx.lbfgs.run(target_evaluator=self)
    self.print_it()

    if calc_data is not None:
      self.calc_correlation(  tmp_miller_array, calc_data )

  def string_it(self,x):
    return str("%4.3f"%(x))

  def calc_correlation(self,obs,calc):
    if calc.is_xray_amplitude_array():
      calc = calc.f_as_f_sq()

    calc = calc.common_set( other=obs )
    calc.use_binning_of( obs  )

    print >> self.out , " The correlation of the calculated F^2 should be similar to "
    print >> self.out , " the estimated values. "
    print >> self.out
    print >> self.out , " Observed correlation between twin related, untwinned calculated F^2"
    print >> self.out , " in resolutiuon ranges, as well as ewstimates D_ncs^2 values:"
    # now loop over all resolution bins, get twin related intensities and get twin related intensities please
    print >> self.out, " Bin    d_max     d_min     CC_obs   D_ncs^2 "
    for i_bin in calc.binner().range_used():
      tmp_array = calc.select( calc.binner().bin_indices() == i_bin )
      tmp_r_class = scaling.twin_r( tmp_array.indices(),
                                    tmp_array.data(),
                                    tmp_array.space_group(),
                                    tmp_array.anomalous_flag(),
                                    self.twin_law )
      low = self.binner.bin_d_range(i_bin)[0]
      high = self.binner.bin_d_range(i_bin)[1]
      d_theory_sq = tmp_r_class.correlation()
      d_ncs_sq = self.x[i_bin]
      d_ncs_sq = self.d_cap/( 1.0+math.exp(-d_ncs_sq) )
      d_ncs_sq = d_ncs_sq*d_ncs_sq

      print >> self.out, "%3i)    %5.4f -- %5.4f :: %6s   %6s"%(i_bin,low, high,
                                                              self.string_it(d_theory_sq),
                                                              self.string_it(d_ncs_sq) )




  def compute_functional_and_gradients(self):
    tmp = self.ml_object.p_tot_given_t_and_coeff( self.x[0],
                                                  self.x[1:] )
    f = tmp[0]
    g = tmp[1:]

    self.f=f
    self.cycle+=1
    if self.cycle%5==0:
      print >> self.out, "%3i "%(self.cycle),
      #self.print_it()
    else:
      print >>self.out, ".",
    if self.cycle%30==0:
      print >>self.out
    self.out.flush()
    return f,flex.double(g)

  def message(self):
    print >> self.out
    print >> self.out, "Estimation of twin fraction, while taking into account the"
    print >> self.out, "effects of possible NCS parallel to the twin axis."
    print >> self.out, "    Zwart, Read, Grosse-Kunstleve & Adams, to be published."
    print >> self.out
    print >> self.out, "  A parameters D_ncs will be estimated as a function of resolution,"
    print >> self.out, "  together with a global twin fraction."
    print >> self.out, "  D_ncs is an estimate of the correlation coefficient between"
    print >> self.out, "  untwinned, error-free, twin related, normalized intensities."
    print >> self.out, "  Large values (0.95) could indicate an incorrect point group."
    print >> self.out, "  Value of D_ncs larger than say, 0.5, could indicate the presence"
    print >> self.out, "  of NCS. The twin fraction should be smaller or similar to other"
    print >> self.out, "  estimates given elsewhere."
    print >> self.out
    print >> self.out, "  The refinement can take some time. "
    print >> self.out, "  For numerical stability issues, D_ncs is limited between 0 and 0.95."
    print >> self.out, "  The twin fraction is allowed to vary between 0 and 0.45."
    print >> self.out, "  Refinement cycle numbers are printed out to keep you entertained."
    print >> self.out


  def print_it(self):
    print >> self.out
    print >> self.out
    print >> self.out, "  Cycle : %3i"%(self.cycle)
    print >> self.out, "  -----------"
    print >> self.out, "  Log[likelihood]: %15.3f"%(self.f)
    tf=self.x[0]
    tf=self.twin_cap/(1+math.exp(-tf))
    print >> self.out, "  twin fraction: %4.3f"%( tf )
    print >> self.out, "  D_ncs in resolution ranges: "
    for ii in xrange(self.x.size()-1):
      d_ncs = self.x[ii+1]
      d_ncs = self.d_cap/( 1.0+math.exp(-d_ncs) )
      low = self.binner.bin_d_range(ii+1)[0]
      high = self.binner.bin_d_range(ii+1)[1]
      print >> self.out, "     %5.4f -- %5.4f :: %4.3f"%(low, high, d_ncs)
    print >> self.out


class ml_murray_rust(object):
  def __init__(self,
               miller_array,
               twin_law,
               out,
               n_points = 4):
    if out == None:
      out = sys.stdout

    assert miller_array.is_xray_intensity_array()
    assert miller_array.sigmas() is not None

    ml_murray_rust_object = scaling.ml_murray_rust( miller_array.data(),
                                                    miller_array.sigmas(),
                                                    miller_array.indices(),
                                                    miller_array.space_group(),
                                                    miller_array.anomalous_flag(),
                                                    twin_law,n_points)
    self.twin_fraction = []
    self.nll = []
    print >> out
    print >> out, "Maximum Likelihood twin fraction determination"
    print >> out, "    Zwart, Read, Grosse-Kunstleve & Adams, to be published."
    for ii in xrange(1,23):
      t=ii/46.0
      self.twin_fraction.append( t )
      self.nll.append( - ml_murray_rust_object.fast_log_p_given_t( t )  )
    print >> out

    i_t_max = flex.min_index( flex.double( self.nll ) )
    self.estimated_alpha = None
    if (i_t_max >= 1) and (i_t_max < len(self.nll)-1 ) :
      tmp_t = [0,0,0]
      tmp_nll = [0,0,0]
      tmp_t[0] = self.twin_fraction[ i_t_max - 1 ]
      tmp_t[1] = self.twin_fraction[ i_t_max     ]
      tmp_t[2] = self.twin_fraction[ i_t_max + 1 ]
      tmp_nll[0] = self.nll[ i_t_max - 1 ]
      tmp_nll[1] = self.nll[ i_t_max ]
      tmp_nll[2] = self.nll[ i_t_max + 1 ]

      tmp_top = (tmp_t[2]**2.0)*(tmp_nll[0] - tmp_nll[1]) + \
                (tmp_t[1]**2.0)*(tmp_nll[1] - tmp_nll[2]) + \
                (tmp_t[0]**2.0)*(tmp_nll[2] - tmp_nll[0])

      tmp_bottom = (tmp_t[2])*(tmp_nll[0] - tmp_nll[1]) + \
                   (tmp_t[1])*(tmp_nll[1] - tmp_nll[2]) + \
                   (tmp_t[0])*(tmp_nll[2] - tmp_nll[0])

      self.estimated_alpha= tmp_top/(2.0*tmp_bottom)

      if (self.estimated_alpha < tmp_t[0]) or ( self.estimated_alpha > tmp_t[2]):
        self.estimated_alpha = self.twin_fraction[ i_t_max ]

    else:
      self.estimated_alpha = self.twin_fraction[ i_t_max ]

    print >> out
    print >> out, "   The estimated twin fraction is equal to %4.3f"%(self.estimated_alpha)
    print >> out


class correlation_analyses(object):
  def __init__(self,
               miller_obs,
               miller_calc,
               twin_law,
               d_weight=0.1,
               out=None):
    self.out=out
    if self.out==None:
      self.out = sys.stdout

    self.twin_law=sgtbx.change_of_basis_op( twin_law )
    self.twin_law= self.twin_law.new_denominators( r_den=miller_calc.space_group().r_den(),
                                                   t_den=miller_calc.space_group().t_den() )

    self.obs=miller_obs.deep_copy()
    self.calc=miller_calc.deep_copy()
    # the incomming data is normalized, normalize the calculated data as well please
    normalizer = absolute_scaling.kernel_normalisation(
        self.calc, auto_kernel=True)
    self.calc = normalizer.normalised_miller.deep_copy()
    self.calc = self.calc.select(self.calc.data()>0)
    self.calc_part2 = self.calc.change_basis( self.twin_law )
    self.calc_part2 = self.calc.customized_copy(
      indices = self.calc_part2.indices(),
      data = self.calc_part2.data(),
      ).map_to_asu()


    # make sure we have common sets of everything
    self.calc, self.calc_part2 = self.calc.common_sets( self.calc_part2 )
    self.obs,  self.calc = self.obs.common_sets( self.calc )
    self.obs,  self.calc_part2 = self.obs.common_sets( self.calc_part2 )

    self.d=d_weight

    self.alpha = []
    self.cc = []

    print >> out, "Perfoming correlation analyses"
    print >> out, "  The supplied calculated data are normalized and artificially twinned"
    print >> out, "  Subsequently a correlation with the observed data is computed."
    print >> out

    if self.obs.data().size()>0:
      for ii in xrange(50):
        alpha=ii/100.0
        self.alpha.append( alpha )
        cc = self.twin_the_data_and_compute_cc(alpha)
        self.cc.append(cc)
    self.max_alpha, self.max_cc = self.find_maximum()
    print >> out, "   Results: "
    print >> out, "     Correlation : ", self.max_cc
    print >> out, "     Estimated twin fraction : ", self.max_alpha
    print >> out

  def twin_the_data_and_compute_cc(self, alpha):
    wx = self.obs.sigmas()
    x = self.obs.data()
    dd = self.obs.d_star_sq().data()
    dd = flex.exp( -7.7*self.d*self.d*dd )
    wy = 1.0-dd*dd
    y = (1-alpha)*self.calc.data() + alpha*self.calc_part2.data()
    y = dd*dd*y
    w_tot = wx*wx + wy*wy
    cc = self.weighted_cc( x,y,w_tot )
    return(cc)

  def weighted_cc(self,x,y,w):
    tw = flex.sum( w )
    meanx = flex.sum(w*x)/tw
    meany = flex.sum(w*y)/tw
    meanxx = flex.sum(w*x*x)/tw
    meanyy = flex.sum(w*y*y)/tw
    meanxy = flex.sum(w*x*y)/tw
    top = meanxy-meanx*meany
    bottom = math.sqrt(  (meanxx - meanx*meanx)*(meanyy - meany*meany) )
    if math.fabs(bottom) > 0:
      return top/bottom
    else:
      return 0.0

  def find_maximum(self):
    max_cc = flex.max( flex.double( self.cc ) )
    max_alpha = flex.max_index( flex.double( self.cc ) )
    max_alpha =  self.alpha[ max_alpha ]
    return max_alpha, max_cc




class twin_law_dependend_twin_tests(object):
  """Twin law dependent test results"""
  def __init__(self,
               twin_law,
               miller_array,
               out=None,
               verbose=0,
               miller_calc=None,
               normalized_intensities=None,
               ncs_test=None,
               n_ncs_bins=None):

    if n_ncs_bins is None:
      n_ncs_bins = 7

    acentric_data = miller_array.select_acentric().set_observation_type(
      miller_array)

    self.twin_law = twin_law

    tmp_detwin_object = scaling.detwin(acentric_data.indices(),
                                       acentric_data.data(),
                                       acentric_data.sigmas(),
                                       acentric_data.space_group(),
                                       acentric_data.anomalous_flag(),
                                       twin_law.operator.as_double_array()[0:9] )

    self.twin_completeness = tmp_detwin_object.completeness()
    self.results_available = False
    self.h_test=None
    self.britton_test=None
    self.r_values=None
    self.correlation=None
    self.ml_murray_rust=None
    self.ml_murry_rust_with_ncs=None

    if self.twin_completeness>0:
      self.results_available = True

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

      if miller_calc is not None:
        ## Magic number to control influence of hiugh resolutiuon reflections
        d_weight=0.25
        try:
          self.correlation = correlation_analyses(normalized_intensities,
                                                  miller_calc,
                                                  twin_law.operator,
                                                  d_weight,
                                                  out)
        except KeyboardInterrupt: raise
        except: pass


      if normalized_intensities.sigmas() is not None:
        self.ml_murray_rust = ml_murray_rust( normalized_intensities,
                                              twin_law.operator.as_double_array()[0:9],
                                              out )

        if ncs_test:
          self.ml_murry_rust_with_ncs = ml_murray_rust_with_ncs( normalized_intensities,
                                                                 twin_law.operator.as_double_array()[0:9],
                                                                 out,
                                                                 n_ncs_bins,
                                                                 miller_calc,
                                                                 start_alpha=self.ml_murray_rust.estimated_alpha)



    else:
      print >> out, "The twin completeness is %3.2f. Twin law dependent test not performed."%(self.twin_completeness)
      print >> out







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
      self.twin_results.twin_law_type.append(
        str(twin_item.twin_law.twin_type) )

      if twin_item.results_available:
        self.twin_results.r_obs.append(
          twin_item.r_values.r_abs_obs)
        self.twin_results.r_calc.append(
          twin_item.r_values.r_abs_calc)

        self.twin_results.britton_alpha.append(
          twin_item.britton_test.estimated_alpha)
        self.twin_results.h_alpha.append(
          twin_item.h_test.estimated_alpha)

        if twin_item.ml_murray_rust is not None:
          self.twin_results.murray_rust_alpha.append(
            twin_item.ml_murray_rust.estimated_alpha)
        else:
          self.twin_results.murray_rust_alpha.append(None)

      else:
        self.twin_results.r_obs.append(None)
        self.twin_results.r_calc.append(None)
        self.twin_results.britton_alpha.append(None)
        self.twin_results.h_alpha.append(None)
        self.twin_results.murray_rust_alpha.append(None)

    if self.twin_results.n_twin_laws > 0 :
      if twin_item.results_available:
        max_tf = -1
        max_index = -1
        for ii, tf in enumerate(self.twin_results.britton_alpha):
          if tf is not None:
            if float(tf) > max_tf:
              max_tf = tf
              max_index = ii

        self.twin_results.most_worrysome_twin_law = max_index
      else:
        tmp_tf = []
        for tf in self.twin_results.britton_alpha:
          if tf is None:
            tmp_tf.append( -100 )
          else:
            tmp_tf.append( tf )
        self.twin_results.most_worrysome_twin_law = flex.max_index(
          flex.double(tmp_tf) )

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
    # Mean vector: 0.48513455414 0.316418789809
    # Variancxe /covriance: 0.000143317086266 0.000192434487707 0.000165215146913
    # Its inverse {{679728., -583582.}, {-583582., 506232.}}
    # These numbers have been obtained from roughly 1200 datasets
    # They were selected using
    # strong high resolution limit (over 85% complete, for i/sigi>3.0) is smaller than 2.0
    # no detected pseudo translation
    # <|L|> > 0.45 (limit decided based on histogram, mainly  for removing some outliers)
    maha_l   = 679728.0 #old value 117820.0
    maha_l2  = 506232.0 #old value 106570
    maha_ll2 = 2.0*-583582.0 #old value -212319
    maha_mean_l = 0.48513455414  #old value: 0.487758242
    mama_mean_l2 = 0.316418789809 #old value: 0.322836996
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
          "are significantly different than is expected from good to reasonable,"
        print >> self.twinning_verdict, \
          "untwinned data."
        if self.twin_results.n_twin_laws > 0:
          print >> self.twinning_verdict, \
            "As there are twin laws possible given the crystal symmetry, twinning could"
          print >> self.twinning_verdict, \
            "be the reason for the departure of the intensity statistics from normality."
          print >> self.twinning_verdict, \
            "It might be worthwhile carrying out refinement with a twin specific target function."
          self.twinning_short=True
          if not self.twin_results.input_point_group == self.twin_results.suspected_point_group:
            print >> self.twinning_verdict
            print >> self.twinning_verdict,\
              "Note that the symmetry of the intensities suggest that the assumed space group"
            print >> self.twinning_verdict, \
              "is too low. As twinning is however suspected, it is not immediately clear if this"
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
            "Overmerging pseudo-symmetric or twinned data, intensity to amplitude conversion problems"
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
           "Show more centric character than is expected for acentric data."
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
          if self.twin_results.patterson_p_value <= self.patterson_p_cut:
            print >> self.twinning_verdict,\
                     "Note however that the presence of translational NCS (and possible rotational pseudo"
            print >> self.twinning_verdict,\
                     "symmetry parallel to the twin axis) can make the detection of twinning difficult."
            print >> self.twinning_verdict,\
                     "Trying various space group and twinning hypothesis in structure refinement might provide an answer"
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
    self.murray_rust_alpha=[]

    self.most_worrysome_twin_law = None

    self.input_point_group = None
    self.suspected_point_group = None
    self.possible_sgs = []

    self.table = None

    self.verdict=None
    self.in_short=None

  def string_it(self, x ):
     if x is None:
       return("None")
     else:
       return( str("%4.3f"%(x)) )

  def make_sym_op_table(self):
    if self.r_calc[0]==None:
      legend = ('Operator',
                'type',
                'R obs.',
                'Britton alpha',
                'H alpha',
                'ML alpha')

      table_data = []
      for item in range( len(self.twin_laws) ):
        tmp = [ self.twin_laws[item],
                self.twin_law_type[item],
                self.string_it( self.r_obs[item]),
                self.string_it( self.britton_alpha[item]),
                self.string_it( self.h_alpha[item]),
                self.string_it( self.murray_rust_alpha[item]) ]
        table_data.append( tmp )


    else:
      legend = ('Operator',
                 'type',
                 'R_abs obs.',
                 'R_abs calc.',
                 'Britton alpha',
                 'H alpha',
                 'ML alpha')

      table_data = []
      for item in range( len(self.twin_laws) ):
        tmp = [ self.twin_laws[item],
                self.twin_law_type[item],
                str("%4.3f"%(self.r_obs[item])),
                str("%4.3f"%(self.r_calc[item])),
                str("%4.3f"%(self.britton_alpha[item])),
                str("%4.3f"%(self.h_alpha[item])),
                str("%4.3f"%(self.murray_rust_alpha[item])) ]
        table_data.append( tmp )




    self.table = table_utils.format( [legend]+table_data,
                                comments=None,
                                has_header=True,
                                separate_rows=False,
                                prefix='| ',
                                postfix=' |')
    self.table_data = table_data



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
    self.independent_stats = """\
  <I^2>/<I>^2 : %5.3f
  <F>^2/<F^2> : %5.3f
  <|E^2-1|>   : %5.3f
  <|L|>, <L^2>: %5.3f, %4.3f
  Multivariate Z score L-test: %5.3f
""" % (self.i_ratio, self.f_ratio, self.e_sq_minus_1, self.l_mean,
       self.l_sq_mean, self.maha_l)
    print >> out, self.independent_stats
    self.z_score_info = """\
 The multivariate Z score is a quality measure of the given
 spread in intensities. Good to reasonable data are expected
 to have a Z score lower than 3.5.
 Large values can indicate twinning, but small values do not
 necessarily exclude it.
"""
    print >> out, self.z_score_info
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
      if self.patterson_p_value <  0.001:
        print >> out, "   (corresponding p value : %8.3e)"%(self.patterson_p_value)
      else:
        print >> out, "   (corresponding p value : %6.5f)"%(self.patterson_p_value)
      print >> out
    print >> out
    print >> out, self.verdict
    print >> out, "-------------------------------------------------------------------------------"



class symmetry_issues(object):
  def __init__(self,
               miller_array,
               max_delta=3.0,
               r_cut=0.05,
               sigma_inflation = 1.25,
               out=None):
    self.sigma_warning = None
    if out == None:
      out = sys.stdout

    self.inflate_sigma = sigma_inflation

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
    self.show(out)
    self.xs_with_pg_choice_in_standard_setting = crystal.symmetry( unit_cell   = self.miller_niggli.unit_cell(),
                                                                   space_group = self.pg_choice,
                                                                   assert_is_compatible_unit_cell=False,
                                                                   ).as_reference_setting()


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
      tmp = self.miller_niggli.f_as_f_sq()
      sigmas = tmp.sigmas()
      if tmp.sigmas() is None:
         sigmas = tmp.d_star_sq().data()
         # lets assume a sigma of 5% at 0 and a sigma of 20% at the high resolution
         Berror = -math.log(0.20/0.05)/flex.max(sigmas)
         sigmas = 0.05*flex.exp( -sigmas*Berror )
         sigmas = tmp.data()*sigmas
         self.sigma_warning = """  ----> WARNING: NO SIGMAS FOUND  <----
          Sigmas now modeled as 0.05*|F|*exp( -Berror/d^2 )
          with Berror=%6.3f"""%(Berror)


      merger = cctbx.xray.merger( tmp.indices(),
                                  tmp.data(),
                                  sigmas*self.inflate_sigma ,
                                  sgtbx.space_group_info(pg).group(),
                                  tmp.anomalous_flag(),
                                  tmp.unit_cell() )
      final_score = -merger.bic()
      self.pg_scores.update( {pg:final_score} )



  def make_pg_r_table(self):
    start = str(sgtbx.space_group_info(group=self.pg_low_prim_set))
    end = str(sgtbx.space_group_info(group=self.explore_sg.pg_high))

    for pg in self.explore_sg.pg_graph.graph.node_objects:
      tot_in, max_in, tot_out, min_out, r_in, r_out = self.get_r_value_total(start,pg)
      self.pg_r_used_table.update( {pg:tot_in} )
      self.pg_max_r_used_table.update( {pg:max_in} )
      self.pg_r_used_split.update( {pg:r_in} )
      self.pg_r_unused_table.update( {pg:tot_out} )
      self.pg_min_r_unused_table.update( {pg:min_out} )
      self.pg_r_unused_split.update( {pg:r_out} )

  def get_r_value_total(self,
                        start_pg,
                        end_pg):
    # do a coset decompostion on start and end
    spg = sgtbx.space_group_info( start_pg ).group()
    epg = sgtbx.space_group_info( end_pg ).group()

    cosets = cctbx.sgtbx.cosets.left_decomposition( epg, spg )
    coset_ops = []
    # get the coset operations in a nice list
    for cs in cosets.partitions:
      for op in cs:
        coset_ops.append(  op.r().as_hkl() )

    r_in = flex.double()
    r_out = flex.double()
    tot_r_in  = [0,0]
    tot_r_out = [0,0]

    tmp_keys = self.ops_and_r_pairs.keys()
    tmp_values = self.ops_and_r_pairs.values()
    tmp = self.ops_and_r_pairs.copy()

    for op in coset_ops:
      if op in tmp_keys:  #tmp.has_key( op ):
        # Work around for absence of pop in python 2.2
        iii = tmp_keys.index( op )
        a = tmp_values[ iii ]
        # now we have to pop these guys from the array
        tmp_tmp_keys = []
        tmp_tmp_values = []
        for ii in xrange( len(tmp_keys) ):
          if ii != iii:
            tmp_tmp_keys.append( tmp_keys[ii] )
            tmp_tmp_values.append( tmp_values[ii] )
        tmp_keys = list(tmp_tmp_keys)
        tmp_values = list(tmp_tmp_values)

        r_in.append( a[0]/max(a[1],1e-8) )
        tot_r_in[0]+=a[0]
        tot_r_in[1]+=a[1]
    for a in tmp_values:
      r_out.append( a[0]/max(a[1],1e-8) )
      tot_r_out[0]+=a[0]
      tot_r_out[1]+=a[1]

    tot_in  = None
    max_in  = None
    tot_out = None
    min_out = None

    if r_in.size() > 0:
      tot_in = tot_r_in[0]/max(1e-8,tot_r_in[1] )
      max_in = flex.max( r_in )
    if r_out.size() > 0:
      tot_out = tot_r_out[0]/max( 1e-8,tot_r_out[1] )
      min_out = flex.min(r_out)

    return tot_in, max_in, tot_out, min_out, r_in, r_out


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
                   'BIC',
                   'choice')

    self.table_data = []
    min_score = 2e+9
    for pg in self.pg_scores:

      tmp = [pg,
             self.string_it(self.pg_r_used_table[pg], "%4.3f"),
             self.string_it(self.pg_max_r_used_table[pg],"%4.3f"),
             self.string_it(self.pg_r_unused_table[pg],"%4.3f"),
             self.string_it(self.pg_min_r_unused_table[pg],"%4.3f"),
             self.string_it(self.pg_scores[pg],"%6.3e"),
             "    " ]
      self.table_data.append( tmp )
      if self.pg_scores[ pg ] < min_score:
        min_score = self.pg_scores[ pg ]
        self.pg_choice = pg

    for row in self.table_data:
      if self.pg_choice == row[0]:
        row[ 6 ] = "<---"

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

  def return_point_groups(self):
    #please return
    # miller_array of niggli stuff
    # Preffered symmetry
    # Highest symmetry
    result = [self.miller_niggli.f_as_f_sq(), self.pg_low_prim_set_name, self.pg_choice,self.pg_lattice_name]
    return result


  def show(self, out=None):
    if out == None:
      out = sys.stdout
    print >> out
    print >> out
    print >> out, "Exploring higher metric symmetry"
    if self.sigma_warning is not None:
      print >> out
      print >> out, self.sigma_warning

    self.caption1 = """
The point group of data as dictated by the space group is %s
The point group in the niggli setting is %s
The point group of the lattice is %s
A summary of R values for various possible point groups follow.
""" % (str(self.pg_input_name),
       str(sgtbx.space_group_info(group=self.pg_low_prim_set)),
       str(self.pg_lattice_name))
    print >> out, self.caption1
    print >> out, self.table
    self.caption2 = """
R_used: mean and maximum R value for symmetry operators *used* in this point group
R_unused: mean and minimum R value for symmetry operators *not used* in this point group
An automated point group suggestion is made on the basis of the BIC (Bayesian information criterion).

"""

    print >> out, self.caption2
    print >> out, "The likely point group of the data is:", self.pg_choice
    print >> out, ""
    print >> out, "Possible space groups in this point group are:"
    for sg in self.sg_possibilities:
      sg[0].show_summary(f=out, prefix= "   ")
      print >> out
    self.caption3 = """
Note that this analysis does not take into account the effects of twinning.
If the data are (almost) perfectly twinned, the symmetry will appear to be
higher than it actually is.
"""
    print >> out, self.caption3




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
      """
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
      """

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
               miller_calc=None,
               additional_parameters=None):

    self.max_delta = 3.0

    symm_issue_table = [0.08, 75, 0.08, 75]
    perform_ncs_analyses=False
    n_ncs_bins=7

    sigma_inflation = 1.25
    if additional_parameters is not None:
      sigma_inflation = additional_parameters.missing_symmetry.sigma_inflation
      perform_ncs_analyses = additional_parameters.twinning_with_ncs.perform_analyses
      n_ncs_bins = additional_parameters.twinning_with_ncs.n_bins

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

    ## Make sure that we actually have some data
    if miller_array.indices().size()==0:
      raise Sorry("No suitable data available after resolution cuts")

    ## Determine possible twin laws
    print >> out, "Determining possible twin laws."
    possible_twin_laws = twin_laws(miller_array,lattice_symmetry_max_delta=self.max_delta)
    possible_twin_laws.show(out=out)
    ##-----------------------------
    self.possible_twin_laws = possible_twin_laws
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

    self.abs_sg_anal = None
    try:
      if miller_array.sigmas() is not None:
        # Look at systematic absences please
        import absences
        self.abs_sg_anal = absences.protein_space_group_choices(
           miller_array = self.normalised_intensities.all,
           threshold = 3.0, out=out, sigma_inflation=sigma_inflation)
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
        y_label = 'P(Z<=z)',
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
      # XXX: for GUI
      self.nz_test_table = data_plots.table_data(
        title = "NZ test",
        column_labels=["z", "Acentric observed", "Acentric untwinned",
          "Centric observed", "Centric untwinned"],
        graph_names=["NZ test"],
        graph_labels=[("Z", "P(Z<=z)")],
        graph_columns=[list(range(5))],
        data=[self.nz_test.z, self.nz_test.ac_obs, self.nz_test.ac_untwinned,
              self.nz_test.c_obs, self.nz_test.c_untwinned])

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
      # XXX: for GUI
      self.l_test_table = data_plots.table_data(
        title="L test, acentric data",
        column_labels=["|l|", "Observed", "Acentric theory",
                       "Acentric theory, perfect twin"],
        graph_names=["L test"],
        graph_labels=[("|l|", "P(L>=l)")],
        graph_columns=[list(range(4))],
        data=[self.l_test.l_values, self.l_test.l_cumul,
              self.l_test.l_cumul_untwinned, self.l_test.l_cumul_perfect_twin])
      ##------------------------

    ##--------------------------






    self.n_twin_laws = len(possible_twin_laws.operators)

    self.twin_law_dependent_analyses = []


    self.twin_law_info = {}
    self.twin_law_names = []
    for ii in range(self.n_twin_laws):
      twin_law_name = possible_twin_laws.operators[ii].operator.r().as_hkl()
      self.twin_law_names.append(twin_law_name)
      print >> out
      print >> out,"---------------------------------------------"
      print >> out," Analysing possible twin law : ", twin_law_name
      print >> out,"---------------------------------------------"

      tmp_twin_law_stuff = twin_law_dependend_twin_tests(
        possible_twin_laws.operators[ii],
        miller_array,
        out=out,
        verbose=verbose,
        miller_calc=miller_calc,
        normalized_intensities=self.normalised_intensities.acentric,
        ncs_test=perform_ncs_analyses,
        n_ncs_bins=n_ncs_bins)

      self.twin_law_dependent_analyses.append( tmp_twin_law_stuff )
      (britton_table, britton_frac, h_test_table, h_frac,ml_table, ml_frac) = \
        (None, None, None, None, None, None)
      if tmp_twin_law_stuff.results_available: # emight phase in a availability test in here

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
        # XXX: for GUI
        britton_frac = tmp_twin_law_stuff.britton_test.estimated_alpha
        britton_table = data_plots.table_data(
          title = "Britton plot for twin law %s" % twin_law_name,
          column_labels=["alpha", "percentage negatives", "fit"],
          graph_names=["Britton plot"],
          graph_labels=[("alpha", "percentage negatives")],
          graph_columns=[[0,1,2]],
          data=[tmp_twin_law_stuff.britton_test.britton_alpha,
                tmp_twin_law_stuff.britton_test.britton_obs,
                tmp_twin_law_stuff.britton_test.britton_fit])
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
        # XXX: for GUI
        h_frac = tmp_twin_law_stuff.h_test.estimated_alpha
        h_test_table = data_plots.table_data(
          title="H test for possible twin law %s" % twin_law_name,
          column_labels=["H", "Observed S(H)", "Fitted S(H)"],
          graph_names=["H test for acentric data"],
          graph_labels=[("H", "S(H)")],
          graph_columns=[[0,1,2]],
          data=[tmp_twin_law_stuff.h_test.h_array,
                tmp_twin_law_stuff.h_test.cumul_obs,
                tmp_twin_law_stuff.h_test.cumul_fit])
        # plot the likelihood profile fro the ML rees test please
        if tmp_twin_law_stuff.ml_murray_rust is not None:
          ml_murray_rust_plot = data_plots.plot_data(
            plot_title = 'Likelihood based twin fraction estimation for possible twin law '\
                          +possible_twin_laws.operators[ii].operator.r().as_hkl(),
            x_label = 'alpha',
            y_label = '-Log[Likelihood]',
            x_data = tmp_twin_law_stuff.ml_murray_rust.twin_fraction,
            y_data = tmp_twin_law_stuff.ml_murray_rust.nll,
            y_legend = 'NLL (acentric data)',
            comments = 'Likelihood based twin fraction estimate')
          if out_plots is not None:
            data_plots.plot_data_loggraph(ml_murray_rust_plot,out_plots)
          # XXX: for GUI
          ml_frac = tmp_twin_law_stuff.ml_murray_rust.estimated_alpha
          ml_table = data_plots.table_data(
            title="Likelihood based twin fraction estimation for twin law %s"%\
                  twin_law_name,
            column_labels=["alpha", "-Log[Likelihood]"],
            graph_names=["Likelihood-based twin fraction estimate"],
            graph_labels=[("alpha", "-Log[Likelihood]")],
            graph_columns=[[0,1]],
            data=[tmp_twin_law_stuff.ml_murray_rust.twin_fraction,
                  tmp_twin_law_stuff.ml_murray_rust.nll])
      self.twin_law_info[twin_law_name] = (britton_table, britton_frac,
                                h_test_table, h_frac,ml_table, ml_frac)
    # now we can check for space group related issues
    self.check_sg = None
    self.suggested_space_group=None
    if self.n_twin_laws > 0:
      self.check_sg = symmetry_issues(
        miller_array,
        self.max_delta,
        out=out,
        sigma_inflation=sigma_inflation)

      nig_data, pg_this_one, pg_choice, pg_high = self.check_sg.return_point_groups()
      xs_choice = crystal.symmetry( nig_data.unit_cell(), pg_choice, assert_is_compatible_unit_cell=False )
      xs_high   = crystal.symmetry( nig_data.unit_cell(), pg_high, assert_is_compatible_unit_cell=False )

      if pg_choice != pg_high:
           merge_data_and_guess_space_groups(miller_array=nig_data, xs=xs_high,out=out,
                                             txt="Merging in *highest possible* point group %s.\n ***** THIS MIGHT NOT BE THE BEST POINT GROUP SYMMETRY *****  "%pg_high  )


      if pg_choice != pg_this_one:
        suggested_space_group = merge_data_and_guess_space_groups(miller_array=nig_data, xs=xs_choice,out=out,
                                                                       txt="Merging in *suggested* point group %s "%pg_choice  )

    ##--------------------------
    self.twin_summary = twin_results_interpretation(
      self.nz_test,
      self.wilson_moments,
      self.l_test,
      self.translation_pseudo_symmetry,
      self.twin_law_dependent_analyses,
      self.check_sg,
      out=out)



def merge_data_and_guess_space_groups(miller_array, txt, xs=None,out=None, sigma_inflation=1.0):
  tmp_ma = miller_array.deep_copy()
  if xs is None:
    xs = tmp_ma.crystal_symmetry()
  tmp_ma = miller_array.customized_copy( crystal_symmetry=xs )
  merge_obj = tmp_ma.change_basis( xs.space_group_info().change_of_basis_op_to_reference_setting() ).merge_equivalents()
  tmp_ma = merge_obj.array()
  r_lin = merge_obj.r_linear()
  normalizer = absolute_scaling.kernel_normalisation(tmp_ma, auto_kernel=True)
  work_array = normalizer.normalised_miller.deep_copy()
  abs_sg_anal = None
  if tmp_ma.sigmas() is not None:
    print >> out
    print >> out
    print >> out, "-"*len(txt)
    print >> out, txt
    print >> out, "-"*len(txt)
    print >> out
    merge_obj.show_summary(out=out)
    print >> out
    print >> out, "Suggesting various space group choices on the basis of systematic absence analyses"
    print >> out
    print >> out
    this_worked=False
    try:
      if miller_array.sigmas() is not None:
          # Look at systematic absences please
          import absences
          abs_sg_anal = absences.protein_space_group_choices(miller_array = work_array,
              threshold = 3.0, out=out, print_all=False, sigma_inflation=sigma_inflation)
          this_worked=True
    except Sorry: pass
    if not this_worked:
      print >> out, "Systematic absence analyses failed"
  return (merge_obj, abs_sg_anal)




def twin_analyses_brief(miller_array,
                        cut_off=4.0,
                        out = None,
                        verbose=0):
  """
  A very brief twin analyses and tries to answer the question whether or
  not the data are twinned.
  possible outputs and the meaning:
  - False: data are not twinned
  - True : data do not behave as expected. One possible explanantion
           is twinning
  - None : data do not behave as expected, and might or might not be
           due to twinning.
           Also gives none when something messes up.
  """
  if(out is None):
    out = sys.stdout
  # first we need to know wheter or not that sigmas make any sense at all
  if (not miller_array.sigmas_are_sensible()):
    #clearly there is something wrong with the sigmas
    #forget about them I would say
    miller_array = miller_array.customized_copy( indices=miller_array.indices(),
                                                 data=miller_array.data(),
                                                 sigmas=None ).set_observation_type( miller_array )

  out_tmp = StringIO()
  out_tmp_plot = StringIO()
  twin_results = None
  twinned=None

  if not miller_array.space_group().is_centric():
    try:
      twin_results = twin_analyses(miller_array,
                                   d_star_sq_low_limit=1.0/100.0,
                                   d_star_sq_high_limit=1.0/(0.001**2.0),
                                   out = out_tmp,
                                   out_plots = out_tmp_plot,
                                   verbose=verbose)
    except Exception, e:
      print >> out, "Twin analysis failed:", str(e)
      return None

  if twin_results is not None:
    if verbose>0:
      twin_results.twin_summary.twin_results.show(out=out)
    if (twin_results.twin_summary.twin_results.maha_l>cut_off):
      if twin_results.twin_summary.twin_results.l_mean <= 0.48:
        twinned = True
    if (twin_results.twin_summary.twin_results.maha_l<=cut_off):
        twinned = False
  return twinned
