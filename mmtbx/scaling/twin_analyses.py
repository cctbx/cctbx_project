
from __future__ import absolute_import, division, print_function
from mmtbx.scaling import absolute_scaling
from mmtbx import scaling
from iotbx import data_plots
from cctbx.sgtbx import pointgroup_tools
from cctbx.array_family import flex
import cctbx.sgtbx.lattice_symmetry
from cctbx import crystal
import cctbx.sgtbx.cosets
from cctbx import maptbx
from cctbx import sgtbx
import cctbx.xray
import scitbx.math
from scitbx import matrix
from libtbx import slots_getstate_setstate, Auto
from libtbx.str_utils import format_value
from libtbx.utils import Sorry, null_out
from libtbx import table_utils
from six.moves import cStringIO as StringIO
import math
import sys
from six.moves import zip
from six.moves import range

# some cutoffs that may need to be adjusted
TWIN_FRAC_SIGNIFICANT = 0.05
PATT_HEIGHT_MISSED_CENTERING = 75
PATT_HEIGHT_TNCS = 0.2

class obliquity(slots_getstate_setstate):
  __slots__ = ["type", "u", "h", "t", "tau", "delta"]
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
  """
  Various scores for a potential twin law given the crystal lattice.
  """
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
    """
    this gives a tensor describing the deformation of the unit cell needed
    to obtain a perfect match. the sum of diagonal elements describes
    the change in volume, off diagonal components measure associated shear.
    """
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
class twin_law(slots_getstate_setstate):
  """
  Basic container for information about a possible twin law, with scores for
  fit to crystal lattice.
  """
  __slots__ = ["operator", "twin_type", "delta_santoro", "delta_le_page",
               "delta_lebedev", "axis_type"]
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

  def __str__(self):
    return str(self.operator.r().as_hkl())

class twin_laws(scaling.xtriage_analysis):
  """
  Container for all possible twin laws given a crystal lattice and space
  group.
  """
  def __init__(self,
               miller_array,
               lattice_symmetry_max_delta=3.0,
               out=None):
    input = miller_array.eliminate_sys_absent(integral_only=True,
        log=out)
    self.change_of_basis_op_to_niggli_cell = \
      input.change_of_basis_op_to_niggli_cell()
    minimum_cell_symmetry = input.crystal_symmetry().change_basis(
      cb_op=self.change_of_basis_op_to_niggli_cell)
    self.lattice_group = sgtbx.lattice_symmetry.group(
      reduced_cell=minimum_cell_symmetry.unit_cell(),
      max_delta=lattice_symmetry_max_delta)
    intensity_symmetry = minimum_cell_symmetry.reflection_intensity_symmetry(
        anomalous_flag=input.anomalous_flag())
    self.euclid = intensity_symmetry.space_group_info().type()\
      .expand_addl_generators_of_euclidean_normalizer(flag_k2l=True,
                                                      flag_l2n=True )
    self.operators = []
    self.m = 0  # number of merohedral twin laws
    self.pm = 0 # number of pseudo-merohedral twin laws
    cb_op = self.change_of_basis_op_to_niggli_cell.inverse()
    for partition in sgtbx.cosets.left_decomposition(
        g = self.lattice_group,
        h = intensity_symmetry.space_group()
            .build_derived_acentric_group()
            .make_tidy()).partitions[1:] :
      if (partition[0].r().determinant() > 0):
        is_pseudo_merohedral=False
        twin_type =str("  M")
        self.m+=1
        euclid_check = sgtbx.space_group( self.euclid )
        try:
          euclid_check.expand_smx( partition[0] )
        except KeyboardInterrupt: raise
        except Exception:
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
        tlq = twin_law_quality(
          xs=miller_array,
          twin_law=cb_op.apply(partition[0]))
        tl = twin_law(
          op=cb_op.apply(partition[0]),
          pseudo_merohedral_flag=str(twin_type),
          axis_type=tlq.delta_le_page()[1],
          delta_santoro=tlq.delta_santoro(),
          delta_le_page=tlq.delta_le_page()[0],
          delta_lebedev=tlq.delta_lebedev())
        self.operators.append(tl)
    table_rows = []
    for tl in self.operators:
      table_rows.append([
        tl.twin_type,
        tl.axis_type,
        "%5.3f"%(tl.delta_santoro),
        "%5.3f"%(tl.delta_le_page),
        "%5.3f"%(tl.delta_lebedev),
        str(tl.operator.r().as_hkl())
      ])
    self.table = table_utils.simple_table(
      column_headers=['Type', 'Axis', 'R metric (%)', 'delta (le Page)',
                      'delta (Lebedev)', 'Twin law'],
      table_rows=table_rows,
      comments="""\
M:  Merohedral twin law
PM: Pseudomerohedral twin law""")
    self.coset_table = sgtbx.cosets.left_decomposition(self.lattice_group,
      intensity_symmetry.space_group().build_derived_acentric_group()\
            .make_tidy())
    self.coset_decomposition = False
    try:
      tmp = self.lattice_group.change_basis(
          self.change_of_basis_op_to_niggli_cell.inverse() )
      self.coset_decomposition = True
    except Exception:
      pass

  def _show_impl(self, out):
    assert (self.m + self.pm)==len(self.operators)
    out.show_sub_header("Twin law identification")
    if (len(self.operators) == 0):
      out.show("\nNo twin laws are possible for this crystal lattice.\n")
    else:
      out.show_paragraph_header("Possible twin laws:")
      out.show_table(self.table, indent=2)
      out.show_lines("""
%(m)-3.0f merohedral twin operators found
%(pm)-3.0f pseudo-merohedral twin operators found
In total, %(n_op)3.0f twin operators were found""" %
        ({ "m" : self.m, "pm" : self.pm, "n_op" : len(self.operators) }))
      out.show("""
 Please note that the possibility of twin laws only means that the lattice
 symmetry permits twinning; it does not mean that the data are actually
 twinned.  You should only treat the data as twinned if the intensity
 statistics are abnormal.""")
    # XXX disabled 2014-06-28
    if False :
      out.show("""
 The presence of twin laws indicates that the symmetry of the lattice (unit
 cell) is higher (has more elements) than the point group of the assigned
 space group.""")
      out.show("""
 There are four likely scenarios associated with the presence of twin laws:""")
      out.show_lines("""
    i.  The assigned space group is incorrect (too low).
   ii.  The assigned space group is correct and the data *are not* twinned.
  iii.  The assigned space group is correct and the data *are* twinned.
   iv.  The assigned space group is not correct (too low) and at the same time,
        the data *are* twinned.""")
      out.show("""
 Xtriage tries to distinguish between these cases by inspecting the intensity
 statistics.  It never hurts to carefully inspect statistics yourself and make
 sure that the automated interpretation is correct.""")

      if self.coset_decomposition :
        out.show_sub_header("Details of automated twin law derivation")
        out.show("""
Below, the results of the coset decomposition are given.  Each coset represents
a single twin law, and all symmetry equivalent twin laws are given.  For each
coset, the operator in (x,y,z) and (h,k,l) notation are given.   The direction
of the axis (in fractional coordinates), the type and possible offsets are
given as well.  Furthermore, the result of combining a certain coset with the
input space group is listed.

This table can be usefull when comparing twin laws generated by xtriage with
those listed in lookup tables.  In the table subgroup H denotes the *presumed
intensity symmetry*. Group G is the symmetry of the lattice.
""")
        coset_out = StringIO()
        self.coset_table.show(out=coset_out,
          cb_op=self.change_of_basis_op_to_niggli_cell.inverse())
        out.show_preformatted_text(coset_out.getvalue())
        out.show("""\
Note that if group H is centered (C,P,I,F), elements corresponding to centering
operators are omitted.  (This is because internally the calculations are done
with the symmetry of the reduced cell.)
""")
      else: # lattice_group.change_basis() failed
        out.show("""
The present crystal symmetry does not allow to have the its lattice symmetry
expressed in the setting desired.  Becasue of this, a full coset table cannot
be produced. Working with the data in the reduced cell will solve this.
""")

def get_twin_laws(miller_array):
  """
  Convenience method for getting a list of twin law operators (as strings)
  """
  return [ str(tl) for tl in twin_laws(miller_array).operators ]

class wilson_normalised_intensities(scaling.xtriage_analysis):
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

  def _show_impl(self, out):
    out.show("Splitting data in centrics and acentrics")
    out.show("  Number of centrics  : %d" % self.centric.data().size())
    out.show("  Number of acentrics : %d", self.acentric.data().size())

########################################################################
# DIAGNOSTIC TESTS
########################################################################

class detect_pseudo_translations(scaling.xtriage_analysis):
  """
  Analyze the Patterson map to identify off-origin peaks that are a significant
  fraction of the origin peak height.
  """
  def __init__(self,
      miller_array, # ideally amplitudes
      low_limit=10.0,
      high_limit=5.0,
      max_sites=100,
      height_cut=0.0,
      distance_cut=15.0,
      p_value_cut=0.05,
      completeness_cut=0.75,
      cut_radius=3.5,
      min_cubicle_edge=5.0,
      completeness_as_non_anomalous=None,
      out=None,
      verbose=0):
    if out is None:
      out=sys.stdout
    if miller_array.is_xray_intensity_array():
      miller_array = miller_array.f_sq_as_f()
    work_array = miller_array.resolution_filter(low_limit,high_limit)
    work_array = work_array.select(work_array.data()>0).set_observation_type(
      miller_array)
    self.space_group = miller_array.space_group()
    self.unit_cell = miller_array.unit_cell()
    self.d_max = low_limit
    self.d_min = high_limit
    self.completeness_in_range = work_array.completeness(
      as_non_anomalous_array = completeness_as_non_anomalous)
    self.xs = crystal.symmetry(unit_cell=miller_array.unit_cell(),
                               space_group=miller_array.space_group())

    if work_array.completeness(
      as_non_anomalous_array = completeness_as_non_anomalous) \
         <completeness_cut:
      print(file=out)
      print(" WARNING (twin_analysis):", file=out)
      print("  The completeness is only %3.2f between %3.1f and %3.1f A."% (
        work_array.completeness(
          as_non_anomalous_array = completeness_as_non_anomalous
          ), low_limit, high_limit), file=out)
      print("  This might not be enough to obtain a good estimate", file=out)
      print("  of the presence or absence of pseudo translational", file=out)
      print("  symmetry.", file=out)
    if work_array.indices().size()==0:
      raise Sorry("No low resolution reflections")

    if work_array.anomalous_flag():
      work_array = work_array.average_bijvoet_mates().set_observation_type(
        miller_array)
    everything_okay = True
    if (work_array.indices().size()<0):
      print("The number of reflection between %3.1f and %3.1f Angstrom" \
         %( low_limit,
            high_limit ), file=out)
      print("is equal to %i" %(work_array.indices().size()), file=out)
      print(" ##  This is not enough to obtain a reasonable estimate of", file=out)
      print(" ##  the presence of translational NCS", file=out)
      everything_okay = False

    if everything_okay:
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
        print(file=out)
        print("No Patterson vectors with a length larger than", file=out)
        print("%5.2f found. removing distance constraint"%(distance_cut), file=out)
        print(file=out)
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
      self.suggested_space_groups_table = None
      if everything_okay:
        self.high_peak = self.suspected_peaks[0][1]
        self.high_peak_distance = self.suspected_peaks[0][3]
        self.high_peak_xyz = self.suspected_peaks[0][0]
        self.high_p_value = self.suspected_peaks[0][2]
        if( self.high_p_value <= self.p_value_cut):
          self.guesstimate_mod_hkl()

  def suggest_new_space_groups(self,t_den=144,out=None):
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
      tmp_space_group = sgtbx.space_group_info( str(tmp_space_group),
        space_group_t_den=t_den)
    except KeyboardInterrupt: raise
    except Exception : pass

    for so in symops:
      sg_str = None
      try:
        new_tmp_sg = tmp_space_group.group().make_tidy()
        smx = sgtbx.rt_mx( so )
        smx = smx.new_denominators( new_tmp_sg.r_den(), new_tmp_sg.t_den() )
        new_tmp_sg.expand_smx( smx )
        sg_str = str( sgtbx.space_group_info( group = new_tmp_sg,
          space_group_t_den=t_den  ) )
        to_ref_set = sgtbx.space_group_info(
          group=new_tmp_sg,
          space_group_t_den=t_den).change_of_basis_op_to_reference_setting()
      except Exception: pass
      if sg_str not in sgs:
        if sg_str is not None:
          sgs.append( sg_str )
          with_operator.append( so )
          new_cell.append( self.unit_cell.change_basis(
            to_ref_set.c_inv().r().as_double() ) )
        else:
          sgs.append( None )
          new_cell.append( self.unit_cell )
    number_of_new_sgs = len(symops)-sgs.count(None)
    if number_of_new_sgs>0:
      table_rows = []
      for sg,op,uc in zip(sgs,with_operator,new_cell):
        if sg is not None:
          table_rows.append([str(sg), str(op),
            "%5.2f, %5.2f, %5.2f,  %5.2f, %5.2f, %5.2f)"% uc.parameters() ])
      self.suggested_space_groups_table = table_utils.simple_table(
        table_rows=table_rows,
        column_headers=("Space group", "Operator",
          "Unit cell of reference setting"))

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
      for trial_num in range(start_num,trial_den):
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

  def _show_impl(self, out):
    out.show_sub_header("Patterson analyses")
    out.show(" Largest Patterson peak with length larger than 15 Angstrom:")
    out.show_preformatted_text(" Frac. coord.              : %8.3f %8.3f %8.3f"
      % self.high_peak_xyz)
    out.show_preformatted_text("""\
 Distance to origin        : %8.3f
 Height relative to origin : %8.3f %%
 p_value(height)           : %12.3e
""" % (self.high_peak_distance,self.high_peak, self.high_p_value))
    out.show_paragraph_header("Explanation")
    out.show("""\
 The p-value, the probability that a peak of the specified height or larger
 is found in a Patterson function of a macromolecule that does not have any
 translational pseudo-symmetry, is equal to %.3e.  p_values smaller than
 0.05 might indicate weak translational pseudo symmetry, or the self vector of
 a large anomalous scatterer such as Hg, whereas values smaller than 1e-3 are
 a very strong indication for the presence of translational pseudo symmetry.
""" % self.high_p_value)
    if (self.high_peak >= 20):
      out.show_text("""\
 Translational pseudo-symmetry is very likely present in these data.  Be
 aware that this will change the intensity statistics and may impact subsequent
 analyses, and in practice may lead to higher R-factors in refinement.
""")
    if (self.high_p_value <= self.p_value_cut):
      table_rows = []
      for ii in range(len(self.suspected_peaks)):
        if self.suspected_peaks[ii][2] > self.p_value_cut:
          break
        table_rows.append([
          "%6.3f,%6.3f,%6.3f" % self.suspected_peaks[ii][0],
          "%8.3f" % self.suspected_peaks[ii][1],
          "%9.3e" % self.suspected_peaks[ii][2] ])
      if (len(table_rows) > 1):
        out.show(" The full list of Patterson peaks is:")
        table = table_utils.simple_table(
          table_rows=table_rows,
          column_headers=["XYZ", "height", "p-value(height)"])
        out.show_table(table, indent=2)
      if (self.suggested_space_groups_table is not None):
        out.show("""\
 If the observed pseudo translationals are crystallographic, the following
 spacegroups and unit cells are possible:
""")
        out.show_table(self.suggested_space_groups_table, indent=2)

class wilson_moments(scaling.xtriage_analysis):
  centric_i_ratio_library= [3.0,2.0]
  centric_f_ratio_library = [0.637,0.785]
  centric_e_sq_minus_one_library = [0.968,0.736]
  acentric_i_ratio_library= [2.0,1.5]
  acentric_f_ratio_library = [0.785,0.885]
  acentric_e_sq_minus_one_library = [0.736, 0.541]
  def __init__(self,
               acentric_z,
               centric_z):
    self.centric_i_ratio = None
    self.centric_f_ratio = None
    self.centric_e_sq_minus_one = None
    self.acentric_i_ratio = None
    self.acentric_f_ratio = None
    self.acentric_abs_e_sq_minus_one = None
    self.compute_ratios(
      acentric_z.data()/acentric_z.epsilons().data().as_double(),
      centric_z.data()/centric_z.epsilons().data().as_double())
    self.centric_present = True
    if centric_z.data().size()==0:
      self.centric_present=False

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

  def _show_impl(self, out):
    out.show_sub_header("Wilson ratio and moments")
    out.show("""Acentric reflections:\n""")
    out.show_preformatted_text("""
   <I^2>/<I>^2    :%4.3f   (untwinned: %4.3f; perfect twin %4.3f)
   <F>^2/<F^2>    :%4.3f   (untwinned: %4.3f; perfect twin %4.3f)
   <|E^2 - 1|>    :%4.3f   (untwinned: %4.3f; perfect twin %4.3f)
""" % (self.acentric_i_ratio,
       self.acentric_i_ratio_library[0],
       self.acentric_i_ratio_library[1],
       self.acentric_f_ratio,
       self.acentric_f_ratio_library[0],
       self.acentric_f_ratio_library[1],
       self.acentric_abs_e_sq_minus_one,
       self.acentric_e_sq_minus_one_library[0],
       self.acentric_e_sq_minus_one_library[1]))
    if self.centric_present:
      out.show("""Centric reflections:\n""")
      out.show_preformatted_text("""
   <I^2>/<I>^2    :%4.3f   (untwinned: %4.3f; perfect twin %4.3f)
   <F>^2/<F^2>    :%4.3f   (untwinned: %4.3f; perfect twin %4.3f)
   <|E^2 - 1|>    :%4.3f   (untwinned: %4.3f; perfect twin %4.3f)
""" % (self.centric_i_ratio,
       self.centric_i_ratio_library[0],
       self.centric_i_ratio_library[1],
       self.centric_f_ratio,
       self.centric_f_ratio_library[0],
       self.centric_f_ratio_library[1],
       self.centric_abs_e_sq_minus_one,
       self.centric_e_sq_minus_one_library[0],
       self.centric_e_sq_minus_one_library[1]))
    # TODO explanation, citation?

class n_z_test(scaling.xtriage_analysis):
  def __init__(self,
               normalised_acentric,
               normalised_centric):
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

  @property
  def table(self):
    return data_plots.table_data(
      title = "NZ test",
      column_labels=["z", "Acentric observed", "Acentric untwinned",
        "Centric observed", "Centric untwinned"],
      graph_names=["NZ test"],
      graph_labels=[("Z", "P(Z<=z)")],
      graph_columns=[list(range(5))],
      # make sure these aren't flex arrays - JSON can't handle them
      data=[list(self.z), list(self.ac_obs), list(self.ac_untwinned),
            list(self.c_obs), list(self.c_untwinned) ])

  def _show_impl(self, out):
    out.show_sub_header("NZ test for twinning and TNCS")
    out.show("""
The NZ test is diagnostic for both twinning and translational NCS.  Note
however that if both are present, the effects may cancel each other out,
therefore the results of the Patterson analysis and L-test also need to be
considered.
""")
    self.sign_ac = '+'
    if self.mean_diff_ac < 0:
      self.sign_ac = '-'
    self.sign_c = '+'
    if self.mean_diff_c < 0:
      self.sign_c = '-'
    out.show_preformatted_text("""
  Maximum deviation acentric      :  %(max_diff_ac)4.3f
  Maximum deviation centric       :  %(max_diff_c)4.3f

  <NZ(obs)-NZ(twinned)>_acentric  : %(sign_ac)1s%(mean_diff_ac)4.3f
  <NZ(obs)-NZ(twinned)>_centric   : %(sign_c)1s%(mean_diff_c)4.3f
""" % { "max_diff_ac" : self.max_diff_ac,
        "max_diff_c" : self.max_diff_c,
        "sign_ac" : self.sign_ac,
        "sign_c" : self.sign_c,
        "mean_diff_ac" : math.fabs(self.mean_diff_ac),
        "mean_diff_c" : math.fabs(self.mean_diff_c), })
    if (out.gui_output):
      out.show_plot(self.table)
      out.show_table(self.table)
    else :
      out.show_table(self.table)

class l_test(scaling.xtriage_analysis):
  """
  Implementation of:

  J. Padilla & T. O. Yeates. A statistic for local intensity differences:
  robustness to anisotropy and pseudo-centering and utility for detecting
  twinning. Acta Crystallogr. D59, 1124-30, 2003.

  This is complementary to the NZ test, but is insensitive to translational
  pseuo-symmetry.
  """
  def __init__(self,
      miller_array,
      parity_h=2,
      parity_k=2,
      parity_l=2):
    acentric_data = miller_array.select_acentric().set_observation_type(
      miller_array)
    if not miller_array.is_xray_intensity_array():
      acentric_data = acentric_data.f_as_f_sq()
    self.parity_h = parity_h
    self.parity_k = parity_k
    self.parity_l = parity_l
    l_stats = scaling.l_test(
      miller_indices=acentric_data.indices(),
      intensity=acentric_data.data()/\
        acentric_data.epsilons().data().as_double(),
      space_group=acentric_data.space_group(),
      anomalous_flag=acentric_data.anomalous_flag(),
      parity_h=parity_h,
      parity_k=parity_k,
      parity_l=parity_l,
      max_delta_h=8);
    self.mean_l = l_stats.mean_l()
    self.mean_l2 = l_stats.mean_l2()
    self.l_cumul = l_stats.cumul()
    self.l_values = flex.double(range(self.l_cumul.size()))/float(
      self.l_cumul.size())
    self.l_cumul_untwinned = self.l_values
    self.l_cumul_perfect_twin = self.l_values*(
      3.0-self.l_values*self.l_values)/2.0
    self.ml_alpha = l_stats.ml_alpha()

  @property
  def table(self):
    return data_plots.table_data(
      title="L test, acentric data",
      column_labels=["|l|", "Observed", "Acentric theory",
                     "Acentric theory, perfect twin"],
      graph_names=["L test"],
      graph_labels=[("|l|", "P(L>=l)")],
      graph_columns=[list(range(4))],
      data=[list(self.l_values), list(self.l_cumul),
            list(self.l_cumul_untwinned), list(self.l_cumul_perfect_twin)])

  def _show_impl(self, out):
    out.show_sub_header("L test for acentric data")
    out.show("""Using difference vectors (dh,dk,dl) of the form:""")
    out.show_preformatted_text(
      """    (%(parity_h)ihp, %(parity_k)ikp, %(parity_l)ilp)""" %
      { "parity_h" : self.parity_h,
        "parity_k" : self.parity_k,
        "parity_l" : self.parity_l })
    out.show("""where hp, kp, and lp are random signed integers such that""")
    out.show_preformatted_text("""    2 <= |dh| + |dk| + |dl| <= 8""")
    out.show_lines("""\
  Mean |L|   :%(mean_l)4.3f  (untwinned: 0.500; perfect twin: 0.375)
  Mean  L^2  :%(mean_l2)4.3f  (untwinned: 0.333; perfect twin: 0.200)""" %
      { "mean_l" : self.mean_l,
        "mean_l2" : self.mean_l2, })
    out.show("""
 The distribution of |L| values indicates a twin fraction of
 %(ml_alpha)3.2f. Note that this estimate is not as reliable as obtained
 via a Britton plot or H-test if twin laws are available.
""" % { "ml_alpha" : self.ml_alpha, })
    if (out.gui_output):
      out.show_plot(self.table)
      out.show_table(self.table)
    else :
      out.show_table(self.table)

    out.show("""\
 Reference:
  J. Padilla & T. O. Yeates. A statistic for local intensity differences:
  robustness to anisotropy and pseudo-centering and utility for detecting
  twinning. Acta Crystallogr. D59, 1124-30, 2003.
""")

########################################################################
# TWIN LAW-DEPENDENT TESTS
########################################################################

class britton_test(scaling.xtriage_analysis):
  def __init__(self,
               twin_law,
               miller_array, # ideally intensities!
               cc_cut_off=0.995,
               verbose=0):
    self.twin_law = twin_law
    self.max_iter=1000
    result = [0.5,1.0,0,0]
    miller_array = miller_array.select(
      miller_array.data() > 0).set_observation_type(miller_array)
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

  @property
  def table(self):
    return data_plots.table_data(
        title = "Britton plot for twin law %s" % str(self.twin_law),
        column_labels=["alpha", "percentage negatives", "fit"],
        graph_names=["Britton plot"],
        graph_labels=[("alpha", "percentage negatives")],
        graph_columns=[[0,1,2]],
        reference_marks=[[self.estimated_alpha]],
        data=[self.britton_alpha, self.britton_obs, self.britton_fit])

  def _show_impl(self, out):
    out.show_paragraph_header("Britton analyses")
    out.show("""
  Extrapolation performed on  %(alpha_cut)3.2f < alpha < 0.495
  Estimated twin fraction: %(estimated_alpha)4.3f
  Correlation: %(correlation)5.4f
""" % { "alpha_cut" : self.alpha_cut,
        "estimated_alpha" : self.estimated_alpha,
        "correlation" : self.correlation })


class h_test(scaling.xtriage_analysis):
  def __init__(self,
               twin_law,
               miller_array, # ideally intensities!
               fraction=0.50):
    self.fraction = fraction
    self.twin_law = twin_law
    miller_array = miller_array.select(
      miller_array.data()>0).set_observation_type(miller_array)
    if not miller_array.is_xray_intensity_array():
      miller_array = miller_array.f_as_f_sq()
    if miller_array.is_real_array():
      acentric_data =  miller_array.select_acentric().set_observation_type(
        miller_array)
      try :
        h_test_object  = scaling.h_test(
          miller_indices=acentric_data.indices(),
          intensity=acentric_data.data(),
          sigma=acentric_data.sigmas(),
          space_group=acentric_data.space_group(),
          anomalous_flag=acentric_data.anomalous_flag(),
          twin_law=twin_law,
          fraction=fraction)
      except ValueError:
        if miller_array.completeness() < 0.05 : # XXX could check for anomalous
          raise Sorry("These data are severely incomplete, which breaks the "+
            "H-test for twinning.  We recommend that you use a full data set "+
            "in Xtriage, otherwise the statistical analyses may be invalid.")
        else :
          raise
      self.mean_h = h_test_object.mean_h()
      self.mean_h2 = h_test_object.mean_h2()
      self.estimated_alpha = h_test_object.alpha()
      self.alpha_from_mean_h = (self.mean_h*2.0-1.0)/-2.0
      self.h_array = h_test_object.h_array()
      self.h_values = h_test_object.h_values()
      self.cumul_obs = h_test_object.h_cumul_obs()
      self.cumul_fit = h_test_object.h_cumul_fit()

  @property
  def table(self):
    return data_plots.table_data(
        title="H test for possible twin law %s" % str(self.twin_law),
        column_labels=["H", "Observed S(H)", "Fitted S(H)"],
        graph_names=["H test for acentric data"],
        graph_labels=[("H", "S(H)")],
        graph_columns=[[0,1,2]],
        reference_marks=[[self.estimated_alpha]],
        data=[self.h_array, self.cumul_obs, self.cumul_fit])

  def _show_impl(self, out):
    out.show_paragraph_header("""H-test on acentric data""")
    out.show_lines("""\
Only %(fraction)3.1f %% of the strongest twin pairs were used.

  mean |H| : %(mean_h)4.3f  (0.50: untwinned; 0.0: 50%% twinned)
  mean H^2 : %(mean_h2)4.3f  (0.33: untwinned; 0.0: 50%% twinned)

Estimation of twin fraction via mean |H|: %(alpha_from_mean_h)4.3f
Estimation of twin fraction via cum. dist. of H: %(estimated_alpha)4.3f
""" % { "fraction" : self.fraction*100.0,
        "mean_h" : self.mean_h,
        "mean_h2" : self.mean_h2,
        "alpha_from_mean_h" : self.alpha_from_mean_h,
        "estimated_alpha" : self.estimated_alpha })

  def __getstate__(self):
    """
    Pickling optimization - see note below
    """
    self.h_values = None
    return self.__dict__

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
    non_zero_sigma_sel = miller_array.sigmas() > 0
    tmp_miller_array = miller_array.deep_copy().select(
      non_zero_sigma_sel).set_observation_type(miller_array)
    # make a binner
    self.binner = tmp_miller_array.setup_binner( n_bins=n_bins )
    self.out = out
    self.f = None
    self.cycle = 0
    self.ml_object = scaling.ml_twin_with_ncs(
      tmp_miller_array.data(),
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

    for ii in range(n_bins):
      self.x.append( -1.0 - ii/20.0 )
    self._show_info(out=out)
    term_parameters = scitbx.lbfgs.termination_parameters(max_iterations=1000)
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self,
      termination_params=term_parameters)
    scitbx.lbfgs.run(target_evaluator=self)
    self.log_likelihood = self.f
    self.twin_fraction = self.twin_cap/(1+math.exp(-(self.x[0])))
    table = []
    for ii in range(self.x.size()-1):
      d_ncs = self.x[ii+1]
      d_ncs = self.d_cap/( 1.0+math.exp(-d_ncs) )
      low = self.binner.bin_d_range(ii+1)[0]
      high = self.binner.bin_d_range(ii+1)[1]
      table.append(["%5.4f - %5.4f" % (low, high), "%4.3f" % d_ncs])
    self.table = table_utils.simple_table(
      table_rows=table,
      column_headers=["Resolution", "D_ncs"])
    self._show_results(out)
    if calc_data is not None:
      self.calc_correlation(
        obs=tmp_miller_array,
        calc=calc_data,
        out=out)

  def string_it(self,x):
    return str("%4.3f"%(x))

  def calc_correlation(self,obs,calc, out=sys.stdout):
    if calc.is_xray_amplitude_array():
      calc = calc.f_as_f_sq()
    calc = calc.common_set( other=obs )
    calc.use_binning_of( obs  )
    print("""
 The correlation of the calculated F^2 should be similar to the estimated
 values.

 Observed correlation between twin related, untwinned calculated F^2
 in resolutiuon ranges, as well as ewstimates D_ncs^2 values:
""", file=out)
    # now loop over all resolution bins, get twin related intensities and get
    # twin related intensities please
    print(" Bin    d_max     d_min     CC_obs   D_ncs^2 ", file=out)
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
      print("%3i)    %5.4f -- %5.4f :: %6s   %6s" % (i_bin, low,
        high, self.string_it(d_theory_sq), self.string_it(d_ncs_sq)), file=out)

  def compute_functional_and_gradients(self):
    tmp = self.ml_object.p_tot_given_t_and_coeff( self.x[0],
                                                  self.x[1:] )
    f = tmp[0]
    g = tmp[1:]

    self.f=f
    self.cycle+=1
    if self.cycle%5==0:
      print("%3i "%(self.cycle), end=' ', file=self.out)
      #self.print_it()
    else:
      print(".", end=' ', file=self.out)
    if self.cycle%30==0:
      print(file=self.out)
    self.out.flush()
    return f,flex.double(g)

  def _show_info(self, out):
    if (not isinstance(out, scaling.xtriage_output)):
      out = scaling.printed_output(out)
    out.show_paragraph_header("Maximum Likelihood twin fraction determination")
    out.show("""
 Estimation of twin fraction, while taking into account the effects of
 possible NCS parallel to the twin axis.
     (Zwart, Read, Grosse-Kunstleve & Adams, to be published.)

 A parameters D_ncs will be estimated as a function of resolution, together
 with a global twin fraction.  D_ncs is an estimate of the correlation
 coefficient between untwinned, error-free, twin related, normalized
 intensities.  Large values (0.95) could indicate an incorrect point group.
 Value of D_ncs larger than say, 0.5, could indicate the presence of NCS. The
 twin fraction should be smaller or similar to other estimates given elsewhere.

 The refinement can take some time.  For numerical stability issues, D_ncs is
 limited between 0 and 0.95.  The twin fraction is allowed to vary between 0
 and 0.45.  Refinement cycle numbers are printed out to keep you entertained.
""")

  def _show_results(self, out):
    if (not isinstance(out, scaling.xtriage_output)):
      out = scaling.printed_output(out)
    out.show("""
  Cycle : %(cycle)3i
  Log[likelihood]: %(log_likelihood)15.3f
  twin fraction: %(twin_fraction)4.3f
  D_ncs in resolution ranges:
"""%{"cycle":self.cycle,"log_likelihood":self.log_likelihood,"twin_fraction":self.twin_fraction})
    out.show_table(self.table, indent=2)

  def _show_impl(self, out):
    self._show_info(out)
    self._show_results(out)

class ml_murray_rust(scaling.xtriage_analysis):
  """
  Maximum-likelihood twin fraction estimation (Zwart, Read, Grosse-Kunstleve &
  Adams, to be published).
  """
  def __init__(self,
               miller_array, # must be intensities
               twin_law,
               n_points = 4):
    assert miller_array.is_xray_intensity_array()
    assert miller_array.sigmas() is not None
    non_zero_sigma_sel = miller_array.sigmas() > 0
    tmp_array = miller_array.select(non_zero_sigma_sel)
    ml_murray_rust_object = scaling.ml_murray_rust(
      z=tmp_array.data(),
      sig_z=tmp_array.sigmas(),
      indices=tmp_array.indices(),
      space_group=tmp_array.space_group(),
      anomalous_flag=tmp_array.anomalous_flag(),
      twin_law=twin_law,
      n_hermite=n_points)
    self.twin_law = twin_law
    self.twin_fraction = []
    self.nll = []
    for ii in range(1,23):
      t=ii/46.0
      self.twin_fraction.append( t )
      self.nll.append( - ml_murray_rust_object.fast_log_p_given_t( t )  )
    i_t_max = flex.min_index( flex.double( self.nll ) )
    self.estimated_alpha = None
    if (i_t_max >= 1) and (i_t_max < len(self.nll)-1 ):
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
      if ((self.estimated_alpha < tmp_t[0]) or
          (self.estimated_alpha > tmp_t[2])):
        self.estimated_alpha = self.twin_fraction[ i_t_max ]
    else:
      self.estimated_alpha = self.twin_fraction[ i_t_max ]

  @property
  def table(self):
    return data_plots.table_data(
      title="Likelihood based twin fraction estimation for twin law %s"%\
            str(self.twin_law),
      column_labels=["alpha", "-Log[Likelihood]"],
      graph_names=["Likelihood-based twin fraction estimate"],
      graph_labels=[("alpha", "-Log[Likelihood]")],
      graph_columns=[[0,1]],
      reference_marks=[[self.estimated_alpha]],
      data=[self.twin_fraction, self.nll])

  def _show_impl(self, out):
    out.show_paragraph_header("Maximum Likelihood twin fraction determination")
    out.show("""\
  Zwart, Read, Grosse-Kunstleve & Adams, to be published.
  The estimated twin fraction is equal to %4.3f
"""%(self.estimated_alpha))

def weighted_cc(x,y,w):
  """
  Utility function for correlation_analyses class.
  """
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

class correlation_analyses(scaling.xtriage_analysis):
  def __init__(self,
               miller_obs,
               miller_calc,
               twin_law,
               d_weight=0.1):
    self.twin_law=sgtbx.change_of_basis_op( twin_law )
    self.twin_law= self.twin_law.new_denominators(
      r_den=miller_calc.space_group().r_den(),
      t_den=miller_calc.space_group().t_den() )
    obs=miller_obs.deep_copy()
    calc=miller_calc.deep_copy()
    # the incomming data is normalized, normalize the calculated data as well please
    normalizer = absolute_scaling.kernel_normalisation(calc, auto_kernel=True)
    calc = normalizer.normalised_miller.deep_copy()
    calc = calc.select(calc.data()>0)
    calc_part2 = calc.change_basis( self.twin_law )
    calc_part2 = calc.customized_copy(
      indices = calc_part2.indices(),
      data = calc_part2.data()).map_to_asu()
    # make sure we have common sets of everything
    calc, calc_part2 = calc.common_sets( calc_part2 )
    obs,  calc = obs.common_sets( calc )
    obs,  calc_part2 = obs.common_sets( calc_part2 )
    self.d = d_weight
    self.alpha = []
    self.cc = []
    def twin_the_data_and_compute_cc(alpha):
      wx = obs.sigmas()
      x = obs.data()
      dd = obs.d_star_sq().data()
      dd = flex.exp( -7.7*self.d*self.d*dd )
      wy = 1.0-dd*dd
      y = (1-alpha)*calc.data() + alpha*calc_part2.data()
      y = dd*dd*y
      w_tot = wx*wx + wy*wy
      cc = weighted_cc( x,y,w_tot )
      return(cc)
    if (obs.data().size() > 0):
      for ii in range(50):
        alpha=ii/100.0
        self.alpha.append( alpha )
        cc = twin_the_data_and_compute_cc(alpha)
        self.cc.append(cc)
    self.max_alpha, self.max_cc = self.find_maximum()
    del self.d
    del self.alpha

  def find_maximum(self):
    max_cc = flex.max( flex.double( self.cc ) )
    max_alpha = flex.max_index( flex.double( self.cc ) )
    max_alpha =  self.alpha[ max_alpha ]
    return max_alpha, max_cc

  def _show_impl(self, out):
    out.show_paragraph_header("Correlation analyses")
    out.show("""
  The supplied calculated data are normalized and artificially twinned;
  subsequently a correlation with the observed data is computed.
""")
    out.show("Correlation : %.3f" % self.max_cc)
    out.show("Estimated twin fraction : %.3f" % self.max_alpha)

#
# WRAPPER CLASS
#
class twin_law_dependent_twin_tests(scaling.xtriage_analysis):
  """Twin law dependent test results"""
  def __init__(self,
      twin_law,
      miller_array,
      out,
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
    tmp_detwin_object = scaling.detwin(
      miller_indices=acentric_data.indices(),
      intensity=acentric_data.data(),
      sigma=acentric_data.sigmas(),
      space_group=acentric_data.space_group(),
      anomalous_flag=acentric_data.anomalous_flag(),
      twin_law=twin_law.operator.as_double_array()[0:9])
    self.twin_completeness = tmp_detwin_object.completeness()
    self.results_available = False
    self.h_test = self.h_test_table = None
    self.britton_test = self.britton_test_table = None
    self.r_values=None
    self.correlation=None
    self.ml_murray_rust = self.ml_murray_rust_table = None
    self.ml_murry_rust_with_ncs=None
    if self.twin_completeness > 0:
      self.results_available = True
      self.h_test = h_test(
        twin_law=twin_law.operator.as_double_array()[0:9],
        miller_array = acentric_data)
      self.britton_test = britton_test(
        twin_law=twin_law.operator.as_double_array()[0:9],
        miller_array=acentric_data,
        verbose=verbose)
      self.r_values = r_values(
        miller_obs=miller_array,
        twin_law=twin_law.operator.as_double_array()[0:9],
        miller_calc=miller_calc)
      if miller_calc is not None:
        ## Magic number to control influence of hiugh resolutiuon reflections
        d_weight=0.25
        try:
          self.correlation = correlation_analyses(
            miller_obs=normalized_intensities,
            miller_calc=miller_calc,
            twin_law=twin_law.operator,
            d_weight=d_weight)
        except KeyboardInterrupt: raise
        except Exception: pass
      if normalized_intensities.sigmas() is not None:
        self.ml_murray_rust = ml_murray_rust(
          miller_array=normalized_intensities,
          twin_law=twin_law.operator.as_double_array()[0:9])
        if ncs_test:
          self.ml_murry_rust_with_ncs = ml_murray_rust_with_ncs(
            miller_array=normalized_intensities,
            twin_law=twin_law.operator.as_double_array()[0:9],
            out=out,
            n_bins=n_ncs_bins,
            calc_data=miller_calc,
            start_alpha=self.ml_murray_rust.estimated_alpha)

  @property
  def h_frac(self):
    if (self.h_test is not None):
      return self.h_test.estimated_alpha

  @property
  def britton_frac(self):
    if (self.britton_test is not None):
      return self.britton_test.estimated_alpha

  @property
  def ml_frac(Self):
    if (self.ml_murray_rust is not None):
      return self.ml_murray_rust.estimated_alpha

  def _show_impl(self, out):
    out.show_sub_header("Analysis of twin law %s" % self.twin_law)
    if (self.results_available):
      self.h_test.show(out)
      self.britton_test.show(out)
      self.r_values.show(out)
      if (self.correlation is not None):
        self.correlation.show(out)
      tables = [self.h_test.table, self.britton_test.table]
      if (self.ml_murray_rust is not None):
        tables.append(self.ml_murray_rust.table)
      out.show_plots_row(tables)
    else:
      out.show(("The twin completeness is only %3.2f. Twin law dependent "+
        "tests not performed.") % (self.twin_completeness))

########################################################################
# OTHER SYMMETRY PROBLEMS
########################################################################

class symmetry_issues(scaling.xtriage_analysis):
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
    self.pg_low_prim_set_name = str(sgtbx.space_group_info(
      group=self.explore_sg.pg_low_prim_set))
    self.pg_lattice_name = str( sgtbx.space_group_info(
      group = self.explore_sg.pg_high ) )
    self.miller_niggli =self.miller_array.change_basis(
      self.xs_input.change_of_basis_op_to_niggli_cell())

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
    # loop over all pg's
    for pg in self.pg_max_r_used_table:
      tmp = self.miller_niggli
      if tmp.is_xray_amplitude_array():
        tmp = tmp.f_as_f_sq()
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
      else :
        non_zero_sigmas_sel = sigmas != 0
        tmp = tmp.select(non_zero_sigmas_sel)
        sigmas = tmp.sigmas()

      merger = cctbx.xray.merger( tmp.indices(),
                                  tmp.data(),
                                  sigmas*self.inflate_sigma ,
                                  sgtbx.space_group_info(pg).group(),
                                  tmp.anomalous_flag(),
                                  tmp.unit_cell() )
      final_score = -merger.bic()
      self.pg_scores.update( {pg:final_score} )
    # create output table
    legend = ['Point group', 'mean R_used', 'max R_used', 'mean R_unused',
              'min R_unused', 'BIC', 'choice']
    table_data = []
    min_score = 2e+9
    def as_string(val, fmt):
      return format_value(fmt, val, replace_none_with="None")
    for pg in self.pg_scores:
      tmp = [pg,
             as_string(self.pg_r_used_table[pg], "%4.3f"),
             as_string(self.pg_max_r_used_table[pg],"%4.3f"),
             as_string(self.pg_r_unused_table[pg],"%4.3f"),
             as_string(self.pg_min_r_unused_table[pg],"%4.3f"),
             as_string(self.pg_scores[pg],"%6.3e"),
             "    " ]
      table_data.append( tmp )
      if self.pg_scores[ pg ] < min_score:
        min_score = self.pg_scores[ pg ]
        self.pg_choice = pg
    for row in self.table_data:
      if self.pg_choice == row[0]:
        row[ 6 ] = "<---"
    self.table = table_utils.simple_table(
      column_headers=legend,
      table_rows=table_data)
    # store the possible spacegroups and the change of basis ops
    # for the most likely pg
    for xs in self.explore_sg.pg_graph.graph.node_objects[
      self.pg_choice ].allowed_xtal_syms:
      self.sg_possibilities.append( xs )
    self.xs_with_pg_choice_in_standard_setting = crystal.symmetry(
      unit_cell   = self.miller_niggli.unit_cell(),
      space_group = self.pg_choice,
      assert_is_compatible_unit_cell=False).as_reference_setting()

  def make_r_table(self):
    tmp_buffer=StringIO()
    # please find all missing sym ops
    start = str(sgtbx.space_group_info(group=self.pg_low_prim_set))
    # FIXME, ordering keys/values changes depending on py2/3
    tmp_key = list(self.explore_sg.pg_graph.graph.edge_objects[ start ].keys())[0]
    tmp_edge = self.explore_sg.pg_graph.graph.edge_objects[ start ][tmp_key]
    for symop in tmp_edge.return_used():
      # get the r value please for this symop
      r_value = r_values(
        miller_obs=self.miller_niggli,
        twin_law=symop.as_double_array()[0:9]).show(out=tmp_buffer)
      top = r_value.r_abs_top_obs
      bottom = r_value.r_abs_bottom_obs
      self.ops_and_r_pairs.update( {symop.r().as_hkl():
                                    (top, bottom)} )
    for symop in tmp_edge.return_unused():
      r_value = r_values(
        miller_obs=self.miller_niggli,
        twin_law=symop.as_double_array()[0:9]).show(out=tmp_buffer)
      top = r_value.r_abs_top_obs
      bottom = r_value.r_abs_bottom_obs
      self.ops_and_r_pairs.update( {symop.r().as_hkl():
                                    (top, bottom)} )

  def make_pg_r_table(self):
    start = str(sgtbx.space_group_info(group=self.pg_low_prim_set))
    end = str(sgtbx.space_group_info(group=self.explore_sg.pg_high))

    for pg in self.explore_sg.pg_graph.graph.node_objects:
      tot_in, max_in, tot_out, min_out, r_in, r_out = \
        self.get_r_value_total(start,pg)
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
    tmp_keys = list(self.ops_and_r_pairs.keys())
    tmp_values = list(self.ops_and_r_pairs.values())
    tmp = self.ops_and_r_pairs.copy()
    for op in coset_ops:
      if op in tmp_keys:  #tmp.has_key( op ):
        # Work around for absence of pop in python 2.2
        iii = tmp_keys.index( op )
        a = tmp_values[ iii ]
        # now we have to pop these guys from the array
        tmp_tmp_keys = []
        tmp_tmp_values = []
        for ii in range( len(tmp_keys) ):
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

  def return_point_groups(self):
    #please return
    # miller_array of niggli stuff
    # Preffered symmetry
    # Highest symmetry
    miller_niggli = self.miller_niggli
    if miller_niggli.is_xray_amplitude_array():
      miller_niggli = miller_niggli.f_as_f_sq()
    return [ miller_niggli, self.pg_low_prim_set_name, self.pg_choice,
             self.pg_lattice_name ]

  def _show_impl(self, out):
    out.show_header("Exploring higher metric symmetry")
    if self.sigma_warning is not None:
      out.show_text(self.sigma_warning)
    if out.gui_output : # XXX GUI only
      out.show_sub_header("Point group and R-factor analysis")
    out.show_lines("""
The point group of data as dictated by the space group is %s
The point group in the niggli setting is %s
The point group of the lattice is %s
A summary of R values for various possible point groups follow.
""" % (str(self.pg_input_name),
       str(sgtbx.space_group_info(group=self.pg_low_prim_set)),
       str(self.pg_lattice_name)))
    out.show_table(self.table, indent=2)
    # Phenix GUI hack
    if hasattr(out, "add_change_symmetry_button"):
      out.add_change_symmetry_button()
    out.show_lines("""
R_used: mean and maximum R value for symmetry operators *used* in this point group
R_unused: mean and minimum R value for symmetry operators *not used* in this point group
""")
    out.show("""
An automated point group suggestion is made on the basis of the BIC (Bayesian
information criterion).
""")
    out.show("The likely point group of the data is: %s\n" % self.pg_choice)
    out.show("Possible space groups in this point group are:")
    for sg in self.sg_possibilities:
      sg_out = StringIO()
      sg[0].show_summary(f=sg_out, prefix= "   ")
      out.show_preformatted_text(sg_out.getvalue())
    out.show("""
Note that this analysis does not take into account the effects of twinning.
If the data are (almost) perfectly twinned, the symmetry will appear to be
higher than it actually is.
""")


class r_values(scaling.xtriage_analysis):
  def __init__(self,
               miller_obs,
               twin_law,
               miller_calc=None,
               n_reflections=400):
    self.obs = miller_obs.deep_copy()
    self.calc=None
    self.n_reflections=n_reflections
    if miller_calc is not None:
      self.calc = miller_calc.deep_copy()
      self.obs, self.calc = self.obs.common_sets( self.calc )
    self.twin_law = twin_law
    self.rvsr_interpretation=None
    self.table=None
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
    del self.obs
    del self.calc

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

  def _show_impl(self, out=sys.stdout):
    out.show_paragraph_header("R vs R statistics")
    out.show("  R_abs_twin = <|I1-I2|>/<|I1+I2|>")
    out.show("    (Lebedev, Vagin, Murshudov. Acta Cryst. (2006). D62, 83-95)")
    out.show("  R_abs_twin observed data   : %4.3f" % (self.r_abs_obs))
    if self.r_abs_calc is not None:
      out.show("   R_abs_twin calculated data : %4.3f" % self.r_abs_calc)
    out.show("  R_sq_twin = <(I1-I2)^2>/<(I1+I2)^2>")
    out.show("  R_sq_twin observed data    : %4.3f" % self.r_sq_obs)
    if self.r_abs_calc is not None :
      out.show("  R_sq_twin calculated data  : %4.3f" % self.r_sq_calc)
    else:
      out.show("  No calculated data available.")
      out.show("  R_twin for calculated data not determined.")

########################################################################
# OVERALL ANALYSES
########################################################################

class twin_results_interpretation(scaling.xtriage_analysis):
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
    # patterson analyses
    self.patterson_height = self.patterson_p_value = None
    if translational_pseudo_symmetry is not None:
      self.patterson_height = translational_pseudo_symmetry.high_peak
      self.patterson_p_value = translational_pseudo_symmetry.high_p_value
    # wilson statistics moments, etc
    self.i_ratio = wilson_ratios.acentric_i_ratio
    self.f_ratio = wilson_ratios.acentric_f_ratio
    self.e_sq_minus_1 = wilson_ratios.acentric_abs_e_sq_minus_one
    # L test
    self.l_mean = l_test.mean_l
    self.l_sq_mean = l_test.mean_l2
    self.compute_maha_l()
    # twin dependent tests
    self.n_twin_laws = len(twin_law_related_test)
    self.twin_laws = []
    self.twin_law_type = []
    self.r_obs = []
    self.r_calc = []
    # twin fraction estimates
    self.britton_alpha = []
    self.h_alpha = []
    self.murray_rust_alpha = []
    have_twin_test_results = False
    for twin_item in twin_law_related_test:
      self.twin_laws.append( twin_item.twin_law.operator.r().as_hkl() )
      self.twin_law_type.append( str(twin_item.twin_law.twin_type) )
      if twin_item.results_available:
        have_twin_test_results = True
        self.r_obs.append(twin_item.r_values.r_abs_obs)
        self.r_calc.append(twin_item.r_values.r_abs_calc)
        self.britton_alpha.append(twin_item.britton_test.estimated_alpha)
        self.h_alpha.append(twin_item.h_test.estimated_alpha)
        if twin_item.ml_murray_rust is not None:
          self.murray_rust_alpha.append(
            twin_item.ml_murray_rust.estimated_alpha)
        else:
          self.murray_rust_alpha.append(None)
      else:
        self.r_obs.append(None)
        self.r_calc.append(None)
        self.britton_alpha.append(None)
        self.h_alpha.append(None)
        self.murray_rust_alpha.append(None)
    if self.n_twin_laws > 0 :
      if have_twin_test_results :
        max_tf = -1
        max_index = -1
        for ii, twin_fraction in enumerate(self.britton_alpha):
          if (twin_fraction is not None) and (float(twin_fraction) > max_tf):
            max_tf = twin_fraction
            max_index = ii
        self.most_worrysome_twin_law = max_index
      else:
        tmp_tf = []
        for twin_fraction in self.britton_alpha:
          if twin_fraction is None:
            tmp_tf.append( -100 )
          else:
            tmp_tf.append( twin_fraction )
        self.most_worrysome_twin_law = flex.max_index( flex.double(tmp_tf) )
      self.suspected_point_group = symmetry_issues.pg_choice
      self.possible_sgs = symmetry_issues.sg_possibilities
      self.input_point_group = symmetry_issues.pg_low_prim_set_name
    else:
      self.input_point_group = None
      self.most_worrysome_twin_law = None
      self.suspected_point_group = None
      self.possible_sgs = None

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
    tmp_l = self.l_mean - maha_mean_l
    tmp_l2 = self.l_sq_mean - mama_mean_l2
    maha_distance_l = tmp_l*tmp_l*maha_l +\
                      tmp_l2*tmp_l2*maha_l2 +\
                      tmp_l*tmp_l2*maha_ll2
    maha_distance_l = math.sqrt(maha_distance_l)
    self.maha_l = maha_distance_l

  def has_pseudo_translational_symmetry(self):
    if (self.patterson_p_value is None):
      return None
    return (self.patterson_p_value <= self.patterson_p_cut)

  def has_abnormal_intensity_statistics(self):
    return ((self.maha_l >= self.maha_l_cut) and (self.l_mean < 0.5))

  def has_twinning(self):
    return self.has_abnormal_intensity_statistics() and (self.n_twin_laws > 0)

  def has_higher_symmetry(self):
    return (self.input_point_group != self.suspected_point_group)

  def patterson_verdict(self):
    verdict = ""
    if self.has_pseudo_translational_symmetry():
      verdict = """\
The analyses of the Patterson function reveals a significant off-origin
peak that is %3.2f %% of the origin peak, indicating pseudo-translational
symmetry.  The chance of finding a peak of this or larger height by random
in a structure without pseudo-translational symmetry is equal to %5.4e.""" % \
        (self.patterson_height, self.patterson_p_value)
      if self.i_ratio > 2:
        verdict += """
The detected translational NCS is most likely also responsible for the
elevated intensity ratio.  See the relevant section of the logfile for more
details."""
    elif (self.patterson_p_value is not None):
      verdict = """\
The largest off-origin peak in the Patterson function is %3.2f%% of the
height of the origin peak. No significant pseudotranslation is detected.""" % \
        self.patterson_height
    return verdict

  def show_verdict(self, out):
    # First, TNCS verdict
    out.show("\n%s\n" % self.patterson_verdict())
    # And now the rest:
    # Case 1: intensity statistics are suspicious
    if self.maha_l >= self.maha_l_cut:
      # Case 1a: looks like twinning...
      if self.l_mean < 0.5 :
        out.show("""\
The results of the L-test indicate that the intensity statistics
are significantly different than is expected from good to reasonable,
untwinned data.""")
        # [1a] ...and we have twin laws!
        if (self.n_twin_laws > 0):
          out.show("""
As there are twin laws possible given the crystal symmetry, twinning could
be the reason for the departure of the intensity statistics from normality.
It might be worthwhile carrying out refinement with a twin specific target
function.

Please note however that R-factors from twinned refinement cannot be directly
compared to R-factors without twinning, as they will always be lower when a
twin law is used.  You should also use caution when interpreting the maps from
refinement, as they will have significantly more model bias.
""")
          self.twinning_short=True
          if self.has_higher_symmetry ():
            out.show("""
Note that the symmetry of the intensities suggest that the assumed space group
is too low. As twinning is however suspected, it is not immediately clear if
this is the case.  Careful reprocessing and (twin)refinement for all cases
might resolve this question.""")
        # [1a] ...but no twin laws possible
        else :
          out.show("""
As there are no twin laws possible given the crystal symmetry, there could be
a number of reasons for the departure of the intensity statistics from
normality.  Overmerging pseudo-symmetric or twinned data, intensity to
amplitude conversion problems as well as bad data quality might be possible
reasons.  It could be worthwhile considering reprocessing the data.""")
          self.twinning_short=None

      # Case 1b: not twinning?
      else:
        self.twinning_short=None
        out.show("""\
The results of the L-test indicate that the intensity statistics show more
centric character than is expected for acentric data.""")
        if self.has_pseudo_translational_symmetry():
          out.show("""\
This behavior might be explained by the presence of the detected pseudo
translation.""")
          self.twinning_short=False
    # Case 2: no twinning suspected
    else:
      self.twinning_short=False
      out.show("""\
The results of the L-test indicate that the intensity statistics behave as
expected. No twinning is suspected.""")
      if self.has_higher_symmetry():
        out.show("""\
The symmetry of the lattice and intensity however suggests that the input
input space group is too low. See the relevant sections of the log file for
more details on your choice of space groups.""")

      if (self.n_twin_laws > 0):
        if (not self.has_higher_symmetry()):
          if (self.most_worrysome_twin_law is not None):
            if (self.britton_alpha[ self.most_worrysome_twin_law ]>
                TWIN_FRAC_SIGNIFICANT):
              out.show("""\
The correlation between the intensities related by the twin law
%(twin_law)s with an estimated twin fraction of %(alpha)3.2f is most likely
due to an NCS axis parallel to the twin axis. This can be verified by
supplying calculated data as well.
""" % { "twin_law" :  self.twin_laws[ self.most_worrysome_twin_law ],
        "alpha" : self.britton_alpha[ self.most_worrysome_twin_law ] })
        else:
          out.show("""\
As the symmetry is suspected to be incorrect, it is advisable to reconsider
data processing.""")
          if self.has_pseudo_translational_symmetry():
            out.show("""\
Note however that the presence of translational NCS (and possible rotational
pseudo symmetry parallel to the twin axis) can make the detection of twinning
difficult.  Trying various space groups and twinning hypotheses in structure
refinement might provide an answer.
""")

  def make_sym_op_table(self):
    def as_string(x):
      return format_value("%4.3f", x).strip()
    legend = None
    table_data = []
    if self.r_calc[0]==None:
      legend = ('Operator', 'type', 'R obs.', 'Britton alpha', 'H alpha',
                'ML alpha')
      for item in range( len(self.twin_laws) ):
        tmp = [ self.twin_laws[item],
                self.twin_law_type[item],
                as_string( self.r_obs[item]),
                as_string( self.britton_alpha[item]),
                as_string( self.h_alpha[item]),
                as_string( self.murray_rust_alpha[item]) ]
        table_data.append( tmp )
    else:
      legend = ('Operator', 'type', 'R_abs obs.', 'R_abs calc.',
                'Britton alpha', 'H alpha', 'ML alpha')
      for item in range( len(self.twin_laws) ):
        tmp = [ self.twin_laws[item],
                self.twin_law_type[item],
                as_string(self.r_obs[item]),
                as_string(self.r_calc[item]),
                as_string(self.britton_alpha[item]),
                as_string(self.h_alpha[item]),
                as_string(self.murray_rust_alpha[item]) ]
        table_data.append( tmp )
    self.table = table_utils.simple_table(
      column_headers=legend,
      table_rows=table_data)

  def _show_impl(self, out):
    out.show_header("Twinning and intensity statistics summary")
    out.show_sub_header("Final verdict")
    self.show_verdict(out=out)
    out.show_sub_header("Statistics independent of twin laws")
    out.show_lines("""\
  <I^2>/<I>^2 : %(i_ratio)5.3f  (untwinned: 2.0, perfect twin: 1.5)
  <F>^2/<F^2> : %(f_ratio)5.3f  (untwinned: 0.785, perfect twin: 0.885)
  <|E^2-1|>   : %(e_sq_minus_1)5.3f  (untwinned: 0.736, perfect twin: 0.541)
  <|L|>       : %(l_mean)5.3f  (untwinned: 0.500; perfect twin: 0.375)
  <L^2>       : %(l_sq_mean)5.3f  (untwinned: 0.333; perfect twin: 0.200)
  Multivariate Z score L-test: %(maha_l)5.3f
""" % { "i_ratio" : self.i_ratio,
        "f_ratio" : self.f_ratio,
        "e_sq_minus_1" : self.e_sq_minus_1,
        "l_mean" : self.l_mean,
        "l_sq_mean" : self.l_sq_mean,
        "maha_l" : self.maha_l })
    out.show("""
 The multivariate Z score is a quality measure of the given spread in
 intensities. Good to reasonable data are expected to have a Z score lower
 than 3.5.  Large values can indicate twinning, but small values do not
 necessarily exclude it.  Note that the expected values for perfect twinning
 are for merohedrally twinned structures, and deviations from untwinned will
 be larger for perfect higher-order twinning.
""")
    if len(self.twin_laws)>0:
      out.show_sub_header("Statistics depending on twin laws")
      self.make_sym_op_table()
      out.show_table(self.table, indent=2)
    else:
      out.show("\nNo (pseudo)merohedral twin laws were found.\n")

  # FIXME this uses the Britton test, but the PDB validation server appears to
  # use the H test.  Which is correct?
  def max_twin_fraction(self):
    if (self.most_worrysome_twin_law is not None):
      return self.britton_alpha[ self.most_worrysome_twin_law ]
    return 0

  def summarize_issues(self):
    issues = []
    bad_data = False
    if self.has_abnormal_intensity_statistics():
      if ((self.n_twin_laws > 0) and
          (self.max_twin_fraction() > TWIN_FRAC_SIGNIFICANT)):
        issues.append((2, "Intensity statistics suggest twinning "+
          "(intensities are significantly different from expected for "+
          "normal data) and one or more twin operators show a significant "+
          "twin fraction.", "Statistics depending on twin laws"))
      else :
        issues.append((1, "The intensity statistics look unusual, but "+
          "twinning is not indicated or not possible in the given space "+
          "group.", "Wilson ratio and moments"))
    else :
      issues.append((0, "The intensity statistics look normal, indicating "+
        "that the data are not twinned.", "Wilson ratio and moments"))
      if (self.max_twin_fraction() > TWIN_FRAC_SIGNIFICANT):
        issues.append((1, "One or more twin operators show a significant "+
          "twin fraction but since the intensity statistics do not indicate "+
          "twinning, you may have an NCS rotation axis parallel to a "+
          "crystallographic axis.", "Statistics depending on twin laws"))
      if self.has_higher_symmetry():
        issues.append((1, "One or more symmetry operators suggest that the "+
          "data has a higher "+
          "crystallographic symmetry (%s)." % str(self.suspected_point_group),
          "Point group and R-factor analysis"))
    if self.patterson_height and self.patterson_height > 75:
      issues.append((2, "Translational NCS is present at a level that might "+
        "be a result of a missed centering operation (one or more peaks "+
        "greater than 75% of the origin).", "Patterson analyses"))
    elif self.patterson_height and self.patterson_height > 20:
      issues.append((2, "Translational NCS is present at a level that may "+
        "complicate refinement (one or more peaks greater than 20% of the "+
        "origin)", "Patterson analyses"))
    else :
      issues.append((0, "Translational NCS does not appear to be present.",
        "Patterson analyses"))
    return issues

class twin_analyses(scaling.xtriage_analysis):
  """Perform various twin related tests"""
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
               additional_parameters=None,
               original_data=None,
               completeness_as_non_anomalous=None,
                ):
    assert miller_array.is_unique_set_under_symmetry()
    if out is None:
      out = sys.stdout
    self.max_delta = 3.0
    symm_issue_table = [0.08, 75, 0.08, 75]
    perform_ncs_analyses=False
    n_ncs_bins=7
    sigma_inflation = 1.25
    if additional_parameters is not None:
      sigma_inflation = additional_parameters.missing_symmetry.sigma_inflation
      perform_ncs_analyses = \
        additional_parameters.twinning_with_ncs.perform_analyses
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
    self.d_max = math.sqrt(1./d_star_sq_low_limit)
    self.d_min = math.sqrt(1./d_star_sq_high_limit)
    self.d_hkl_for_l_test = d_hkl_for_l_test

    miller_array = miller_array.resolution_filter(
      math.sqrt(1.0/d_star_sq_low_limit),
      math.sqrt(1.0/d_star_sq_high_limit))
    ## Make sure that we actually have some data
    if (miller_array.indices().size() == 0):
      raise Sorry("No suitable data available after resolution cuts")
    if miller_array.observation_type() is None:
      raise RuntimeError("Observation type unknown")

    # minimize the number of conversions between amplitude and intensities
    # by preparing both in advance
    f_obs = i_obs = miller_array
    if miller_array.is_real_array():
      if miller_array.is_xray_intensity_array():
        f_obs = miller_array.f_sq_as_f()
      else :
        i_obs = miller_array.f_as_f_sq()
    else:
      raise RuntimeError("Observations should be a real array.")
    # now do the same for calculated data
    f_calc = i_calc = miller_calc
    if (miller_calc is not None):
      if miller_calc.is_xray_amplitude_array():
        i_calc = miller_calc.f_as_f_sq()
      else :
        f_calc = miller_calc.f_sq_as_f()

    ## Determine possible twin laws
    possible_twin_laws = twin_laws(
      miller_array=f_obs,
      lattice_symmetry_max_delta=self.max_delta)
    ##-----------------------------
    self.possible_twin_laws = possible_twin_laws
    self.normalised_intensities = wilson_normalised_intensities(
      miller_array=i_obs,
      normalise=normalise,
      out=out,
      verbose=verbose)
    ## Try to locat e pseudo translational symm.
    ## If no refls are available at low reso,
    ## an exception is thrown and caught here not to disturb things too much
    self.translational_pseudo_symmetry = None
    try:
      self.translational_pseudo_symmetry = detect_pseudo_translations(
        miller_array=f_obs,
        out=out,
        completeness_as_non_anomalous=completeness_as_non_anomalous,
        verbose=verbose)
    except Sorry: pass

    self.abs_sg_anal = None
    try:
      if miller_array.sigmas() is not None:
        # Look at systematic absences please
        #from mmtbx.scaling import absences
        from mmtbx.scaling import absences
        self.abs_sg_anal = absences.protein_space_group_choices(
          miller_array = self.normalised_intensities.all,
          threshold = 3.0,
          sigma_inflation=sigma_inflation,
          original_data=original_data)#.show(out)
    except Sorry: pass

    centric_cut = self.normalised_intensities.centric
    acentric_cut = self.normalised_intensities.acentric
    self.wilson_moments = wilson_moments(
      acentric_cut,
      centric_cut)
    self.nz_test = n_z_test(
      normalised_acentric=acentric_cut,
      normalised_centric=centric_cut)
    self.l_test=None
    parity_h = parity_k = parity_l = 2
    if (not d_hkl_for_l_test in [None, Auto]):
      parity_h, parity_k, parity_l = d_hkl_for_l_test
    if self.translational_pseudo_symmetry is not None:
      if (d_hkl_for_l_test in [None, Auto]):
        parity_h = self.translational_pseudo_symmetry.mod_h
        parity_k = self.translational_pseudo_symmetry.mod_k
        parity_l = self.translational_pseudo_symmetry.mod_l
      self.l_test = l_test(
        miller_array=acentric_cut,
        parity_h=parity_h,
        parity_k=parity_k,
        parity_l=parity_l)
    else:
      self.l_test = l_test(
        miller_array=acentric_cut,
        parity_h=parity_h,
        parity_k=parity_k,
        parity_l=parity_l)
    ##--------------------------
    self.n_twin_laws = len(possible_twin_laws.operators)
    self.twin_law_dependent_analyses = []
    self.twin_law_names = []
    for ii in range(self.n_twin_laws):
      twin_law_name = possible_twin_laws.operators[ii].operator.r().as_hkl()
      self.twin_law_names.append(twin_law_name)
      tmp_twin_law_stuff = twin_law_dependent_twin_tests(
        twin_law=possible_twin_laws.operators[ii],
        miller_array=i_obs,
        out=out,
        verbose=verbose,
        miller_calc=i_calc,
        normalized_intensities=self.normalised_intensities.acentric,
        ncs_test=perform_ncs_analyses,
        n_ncs_bins=n_ncs_bins)
      self.twin_law_dependent_analyses.append( tmp_twin_law_stuff )
    # now we can check for space group related issues
    self.check_sg = None
    self.suggested_space_group=None
    if self.n_twin_laws > 0:
      self.check_sg = symmetry_issues(
        miller_array=i_obs,
        max_delta=self.max_delta,
        out=out,
        sigma_inflation=sigma_inflation)
      nig_data, pg_this_one, pg_choice, pg_high = \
        self.check_sg.return_point_groups()
      xs_choice = crystal.symmetry(
        unit_cell=nig_data.unit_cell(),
        space_group_symbol=pg_choice,
        assert_is_compatible_unit_cell=False )
      xs_high   = crystal.symmetry(
        unit_cell=nig_data.unit_cell(),
        space_group_symbol=pg_high,
        assert_is_compatible_unit_cell=False )
      if pg_choice != pg_high:
        # FIXME
        merge_data_and_guess_space_groups(
          miller_array=nig_data,
          xs=xs_high,
          out=out,
          txt="Merging in *highest possible* point group %s.\n ***** THIS MIGHT NOT BE THE BEST POINT GROUP SYMMETRY *****  "%pg_high,
          check_absences=False)

      if pg_choice != pg_this_one:
        suggested_space_group = merge_data_and_guess_space_groups(
          miller_array=nig_data,
          xs=xs_choice,
          out=out,
          txt="Merging in *suggested* point group %s "%pg_choice,
          check_absences=False)
    ##--------------------------
    self.twin_summary = twin_results_interpretation(
      nz_test=self.nz_test,
      wilson_ratios=self.wilson_moments,
      l_test=self.l_test,
      translational_pseudo_symmetry=self.translational_pseudo_symmetry,
      twin_law_related_test=self.twin_law_dependent_analyses,
      symmetry_issues=self.check_sg)

  def _show_impl(self, out):
    if (self.abs_sg_anal):
      self.abs_sg_anal.show(out)
    out.show_header("Diagnostic tests for twinning and pseudosymmetry")
    out.show("Using data between %4.2f to %4.2f Angstrom." % (self.d_max,
      self.d_min))
    if (self.translational_pseudo_symmetry is not None):
      self.translational_pseudo_symmetry.show(out)
    self.wilson_moments.show(out)
    self.nz_test.show(out)
    if (self.l_test is not None):
      self.l_test.show(out)
    out.show_header("Twin laws")
    self.possible_twin_laws.show(out=out)
    if (self.n_twin_laws > 0):
      out.show_sub_header("Twin law-specific tests")
      out.show("""\
 The following tests analyze the input data with each of the possible twin
 laws applied.  If twinning is present, the most appropriate twin law will
 usually have a low R_abs_twin value and a consistent estimate of the twin
 fraction (significantly above 0) from each test.  The results are also
 compiled in the summary section.

 WARNING: please remember that the possibility of twin laws, and the results
 of the specific tests, does not guarantee that twinning is actually present
 in the data.  Only the presence of abnormal intensity statistics (as judged
 by the Wilson moments, NZ-test, and L-test) is diagnostic for twinning.
""")
      for twin_tests in self.twin_law_dependent_analyses :
        twin_tests.show(out=out)
      if (self.check_sg is not None):
        self.check_sg.show(out=out)
    self.twin_summary.show(out=out)

  # XXX grossness ahead, beware
  # Objects of this class are relatively bulky due to the storage of multiple
  # derived Miller arrays that may be used elsewhere.  Since we do not need
  # these arrays for simply displaying the results (in the Phenix GUI or
  # elsewhere), they are deleted prior to pickling to reduce the amount of
  # data that needs to be transfered or saved.  It is not necessary to
  # implement __setstate__, since we are still just pickling self.__dict__.
  def __getstate__(self):
    """
    Pickling function with storage efficiency optimizations.
    """
    if (self.abs_sg_anal is not None):
      self.abs_sg_anal.miller_array = None
      self.abs_sg_anal.absences_table.miller_array = None
    if (self.check_sg is not None):
      self.check_sg.miller_niggli = None
      self.check_sg.miller_array = None
    self.normalised_intensities = None
    return self.__dict__

def merge_data_and_guess_space_groups(miller_array, txt, xs=None,out=None,
    sigma_inflation=1.0, check_absences=True):
  tmp_ma = miller_array.deep_copy()
  if xs is None:
    xs = tmp_ma.crystal_symmetry()
  tmp_ma = miller_array.customized_copy( crystal_symmetry=xs )
  merge_obj = tmp_ma.change_basis(
      xs.space_group_info().change_of_basis_op_to_reference_setting()
    ).merge_equivalents()
  tmp_ma = merge_obj.array()
  r_lin = merge_obj.r_linear()
  normalizer = absolute_scaling.kernel_normalisation(tmp_ma, auto_kernel=True)
  work_array = normalizer.normalised_miller.deep_copy()
  abs_sg_anal = None
  if tmp_ma.sigmas() is not None and (check_absences):
    print(file=out)
    print(file=out)
    print("-"*len(txt), file=out)
    print(txt, file=out)
    print("-"*len(txt), file=out)
    print(file=out)
    merge_obj.show_summary(out=out)
    print(file=out)
    print("Suggesting various space group choices on the basis of systematic absence analyses", file=out)
    print(file=out)
    print(file=out)
    this_worked=False
    try:
      if miller_array.sigmas() is not None:
        # Look at systematic absences please
        #from mmtbx.scaling import absences
        from mmtbx.scaling import absences
        abs_sg_anal = absences.protein_space_group_choices(
          miller_array = work_array,
          threshold = 3.0,
          print_all=False,
          sigma_inflation=sigma_inflation).show(out)
    except Sorry:
      print("Systematic absence analyses failed", file=out)
  return (merge_obj, abs_sg_anal)

########################################################################
# MILLER ARRAY EXTENSIONS
# Injector class to extend the Miller array class with the intensity analyses
# contained in this module.
def analyze_intensity_statistics(self, d_min=2.5,
    completeness_as_non_anomalous=None, log=None):
  """
  Detect translational pseudosymmetry and twinning.  Returns a
  twin_law_interpretation object.
  """
  if (log is None) : log = null_out()
  if self.space_group().is_centric():
    return None
  tmp_array = self.resolution_filter(d_min=d_min)
  if (not self.sigmas_are_sensible()):
    tmp_array = tmp_array.customized_copy(
      indices=tmp_array.indices(),
      data=tmp_array.data(),
      sigmas=None).set_observation_type( tmp_array )
  twin_results = twin_analyses(
    miller_array=tmp_array,
    d_star_sq_low_limit=1.0/100.0, # XXX need to confirm
    d_star_sq_high_limit=1.0/(0.001**2.0), # XXX need to confirm
    completeness_as_non_anomalous=completeness_as_non_anomalous,
    out = null_out(),
    out_plots = null_out(),
    verbose=False)
  summary = twin_results.twin_summary
  summary.show(out=log)
  return summary

def twin_analyses_brief(miller_array,
                        cut_off=2.5,
                        completeness_as_non_anomalous=None,
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
  if (out is None) and (verbose != 0):
    out = sys.stdout
  summary = miller_array.analyze_intensity_statistics(d_min=cut_off,
    completeness_as_non_anomalous=completeness_as_non_anomalous,
    log=out, )
  return summary.has_twinning()
