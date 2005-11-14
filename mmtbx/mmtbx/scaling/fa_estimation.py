from cctbx import maptbx
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx import adptbx
import cctbx.sgtbx.lattice_symmetry
import cctbx.sgtbx.cosets
from cctbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
import mmtbx.scaling
from mmtbx.scaling import absolute_scaling, relative_scaling
from mmtbx.scaling import matthews, twin_analyses
from mmtbx.scaling import basic_analyses, data_statistics,pair_analyses
import libtbx.phil.command_line
from cStringIO import StringIO
from scitbx.python_utils import easy_pickle
import sys, os


class combined_scaling(object):
  def __init__(self,
               miller_array_x1,
               miller_array_x2,
               options=None,
               out=None):
    ## miller_array_x1 : 'reference'
    ## miller_array_x2 : 'derivative' or something similar
    ## options : a phil file with scaling options
    ## convergence : rescale after rejection?
    ##
    ##
    ## These are the tasks ahead
    ## 0) map to asu, common sets
    ## 1) least squares scaling of x1 and x2
    ## 2) local scaling of x1 and x2
    ## 3) outlier rejections
    ## 4) least squares scaling of x1 and x2
    ## 5) local scaling of x1 and x2
    ## 6) generation of delta F's
    ##
    ## we will perform operations on this arrays
    self.x1 = miller_array_x1.deep_copy().map_to_asu()
    self.x2 = miller_array_x2.deep_copy().map_to_asu()
    ## Get the scaling and rejection options
    self.lsq_options = options.least_squares_options
    self.loc_options = options.local_scaling_options
    self.out_options = options.outlier_rejection_options
    self.overall_protocol = options.target
    self.cut_level_rms = None
    ## cycling options
    self.auto_cycle=False
    if options.iterations=='auto':
      self.auto_cycle=True

    self.max_iterations = options.max_iterations
    if not self.auto_cycle:
      assert self.max_iterations > 0

    ## output options
    self.out=out
    if self.out==None:
      self.out=sys.stdout
    ## get common sets, and make a reference copy please
    self.x1, self.x2 = self.x1.common_sets( self.x2 )
    self.s1 = self.x1.deep_copy().map_to_asu()
    self.s2 = self.x2.deep_copy().map_to_asu()


    scaling_tasks={'lsq':True, 'local':True }

    if self.overall_protocol=='ls':
      scaling_tasks['lsq']=True
      scaling_tasks['local']=False

    if self.overall_protocol=='loc':
      scaling_tasks['lsq']=False
      scaling_tasks['local']=True

    if self.overall_protocol=='ls_and_loc':
      scaling_tasks['lsq']=True
      scaling_tasks['local']=True

    assert ( scaling_tasks['lsq'] or scaling_tasks['local'] )

    self.convergence=False
    counter = 0

    print >> self.out, "=========================================="
    print >> self.out, "=             Relative scaling           ="
    print >> self.out, "=========================================="



    while not self.convergence:
      print >> self.out
      print >> self.out, "--------------------------"
      print >> self.out, "    Scaling cycle %i   "%(counter)
      print >> self.out, "--------------------------"

      if counter == 0:
        self.cut_level_rms = self.out_options.cut_level_rms_primary
      else:
        self.cut_level_rms = self.out_options.cut_level_rms_secondary

      ## take the common sets
      self.s1 = self.s1.common_set(  self.x1 )
      self.s2 = self.s2.common_set(  self.x2 )
      self.x1, self.x2 = self.s1.common_sets( self.s2 )


      if scaling_tasks['lsq']:
        self.perform_least_squares_scaling()
      if scaling_tasks['local']:
        self.perform_local_scaling()
      num_reject = self.perform_outlier_rejection()
      if num_reject==0:
        self.convergence=True
      if not self.auto_cycle:
        if counter==self.max_iterations:
          self.convergence=True
      counter+=1

    ## Now the datasets have been scaled, we can return bnotyh datasets
    ## They can be used by the usewr to computer either
    ## - isomorphous differences
    ## - be used in an F1-F2 FFT
    del self.s1
    del self.s2




  def perform_least_squares_scaling(self):
    print >> self.out
    print >> self.out, "Least squares scaling"
    print >> self.out, "---------------------"
    print >> self.out

    ##-----Get options-----
    use_exp_sigmas=self.lsq_options.use_experimental_sigmas

    use_int = True
    if self.lsq_options.scale_data=='amplitudes':
      use_int=False

    use_wt=True
    if self.lsq_options.scale_target == 'basic':
      use_wt=False
    ##-------------------------------------------
    ls_scaling = relative_scaling.ls_rel_scale_driver(
        self.x1,
        self.x2,
        use_intensities=use_int,
        scale_weight=use_wt,
        use_weights=use_exp_sigmas)
    ls_scaling.show(out=self.out)
    ##----- Update the miller arrays please-------
    self.x1 = ls_scaling.native.deep_copy()
    self.x2 = ls_scaling.derivative.deep_copy()

  def perform_local_scaling(self):
    print >> self.out
    print >> self.out, "Local scaling"
    print >> self.out, "-------------"
    print >> self.out

    ##-----Get options-----
    use_exp_sigmas=self.loc_options.use_experimental_sigmas
    use_int = True

    if self.loc_options.scale_data=='amplitudes':
      use_int=False

    moment_based_local_scale=False
    if self.loc_options.scale_target=='local_moment':
      moment_based_local_scale=True
    ##--------------------
    local_scaling = relative_scaling.local_scaling_driver(
      self.x1,
      self.x2,
      use_intensities=use_int,
      use_weights=use_exp_sigmas,
      moment_based=moment_based_local_scale,
      max_depth=self.loc_options.max_depth,
      target_neighbours=self.loc_options.target_neighbours,
      sphere=self.loc_options.neighbourhood_sphere,
      out=self.out)

    self.x1 = local_scaling.native.deep_copy()
    self.x2 = local_scaling.derivative.deep_copy()

  def perform_outlier_rejection(self):
    print >> self.out
    print >> self.out, "Outlier rejections"
    print >> self.out, "------------------"
    print >> self.out
    print >> self.out, " sigma criterion : %4.1f "%(self.out_options.cut_level_sigma)
    print >> self.out, " rms criterion   : %4.1f "%(self.cut_level_rms)
    print >> self.out, " protocol        :", self.out_options.protocol
    print >> self.out

    outlier_protocol ={'solve':False,'rms':False, 'rms_and_sigma':False }
    ## Please set the protocol to what is specified
    outlier_protocol[ self.out_options.protocol ]=True

    outlier_rejection = pair_analyses.outlier_rejection(
      self.x1,
      self.x2,
      cut_level_rms=self.cut_level_rms,
      cut_level_sigma=self.out_options.cut_level_sigma
    )

    number_of_outliers = self.x1.size() - outlier_rejection.nat.size()
    number_of_outliers += self.x2.size() - outlier_rejection.der.size()

    self.x1 = outlier_rejection.nat.deep_copy()
    self.x2 = outlier_rejection.der.deep_copy()

    self.x1, self.x2 = self.x1.common_sets( self.x2 )
    return( number_of_outliers )




class ano_scaling(object):
  def __init__(self,
               miller_array_x1,
               options,
               out=None):
    ## These are the tasks ahead
    ##
    ## 1) splitting up x1 in hemsispheres x1p x1n
    ## 2) sumbiut the two halve data sets to the combined scaler
    ##

    assert miller_array_x1.anomalous_flag()
    assert miller_array_x1.indices().size() > 0

    self.s1p, self.s1n = miller_array_x1.hemispheres_acentrics()

    self.s1p = self.s1p.set_observation_type( miller_array_x1 )

    self.s1n = self.s1n.customized_copy( indices=-self.s1n.indices() )
    self.s1n = self.s1n.set_observation_type( miller_array_x1 )


    assert self.s1p.indices().size() == self.s1n.indices().size()

    self.x1p = self.s1p.deep_copy()
    self.x1n = self.s1n.deep_copy()

    ## Now we have a 'native' and a 'derivative'
    ##
    ## Submit these things too the combined scaler
    ano_scaler=combined_scaling(
      self.x1p,
      self.x1n,
      options,
      out)
