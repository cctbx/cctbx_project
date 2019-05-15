from __future__ import absolute_import, division, print_function
from six.moves import range
from scipy.optimize import leastsq # special import
import math,copy
import scitbx.math
from cctbx.sgtbx import tensor_rank_2_constraints
from scitbx.array_family import flex
from xfel.metrology.mark0 import correction_vectors
from xfel.metrology.mark3 import fit_translation2
from xfel.metrology import mark3_collect_data
from xfel import get_radial_tangential_vectors
from scitbx.lbfgs.tst_curvatures import lbfgs_with_curvatures_mix_in
from rstbx.bandpass import parameters_bp3
from scitbx.matrix import col
from cctbx.crystal import symmetry
from math import pi
from scitbx import matrix
from xfel.metrology.legacy_scale import mark5_iteration,vector_collection
from xfel.metrology.legacy_scale import bandpass_gaussian
from rstbx.symmetry.constraints import AGconvert
#XXX remove the xfel/metrology/legacy_scale wrapping for AGconvert; not needed

class fit_translation4(mark5_iteration,fit_translation2):
  def parameter_based_model(self,params):
    PIXEL_SZ = 0.11 # mm/pixel
    all_model = mark3_collect_data(self.frame_id, self.HKL)
    self.FRAMES["refined_detector_origin"] = flex.vec3_double(len(self.FRAMES["frame_id"]))

    for iframe in range(len(self.FRAMES["frame_id"])):
      frame_id = self.FRAMES["frame_id"][iframe]
      self.frame_id_to_param_no[frame_id] = iframe

      detector_origin = self.parameter_based_model_one_frame_detail(frame_id,iframe,all_model)
      try:
        self.bandpass_models[frame_id].gaussian_fast_slow()
      except Exception as e:
        print("Exception from picture",e)
        raise e

      try:
       all_model.collect_mean_position(self.bandpass_models[frame_id].mean_position,
                        self.bandpass_models[frame_id].observed_flag,
                        frame_id);
       all_model.collect_distance(self.bandpass_models[frame_id].part_distance,frame_id)
       self.FRAMES["refined_detector_origin"][iframe] = detector_origin/(-PIXEL_SZ)

      except Exception as e:
          print("Exception from collect",e)
          raise e

      #print "HATTNE check 0", list(self.bandpass_models[frame_id].mean_position), len(self.bandpass_models[frame_id].mean_position)
      #print "HATTNE check 1", list(self.bandpass_models[frame_id].observed_flag), len(self.bandpass_models[frame_id].observed_flag)


    sq_displacements = ((all_model.cx - self.spotcx)*(all_model.cx - self.spotcx) +
                        (all_model.cy - self.spotcy)*(all_model.cy - self.spotcy))


#    print "In parameter_based_model, HATTNE got sqdisplacements", flex.mean(sq_displacements)
#    print "  * am, am", flex.mean(all_model.cx), flex.mean(all_model.cy), len(all_model.cx), len(all_model.cy)
#    print "  * cx, cy", flex.mean(self.spotcx), \
#                        flex.mean(self.spotcy), \
#                        len(self.spotcx), len(self.spotcy)
#    print "  * diffs ", flex.mean(all_model.cx - self.spotcx), flex.mean(all_model.cy - self.spotcy), flex.max(flex.abs(all_model.cx - self.spotcx)), flex.max(flex.abs(all_model.cy - self.spotcy))

#    from matplotlib import pyplot as plt
#    plt.plot(self.spotfx, self.spotfy, "r.")
##    plt.plot(self.spotcx, self.spotcy, "g.")
#    plt.plot(all_model.cx, all_model.cy, "g.")
#    plt.show()


    #print list(all_model.cx - self.spotcx)



    selected_sq_displacements = sq_displacements.select( all_model.flags == True )
    #print "Root Mean squared displacement all spots      %8.3f"%math.sqrt(
    #  flex.sum(selected_sq_displacements)/len(selected_sq_displacements))
    return all_model.cx, all_model.cy, all_model.flags, all_model.part_distance

  def detector_origin_analysis(self):
    self.FRAMES["detector_origin_x_refined"]=flex.double()
    self.FRAMES["detector_origin_y_refined"]=flex.double()
    self.FRAMES["distance_refined"]=flex.double()
    for iframe in range(len(self.FRAMES["frame_id"])):
      if iframe < self.n_refined_frames:
        SIGN = -1.
        PIXEL_SZ = 0.11 # mm/pixel
        detector_origin = col((-self.FRAMES["beam_x"][iframe]
                               + SIGN * PIXEL_SZ * self.frame_translations.x[2*iframe],
                               -self.FRAMES["beam_y"][iframe]
                               + SIGN * PIXEL_SZ * self.frame_translations.x[1+2*iframe],
                               0.))
        self.FRAMES["detector_origin_x_refined"].append(detector_origin[0])
        self.FRAMES["detector_origin_y_refined"].append(detector_origin[1])
        self.FRAMES["distance_refined"].append(
           self.frame_distances.x[iframe] +
           self.FRAMES["distance"][iframe]
        )

    xm = flex.mean_and_variance(self.FRAMES["detector_origin_x_refined"])
    ym = flex.mean_and_variance(self.FRAMES["detector_origin_y_refined"])
    print("Beam x mean %7.3f sigma %7.3f mm"%(
      xm.mean(), xm.unweighted_sample_standard_deviation()))
    print("Beam y mean %7.3f sigma %7.3f mm"%(
      ym.mean(), ym.unweighted_sample_standard_deviation()))

    time_series = False
    import os
    files = [os.path.basename(f) for f in self.FRAMES["unique_file_name"]]
    longs = [long("".join([a for a in name if a.isdigit()]))//1000 for name in files]
    floats = flex.double([float(L) for L in longs])[
      :len(self.FRAMES["detector_origin_x_refined"])]

    order = flex.sort_permutation(floats)
    time_sorted_x_beam = self.FRAMES["detector_origin_x_refined"].select(order)
    time_sorted_y_beam = self.FRAMES["detector_origin_y_refined"].select(order)

    if time_series:
      from matplotlib import pyplot as plt
      plt.plot(range(len(order)),time_sorted_x_beam,"r-")
      plt.plot(range(len(order)),time_sorted_y_beam,"b-")
      plt.show()

    for item in order:
      print(files[item], "%8.3f %8.3f dist %8.3f"%(
        self.FRAMES["detector_origin_x_refined"][item],
        self.FRAMES["detector_origin_y_refined"][item],
        self.FRAMES["distance_refined"][item]))

  def __init__(self,params):
    self.optional_params = None
    self.frame_id_to_param_no = {}
    correction_vectors.__init__(self)
    # XXX Set in mark0.correction_vectors.read_data()
    #self.nominal_tile_centers(params.effective_tile_boundaries)
    self.params = params
    self.read_data(params)
    # for eventual write-out, deep copy the input:
    self.OUTPUT = copy.deepcopy(self.FRAMES)

    if params.max_frames==None:
      FR_LIMIT = len(self.FRAMES["frame_id"]) # don't adjust for degrees of freedom
      self.n_refined_frames = FR_LIMIT
    else:
      self.n_refined_frames = params.max_frames

    # XXX CHECK THIS BLOCK MORE CAREFULLY!
    mark5_iteration.__init__(self)
    self.tile_translations=self.register(
      "tile_trans",ndata=len(self.tiles) // 4,kdim=2,data=[0.]*(len(self.tiles) // 2)) # x & y displacements for each of 64 tiles
    self.tile_rotations=self.register("tile_rot",ndata=len(self.tiles) // 4,data=[0.]*(len(self.tiles) // 4))
    self.frame_translations=self.register(
      "frame_trans",ndata=self.n_refined_frames,kdim=2,data=[0.]*(2*self.n_refined_frames))
    self.frame_distances=self.register("frame_distance",
      ndata=self.n_refined_frames,data=[0.]*self.n_refined_frames)
    self.frame_rotz=self.register("frame_rotz",
      ndata=self.n_refined_frames,data=[0.]*self.n_refined_frames)

    self.frame_roty=self.register_local("frame_roty",
      ndata=self.n_refined_frames,data=[0.]*self.n_refined_frames)
    self.frame_rotx=self.register_local("frame_rotx",
      ndata=self.n_refined_frames,data=[0.]*self.n_refined_frames)
    self.mean_energy_factor=self.register_local("mean_energy_factor",
      ndata=self.n_refined_frames,data=[1.]*self.n_refined_frames)
    self.bandpass_logfac=self.register_local("bandpass_logfac",
      ndata=self.n_refined_frames,data=[0.]*self.n_refined_frames)
    self.mosaicity_factor=self.register_local("mosaicity_factor",
      ndata=self.n_refined_frames,data=[1.]*self.n_refined_frames)
    self.g_factor=self.register_local("g_factor",
      ndata=self.n_refined_frames,kdim=6,data=[1.]*6*self.n_refined_frames)

    self.vector_data = vector_collection(self.frame_id, self.HKL)
    self.jacobian_frame_roty = self.vector_data.register(tag="frame_roty")
    self.jacobian_frame_rotx = self.vector_data.register(tag="frame_rotx")
    """0. Read in only 500 frames, refine only 20
       1. Define the new parameter here.
       2. Define its new behavior in parameter_based_function()
       3. set use_curvatures=False
       4. implement refinement with finite differences:
            include the array name in setting compute_finite_difference_gradients_if_requested()
            fix the requirement for "zip"
            verify fd_gradients go to zero upon convergence
            verify parameter vector has reasonable values
       5. implement analytic gradients, list out comparison columns
       6. document the gradients in a tex file
       7. implement curvatures & flag them in here
    """
    self.x = self.as_x_array()





  def compute_finite_difference_gradients_if_requested(self):
    for param_array,tag in zip([self.frame_rotz],["frame_rotz"]):
      result = flex.double()
      for np in range(len(param_array.x)):
        reserve = param_array.x[np]
        epsilon = 0.00001
        param_array.x[np] = reserve + epsilon
        hicx, hicy, hiflags, hipart_distance = self.parameter_based_model(self.params)
        parameter_nos = flex.int([self.frame_id_to_param_no[f] for f in self.frame_id])

        functional_hi = self.compute_functional_only(
          tox = self.To_x, toy = self.To_y,
          spotcx = hicx, spotcy = hicy,
          spotfx = self.spotfx, spotfy = self.spotfy,
          master_tiles = self.master_tiles,
          frames = parameter_nos,
          part_distance = hipart_distance)

        param_array.x[np] = reserve - epsilon
        locx, locy, loflags, lopart_distance = self.parameter_based_model(self.params)

        functional_lo = self.compute_functional_only(
          tox = self.To_x, toy = self.To_y,
          spotcx = locx, spotcy = locy,
          spotfx = self.spotfx, spotfy = self.spotfy,
          master_tiles = self.master_tiles,
          frames = parameter_nos,
          part_distance = lopart_distance)

        #print "GRADIENT ",dist_parm," finite diff",(plus_funct-nega_funct)/(2.*epsilon),"analytical",newgrad[dist_parm]
        param_array.x[np]=reserve
        result.append( ( functional_hi - functional_lo )/(2.*epsilon) )

      print("Setting fd grads", end=' ')
      for a in result:  print("%8.3f"%a, end=' ')
      print()
      print("       on values", end=' ')
      for a in param_array.x:  print("%8.3f"%a, end=' ')
      print()
      print("   cpp gradients", end=' ')
      for a in param_array.gradients: print("%8.3f"%a, end=' ')
      print()
      print()
      self.set_gradient_array(tag,result)

  def compute_functional_and_gradients(self):

    #print "HATTNE in compute_functional_and_gradients"

    self.from_x_array(self.x)
    self.cx, self.cy, self.flags, part_distance = self.parameter_based_model(self.params)

    #print "HATTNE in compute_functional_and_gradients", flex.mean(self.cx), flex.mean(self.cy), len(self.cx)


    #import sys
    #sys.exit(0) # HATTNE





    parameter_nos = flex.int([self.frame_id_to_param_no[f] for f in self.frame_id])
    self.rezero_gradients_curvatures()
    self.set_refined_origins_to_c(self.FRAMES["refined_detector_origin"])
    functional = self.compute_target(
       tox = self.To_x, toy = self.To_y,
       spotcx = self.cx.deep_copy(), spotcy = self.cy.deep_copy(),
       spotfx = self.spotfx, spotfy = self.spotfy,
       master_tiles = self.master_tiles,
       frames = parameter_nos,
       part_distance = part_distance)

    #self.compute_finite_difference_gradients_if_requested()

    gradients = self.get_gradient_array()
    self.c_curvatures = self.get_curvature_array()

#    print "CURVATURES", len(self.c_curvatures)
#    lst = list(self.c_curvatures)
#    print " block 1: ", lst[0:20]
#    print " block 2: ", lst[20:20 + 58]
#    print " block 3: ", lst[20 + 58:20 + 2 * 58]
#    print " block 4: ", lst[20 + 2* 58:20 + 3 * 58]

    print("Functional ",functional, math.sqrt(functional/self.cx.size()))
    return functional,gradients

  def radial_transverse_analysis(self):
    # very time consuming; should be recoded (C++?)
    SIGN = -1.
    PIXEL_SZ = 0.11 # mm/pixel
    zunit = col([0.,0.,1.])
    self.frame_del_radial = {}
    self.frame_del_azimuthal = {}
    for x in range(len(self.frame_id)):
      frame_id = self.frame_id[x]
      frame_param_no = self.frame_id_to_param_no[frame_id]
      if frame_param_no >= self.n_refined_frames: continue
      if frame_id not in self.frame_del_radial:
        self.frame_del_radial[frame_id] = flex.double()
        self.frame_del_azimuthal[frame_id] = flex.double()

      correction_vector = col([self.correction_vector_x[x],self.correction_vector_y[x],0.])
      position_vector = col([self.spotfx[x] -
                      SIGN * self.FRAMES["detector_origin_x_refined"][frame_param_no]/PIXEL_SZ,
                             self.spotfy[x] -
                      SIGN * self.FRAMES["detector_origin_y_refined"][frame_param_no]/PIXEL_SZ,
                             0.])
      position_unit = position_vector.normalize()
      transverse_unit = position_unit.cross(zunit)

      self.frame_del_radial[frame_id].append(correction_vector.dot(position_unit))
      self.frame_del_azimuthal[frame_id].append(correction_vector.dot(transverse_unit))

  def print_table(self):
    from libtbx import table_utils
    from libtbx.str_utils import format_value
    table_header = ["Tile","Dist","Nobs","aRmsd","Rmsd","delx","dely","disp","rotdeg",
                    "Rsigma","Tsigma","Transx","Transy","DelRot","Rotdeg"]
    table_data = []
    table_data.append(table_header)
    sort_radii = flex.sort_permutation(flex.double(self.radii))
    tile_rmsds = flex.double()
    radial_sigmas = flex.double(len(self.tiles) // 4)
    tangen_sigmas = flex.double(len(self.tiles) // 4)

    wtaveg = [0.]*(len(self.tiles) // 4)
    for x in range(len(self.tiles) // 4):
      if self.tilecounts[x] >= 3:
        wtaveg[x] = self.weighted_average_angle_deg_from_tile(x, self.post_mean_cv[x], self.correction_vector_x,
          self.correction_vector_y)

    for idx in range(len(self.tiles) // 4):
      x = sort_radii[idx]
      if self.tilecounts[x] < 3:
        radial = (0,0)
        tangential = (0,0)
        rmean,tmean,rsigma,tsigma=(0,0,1,1)
      else:
        radial,tangential,rmean,tmean,rsigma,tsigma = get_radial_tangential_vectors(self,x,
          self.post_mean_cv[x],
          self.correction_vector_x, self.correction_vector_y,
          self.model_calcx-self.refined_cntr_x,
          self.model_calcy-self.refined_cntr_y)

      # paired rotations of two ASICS on the same sensor
      if x%2==0:
        # previous method: delrot = "%5.2f"%(wtaveg[x]-wtaveg[x+1])
        delrot = "%5.2f"%(self.tile_rotations.x[x] - self.tile_rotations.x[1+x])
      else:
        delrot = ""

      radial_sigmas[x]=rsigma
      tangen_sigmas[x]=tsigma
      table_data.append(  [
        format_value("%3d",   x),
        format_value("%7.2f", self.radii[x]),
        format_value("%6d",  self.tilecounts[x]),
        format_value("%5.2f", self.asymmetric_tile_rmsd[x]),
        format_value("%5.2f", self.tile_rmsd[x]),
        format_value("%5.2f", self.post_mean_cv[x][0]),
        format_value("%5.2f", self.post_mean_cv[x][1]),
        format_value("%5.2f", matrix.col(self.post_mean_cv[x]).length()),
        format_value("%6.2f", wtaveg[x]),
        format_value("%6.2f", rsigma),
        format_value("%6.2f", tsigma),
        format_value("%5.2f", self.tile_translations.x[2*x]),
        format_value("%5.2f", self.tile_translations.x[2*x+1]),
        copy.copy(delrot),
        format_value("%5.2f", self.tile_rotations.x[x])
      ])
    table_data.append([""]*len(table_header))
    rstats = flex.mean_and_variance(radial_sigmas,self.tilecounts.as_double())
    tstats = flex.mean_and_variance(tangen_sigmas,self.tilecounts.as_double())
    table_data.append(  [
        format_value("%3s",   "ALL"),
        format_value("%s", ""),
        format_value("%6d",  self.overall_N),
        format_value("%5.2f", math.sqrt(flex.mean(self.delrsq))),
        format_value("%5.2f", self.overall_rmsd),
        format_value("%5.2f", self.overall_cv[0]),
        format_value("%5.2f", self.overall_cv[1]),
        format_value("%5.2f", flex.mean(flex.double([cv.length() for cv in self.post_mean_cv]))),
        format_value("%s", ""),
        format_value("%6.2f", rstats.mean()),
        format_value("%6.2f", tstats.mean()),
        format_value("%s", ""),
        format_value("%s", ""),
        #root mean squared difference in same-sensor (adjacent)-ASIC rotations, weighted by minimum # of observations on either ASIC of the sensor
        format_value("%5.2f", math.sqrt(
           flex.sum(
             flex.double([
               (min([self.tilecounts[2*isen],self.tilecounts[2*isen+1]])) *
                    (self.tile_rotations.x[2*isen] - self.tile_rotations.x[1+2*isen])**2
               for isen in range(len(self.tiles) // 8)]
             )
           )/
           flex.sum(
             flex.double(
               [(min([self.tilecounts[2*isen],self.tilecounts[2*isen+1]])) for isen in range(len(self.tiles) // 8)]
             )
           )
        )),
        format_value("%s", ""),
    ])

    print()
    print(table_utils.format(table_data,has_header=1,justify='center',delim=" "))

  def print_table_2(self):

    from libtbx import table_utils
    from libtbx.str_utils import format_value
    table_header = ["Tile","Dist","Nobs","aRmsd","Rmsd","delx","dely","disp","rotdeg",
                    "Rsigma","Tsigma","Transx","Transy","DelRot","Rotdeg"]
    table_data = []
    table_data.append(table_header)
    sort_radii = flex.sort_permutation(flex.double(self.radii))
    tile_rmsds = flex.double()
    radial_sigmas = flex.double(len(self.tiles) // 4)
    tangen_sigmas = flex.double(len(self.tiles) // 4)

    wtaveg = [0.]*(len(self.tiles) // 4)
    for x in range(len(self.tiles) // 4):
      if self.tilecounts[x] >= 3:
        wtaveg[x] = self.weighted_average_angle_deg_from_tile(x, self.post_mean_cv[x], self.correction_vector_x,
          self.correction_vector_y)

    def add_line_to_table(idx):
      x = sort_radii[idx]
      if self.tilecounts[x] < 3:
        radial = (0,0)
        tangential = (0,0)
        rmean,tmean,rsigma,tsigma=(0,0,1,1)
      else:
        radial,tangential,rmean,tmean,rsigma,tsigma = get_radial_tangential_vectors(self,x,
          self.post_mean_cv[x],
          self.correction_vector_x, self.correction_vector_y,
          self.model_calcx-self.refined_cntr_x,
          self.model_calcy-self.refined_cntr_y)

      table_data.append(  [
        format_value("%3d",   x),
        format_value("%7.2f", self.radii[x]),
        format_value("%6d",  self.tilecounts[x]),
        format_value("%5.2f", self.asymmetric_tile_rmsd[x]),
        format_value("%5.2f", self.tile_rmsd[x]),
        format_value("%5.2f", self.post_mean_cv[x][0]),
        format_value("%5.2f", self.post_mean_cv[x][1]),
        format_value("%5.2f", matrix.col(self.post_mean_cv[x]).length()),
        format_value("%6.2f", wtaveg[x]),
        format_value("%6.2f", rsigma),
        format_value("%6.2f", tsigma),
        format_value("%5.2f", self.tile_translations.x[2*x]),
        format_value("%5.2f", self.tile_translations.x[2*x+1]),
        "",
        format_value("%5.2f", self.tile_rotations.x[x])
      ])

    # order the printout by sensor, starting from innermost
    new_order = []
    mutable = list(sort_radii)
    idx = 0
    unit_translation_increments = flex.double(len(mutable)*2)
    while 1:
      if idx >= len(mutable): break
      if self.radii[mutable[idx]]==0.0:
        idx+=1; continue
      tile_select = mutable[idx]
      if tile_select%2 == 0:
        # even
        sensor_tiles = (tile_select, tile_select+1)
        sensor_ptrs = (idx, mutable.index(tile_select+1))
      else:
        # odd
        sensor_tiles = (tile_select-1, tile_select)
        sensor_ptrs = ( mutable.index(tile_select-1), idx)

      if self.tilecounts[mutable[sensor_ptrs[0]]] + self.tilecounts[mutable[sensor_ptrs[1]]] < \
         self.params.min_count:
         idx+=1
         continue

      sum_weight = 0.0
      sum_wt_x = 0.0
      sum_wt_y = 0.0
      for iptr, ptr in enumerate(sensor_ptrs):
        if ptr in new_order: break
        if self.tilecounts[mutable[ptr]] > 0:
          #print mutable[ptr]
          add_line_to_table (ptr)
          sum_weight += self.tilecounts[mutable[ptr]]
          sum_wt_x += self.tilecounts[mutable[ptr]] * self.tile_translations.x[2*mutable[ptr]]
          sum_wt_y += self.tilecounts[mutable[ptr]] * self.tile_translations.x[2*mutable[ptr]+1]
        new_order.append(ptr)
        if iptr==1:
          #print
          sensor_line = [""]*len(table_header)
          sensor_line[2]="%6d"%sum_weight
          sensor_line[11]="%5.2f"%round(sum_wt_x/sum_weight,0)
          sensor_line[12]="%5.2f"%round(sum_wt_y/sum_weight,0)
          unit_translation_increments[2*mutable[ptr]-2] = round(sum_wt_x/sum_weight,0)
          unit_translation_increments[2*mutable[ptr]-1] = round(sum_wt_y/sum_weight,0)
          unit_translation_increments[2*mutable[ptr]] = round(sum_wt_x/sum_weight,0)
          unit_translation_increments[2*mutable[ptr]+1] = round(sum_wt_y/sum_weight,0)
          table_data.append(sensor_line)
          table_data.append([""]*len(table_header))
      idx+=1
      if idx>=len(mutable): break

    print("Grouped by sensor, listing lowest Q-angle first:")
    print(table_utils.format(table_data,has_header=1,justify='center',delim=" "))
    return unit_translation_increments

  def same_sensor_table(self,verbose=True):
    radii = flex.double() # from-instrument-center distance in pixels
    delrot= flex.double() # delta rotation in degrees
    meanrot=flex.double() # mean rotation in degrees
    weight= flex.double() #
    displacement = [] # vector between two same-sensor ASICS in pixels
    for x in range(len(self.tiles) // 8):
      delrot.append(self.tile_rotations.x[2*x] - self.tile_rotations.x[1+2*x])
      meanrot.append(0.5*(self.tile_rotations.x[2*x] + self.tile_rotations.x[1+2*x]))
      radii.append((self.radii[2*x]+self.radii[2*x+1])/2)
      weight.append(min([self.tilecounts[2*x],self.tilecounts[2*x+1]]))
      displacement.append(   col((self.To_x[2*x+1], self.To_y[2*x+1]))
                            -col((self.tile_translations.x[2*(2*x+1)], self.tile_translations.x[2*(2*x+1)+1]))
                            -col((self.To_x[2*x], self.To_y[2*x]))
                            +col((self.tile_translations.x[2*(2*x)], self.tile_translations.x[2*(2*x)+1]))  )

    unrotated_displacement = [] # same, except correct for the off-square rotation of this sensor
    for i,a in enumerate(displacement):
      corrected_a = a.rotate_2d(angle=meanrot[i],deg=True)
      while (corrected_a[0] < 0. or abs(corrected_a[1]) > abs(corrected_a[0])):
        corrected_a = corrected_a.rotate_2d(angle=90., deg=True)
      unrotated_displacement.append( corrected_a )

    order = flex.sort_permutation(radii)
    if verbose:
      for x in order:
        print("%02d %02d %5.0f"%(2*x,2*x+1,weight[x]), end=' ')
        print("%6.1f"%radii[x], end=' ')
        print("%5.2f"%(delrot[x]), end=' ')
        print("%6.3f"%(displacement[x].length()-194.), end=' ') # ASIC is 194; just print gap
        #print "  %6.3f"%(self.tile_distances.x[x])
        print("lateral %7.3f transverse %7.3f pix"%(unrotated_displacement[x][0], unrotated_displacement[x][1]))
    stats = flex.mean_and_variance(flex.double([t.length()-194. for t in displacement]),weight)
    print(" sensor gap is %7.3f px +/- %7.3f"%(stats.mean(), stats.gsl_stats_wsd()))
    stats = flex.mean_and_variance(flex.double([t[0] for t in unrotated_displacement]),weight)
    print("lateral gap is %7.3f px +/- %7.3f"%(stats.mean(), stats.gsl_stats_wsd()))
    stats = flex.mean_and_variance(flex.double([t[1] for t in unrotated_displacement]),weight)
    print("transverse gap is %7.3f px +/- %7.3f"%(stats.mean(), stats.gsl_stats_wsd()))

  @staticmethod
  def print_unit_translations(data, params, optional):
    from scitbx.array_family import flex
    def pretty_format(data):
      out = """"""
      for quad in [0,1,2,3]:
        for blockof2 in [0,1,2,3]:
          format = "%3.0f,"*8
          if quad==3 and blockof2==3:
            format = format[0:-1]
          format = "     %s"%format
          out+=format+"""
"""
        if quad<3: out += """
"""
      # reality check here in case the PAD is not 64-tiles
      Nparam = len(data)
      decoration = out.split("%3.0f,")
      if len(decoration)>Nparam+1:
        decoration = decoration[len(decoration)-Nparam:]
        return "%3.0f,".join(decoration)

      return out

    from spotfinder.applications.xfel.cxi_phil import cxi_versioned_extract
    if params.detector_format_version is not None:
      base_arguments = ["distl.detector_format_version=%s"%params.detector_format_version]

      if optional is not None:

        if optional.distl.quad_translations is not None:
          base_arguments.append("distl.quad_translations="+",".join([str(s) for s in optional.distl.quad_translations]))
        if optional.distl.tile_translations is not None:
          base_arguments.append("distl.tile_translations="+",".join([str(s) for s in optional.distl.tile_translations]))

      stuff = cxi_versioned_extract(base_arguments)
      old = flex.double(stuff.distl.tile_translations)
      print("cctbx already defines unit pixel translations for detector format version %s:"%params.detector_format_version)
      print(pretty_format(old)%tuple(old))
      print("new unit pixel increments will be SUBTRACTED off these to get final translations")
      print()

    else:
      print("no pre-existing translations were input")
      print()
      old = flex.double(128)

    print("Unit translations to be pasted into spotfinder/applications/xfel/cxi_phil.py:")

    new = old - flex.double(data)
    #overall_format = """    working_extract.distl.tile_translations = [
    overall_format = '''distl {\n  tile_translations = """
'''+pretty_format(new)+'''    """}'''

    print(overall_format%tuple(new))


  def run_cycle_a(self):
    self.bandpass_models = {}


    #print "HATTNE MARKER #0"
    lbfgs_with_curvatures_mix_in.__init__(self,
      min_iterations=0,
      max_iterations=1000,
      traditional_convergence_test_eps=1.0, # new for LD91; to shorten the run time
      use_curvatures=True)
    #print "HATTNE MARKER #1"

    C = self
    C.post_min_recalc()

    #import sys
    #sys.exit(0) # HATTNE

    C.print_table()
    unit_translations = C.print_table_2()
    self.print_unit_translations(unit_translations, self.params, self.optional_params)

    #import sys
    #sys.exit(0) # HATTNE

    sum_sq = 0.
    for key in C.frame_delx.keys():
      param_no = C.frame_id_to_param_no[key]
      if param_no >= C.n_refined_frames:continue
      mn_x = flex.mean(C.frame_delx[key])
      mn_y = flex.mean(C.frame_dely[key])
      print("frame %4d count %4d delx %7.2f  dely %7.2f"%(key,
        len(C.frame_delx[key]),
        mn_x,
        mn_y ), end=' ')
      sum_sq += mn_x*mn_x + mn_y*mn_y
      if param_no<C.n_refined_frames:
        print("%7.2f %7.2f"%(C.frame_translations.x[2*param_no],C.frame_translations.x[1+2*param_no]))
      else:  print("N/A")
    displacement = math.sqrt(sum_sq / len(C.frame_delx))
    print("rms displacement of frames %7.2f"%displacement)

    C.detector_origin_analysis()
    C.radial_transverse_analysis()
    sum_sq_r = 0.
    sum_sq_a = 0.
    for key in C.frame_del_radial.keys():
      mn_r = flex.mean(C.frame_del_radial[key])
      mn_a = flex.mean(C.frame_del_azimuthal[key])
      print("frame %4d count %4d delrad %7.2f  delazi %7.2f"%(key,
        len(C.frame_del_radial[key]),
        mn_r,
        mn_a ), end=' ')
      sum_sq_r += mn_r*mn_r
      sum_sq_a += mn_a*mn_a
      param_no = C.frame_id_to_param_no[key]
      print("deldist %7.3f mm"%C.frame_distances.x[param_no], end=' ')

      print("distance %7.3f mm"%(
        C.frame_distances.x[param_no] +
        C.FRAMES["distance"][param_no]), end=' ')

      print("rotz_deg=%6.2f"%( (180./math.pi)*C.frame_rotz.x[param_no] ))
    print("rms radial displacement %7.2f"%(math.sqrt(sum_sq_r/len(C.frame_del_radial))))
    print("rms azimut displacement %7.2f"%(math.sqrt(sum_sq_a/len(C.frame_del_radial))))
    C.same_sensor_table()

    print("Cycle A translations & rotations")
    print({"translations": list(self.tile_translations.x),
           "rotations": list(self.tile_rotations.x)})
    print("integration {")
    print("  subpixel_joint_model{")
    print("rotations= \\")
    for irot in range(len(self.tile_rotations.x)):
      print(" %11.8f "%self.tile_rotations.x[irot], end=' ')
      if irot%4==3 and irot!=len(self.tile_rotations.x)-1: print("\\")
    print()
    print("translations= \\")
    for irot in range(len(self.tile_translations.x)):
      print(" %11.8f"%self.tile_translations.x[irot], end=' ')
      if irot%4==3 and irot!=len(self.tile_translations.x)-1: print("\\")
    print()
    print("  }")
    print("}")
    print()



  def run_cycle_b(self,iteration):
    for iframe in range(len(self.FRAMES["frame_id"])):
      frame_id = self.FRAMES["frame_id"][iframe]
      self.frame_id_to_param_no[frame_id] = iframe

    self.bandpass_models = {}
    all_model = mark3_collect_data(self.frame_id, self.HKL)
    for iframe in range(min(self.n_refined_frames,len(self.FRAMES["frame_id"]))):
      frame_id = self.FRAMES["frame_id"][iframe]
      file_name = self.FRAMES["unique_file_name"][iframe]

      self.parameter_based_model_one_frame_detail(frame_id,iframe,all_model)
      #instantiate a helper class for holding the per-frame data
      class per_frame_helper:
        def __init__(pfh):
          pfh.master_tiles = flex.int()
          pfh.spotfx = flex.double()
          pfh.spotfy = flex.double()
          pfh.cosine = [math.cos(self.tile_rotations.x[tidx]*(math.pi/180))
                        for tidx in range(len(self.To_x))]
          pfh.sine   = [math.sin(self.tile_rotations.x[tidx]*(math.pi/180))
                        for tidx in range(len(self.To_x))]
          self.parameter_based_model_one_frame_detail(frame_id,iframe,all_model)
          self.bandpass_models[frame_id].gaussian_fast_slow()
          mean_position = self.bandpass_models[frame_id].mean_position
          first_index = all_model.get_first_index(frame_id)
          for ridx in range(mean_position.size()):
            pfh.master_tiles.append( self.master_tiles[first_index + ridx] );
            pfh.spotfx.append( self.spotfx[first_index + ridx] );
            pfh.spotfy.append( self.spotfy[first_index + ridx] );
        def fvec_callable(pfh,current_values,frame_id,iframe,all_model):

          #print "HATTNE entering fvec_callable()"

          self.frame_roty.x[iframe] = current_values[0]
          self.frame_rotx.x[iframe] = current_values[1]
#          self.mean_energy_factor.x[iframe] = current_values[2]
#          self.bandpass_logfac.x[iframe] = current_values[3]
#          self.mosaicity_factor.x[iframe] = current_values[4]
          for iparam in range(self.bandpass_models["n_independent"]):
            self.g_factor.x[iparam+6*iframe] = current_values[2+iparam]
          try:
            self.parameter_based_model_one_frame_detail(frame_id,iframe,all_model)
          except Exception as e:
            print("Failed on", iframe, self.FRAMES["unique_file_name"][iframe])
            raise e

          #print "  HATTNE marker #0", frame_id, type(self.bandpass_models), self.bandpass_models.keys(), type(self.bandpass_models[frame_id])

          self.bandpass_models[frame_id].gaussian_fast_slow()
          #print "  HATTNE marker #1"
          mean_position = self.bandpass_models[frame_id].mean_position

          translations    = self.tile_translations.x
          fval = []
          for ridx in range(mean_position.size()):
            itile = pfh.master_tiles[ridx];

            calc_minus_To_x = mean_position[ridx][1] - self.To_x[itile];
            calc_minus_To_y = mean_position[ridx][0] - self.To_y[itile];

            rotated_o_x = calc_minus_To_x * pfh.cosine[itile]\
                        - calc_minus_To_y * pfh.sine[itile];
            rotated_o_y = calc_minus_To_x * pfh.sine[itile]\
                        + calc_minus_To_y * pfh.cosine[itile];

            model_calcx = rotated_o_x + (self.To_x[itile] + translations[2*itile]);
            model_calcy = rotated_o_y + (self.To_y[itile] + translations[2*itile+1]);

            delx = model_calcx - pfh.spotfx[ridx];
            dely = model_calcy - pfh.spotfy[ridx];

            fval.append(delx)
            fval.append(dely)

          cum = 0.
          for item in fval:
            cum += item*item

          rmsd = math.sqrt( cum / (len(fval)/2) )
          print("rmsd", rmsd)
          return fval

      self.helper = per_frame_helper()

      print("Trying iframe",iframe,"independent=%d"%self.bandpass_models["n_independent"], end=' ')

      results = leastsq(
        func = self.helper.fvec_callable,
        x0 = [self.frame_roty.x[iframe],self.frame_rotx.x[iframe],
#              self.mean_energy_factor.x[iframe],self.bandpass_logfac.x[iframe],
#              self.mosaicity_factor.x[iframe]
             ] + [self.g_factor.x[n+6*iframe]
                  for n in range(self.bandpass_models["n_independent"])],
        args = (frame_id,iframe,all_model),
        Dfun = None, #estimate the Jacobian
        full_output = True)

      print("with %d reflections"%self.bandpass_models[frame_id].mean_position.size(), end=' ')
      print("result %6.2f degrees"%(results[0][0]*180./math.pi), end=' ')
      print("result %6.2f degrees"%(results[0][1]*180./math.pi), end=' ')
      #modify this line to deal with cubic space group
      print("energy factor %6.4f"%(results[0][2]),"metrical%1d %7.5f %7.5f"%(iteration,results[0][2],results[0][3]), end=' ')
      fff = results[2]["fvec"]
      cum = 0.
      for item in fff:
        cum += item*item
      rmsd = math.sqrt( cum / (len(fff)/2) )
      print("rmsd", rmsd, end=' ')
      print(file_name)
      #print results[4]

  def parameter_based_model_one_frame_detail(self,frame_id,iframe,all_model):
      PIXEL_SZ = 0.11 # mm/pixel
      SIGN = -1.
      if iframe < self.n_refined_frames:
        detector_origin = col((-self.FRAMES["beam_x"][iframe]
                             + SIGN * PIXEL_SZ * self.frame_translations.x[2*iframe],
                             -self.FRAMES["beam_y"][iframe]
                             + SIGN * PIXEL_SZ * self.frame_translations.x[1+2*iframe],
                             0.))
        self.OUTPUT["beam_x"][iframe] = -detector_origin[0]
        self.OUTPUT["beam_y"][iframe] = -detector_origin[1]
      else:
        detector_origin = col((-self.FRAMES["beam_x"][iframe],-self.FRAMES["beam_y"][iframe],0.))

      if frame_id not in self.bandpass_models:

        reserve_orientation = self.FRAMES["orientation"][iframe]
        effective_orientation = reserve_orientation

        #Not necessary to apply the 3 offset rotations; they have apparently
        #  been applied already.\
        #  .rotate_thru((1,0,0),self.FRAMES["rotation100_rad"][iframe]
        # ).rotate_thru((0,1,0),self.FRAMES["rotation010_rad"][iframe]
        # ).rotate_thru((0,0,1),self.FRAMES["rotation001_rad"][iframe])

        crystal = symmetry(unit_cell=effective_orientation.unit_cell(),space_group = "P1")
        indices = all_model.frame_indices(frame_id)

        parameters = parameters_bp3(
           indices=indices, orientation=effective_orientation,
           incident_beam=col(correction_vectors.INCIDENT_BEAM),
           packed_tophat=col((1.,1.,0.)),
           detector_normal=col(correction_vectors.DETECTOR_NORMAL),
           detector_fast=col((0.,1.,0.)),detector_slow=col((1.,0.,0.)),
           pixel_size=col((PIXEL_SZ,PIXEL_SZ,0)),
           pixel_offset=col((0.,0.,0.0)),
           distance=self.FRAMES["distance"][iframe],
           detector_origin=detector_origin
        )

        #print "PARAMETER check   ", effective_orientation
        #print "PARAMETER distance", self.FRAMES['distance'][iframe]
        #print "PARAMETER origin  ", detector_origin

        ucbp3 = bandpass_gaussian(parameters=parameters)
        ucbp3.set_active_areas( self.tiles ) #self.params.effective_tile_boundaries
        integration_signal_penetration=0.0 # easier to calculate distance derivatives

        ucbp3.set_sensor_model( thickness_mm = 0.5, mu_rho = 8.36644, # CS_PAD detector at 1.3 Angstrom
          signal_penetration = integration_signal_penetration)
        #ucbp3.set_subpixel( flex.double(tp038_trans_values) ) #back off this; let minimizer figure it out.

        half_mosaicity_rad = self.FRAMES["half_mosaicity_deg"][iframe] * pi/180.
        ucbp3.set_mosaicity(half_mosaicity_rad)
        ucbp3.set_bandpass(self.FRAMES["wave_HE_ang"][iframe],self.FRAMES["wave_LE_ang"][iframe])
        ucbp3.set_orientation(effective_orientation)
        ucbp3.set_domain_size(self.FRAMES["domain_size_ang"][iframe])
        ucbp3.set_vector_output_pointers(self.vector_data,
                                         frame_id,iframe<self.n_refined_frames)

        if "best_index" not in self.bandpass_models:
          from labelit.dptbx import lepage
          M = lepage.character(effective_orientation)
          s = len(M.best())
          for index in M.best():
            index['counter'] = s
            s-=1
            if index["max_angular_difference"]==0.0:
              best_index = index
              break

          self.bandpass_models["best_index"] = best_index
          self.bandpass_models["constraints"] = tensor_rank_2_constraints(space_group=best_index['reduced_group'],reciprocal_space=True)
          self.bandpass_models["n_independent"] = self.bandpass_models["constraints"].n_independent_params()

        self.bandpass_models[frame_id]=ucbp3

      if iframe < self.n_refined_frames:
        self.bandpass_models[frame_id].set_detector_origin(detector_origin)
        self.bandpass_models[frame_id].set_distance(
          self.FRAMES["distance"][iframe] + self.frame_distances.x[iframe])
        self.OUTPUT["distance"][iframe] = self.FRAMES["distance"][iframe] + self.frame_distances.x[iframe]
        #half_mosaicity_rad = self.FRAMES["half_mosaicity_deg"][iframe] * pi/180. + \
        #                     self.half_mosaicity_rad.x[iframe]
        #self.bandpass_models[frame_id].set_mosaicity(half_mosaicity_rad)
        reserve_orientation = self.FRAMES["orientation"][iframe]
        effective_orientation =   reserve_orientation.rotate_thru((0,0,1),self.frame_rotz.x[iframe])
        effective_orientation = effective_orientation.rotate_thru((0,1,0),self.frame_roty.x[iframe])
        effective_orientation = effective_orientation.rotate_thru((1,0,0),self.frame_rotx.x[iframe])

        convert = AGconvert()
        convert.forward(effective_orientation)
        u_independent = list(self.bandpass_models["constraints"].independent_params(all_params=convert.G))
        for x in range(self.bandpass_models["n_independent"]):
          u_independent[x] *= self.g_factor.x[x+6*iframe]
        u_star = self.bandpass_models["constraints"].all_params(independent_params=tuple(u_independent))
        convert.validate_and_setG(u_star)
        effective_orientation = convert.back_as_orientation()
        self.OUTPUT["orientation"][iframe]=effective_orientation
        self.bandpass_models[frame_id].set_orientation(effective_orientation)
        mean_wave = (self.FRAMES["wave_HE_ang"][iframe] + self.FRAMES["wave_LE_ang"][iframe])/2.
        #mean_wave *= self.mean_energy_factor.x[iframe]
        bandpassHW =(self.FRAMES["wave_LE_ang"][iframe] - self.FRAMES["wave_HE_ang"][iframe])/2.
        self.bandpass_models[frame_id].set_bandpass(mean_wave - bandpassHW, mean_wave + bandpassHW)

      return detector_origin

if (__name__ == "__main__"):

  result = run(args=[#"mysql.runtag=L785v_120corner","mysql.passwd=sql789",
                     #"mysql.user=nick","mysql.database=xfelnks",
                     #"show_plots=True",
                     #"show_consistency=True",
                     "max_frames=1001",
#                     "effective_tile_boundaries=713, 437, 907, 622, 516, 438, 710, 623, 713, 650, 907, 835, 516, 650, 710, 835, 509, 18, 694, 212, 507, 215, 692, 409, 721, 19, 906, 213, 720, 215, 905, 409, 86, 231, 280, 416, 283, 231, 477, 416, 85, 19, 279, 204, 283, 19, 477, 204, 106, 444, 291, 638, 106, 640, 291, 834, 318, 443, 503, 637, 318, 640, 503, 834, 434, 849, 619, 1043, 436, 1046, 621, 1240, 647, 848, 832, 1042, 648, 1045, 833, 1239, 18, 1066, 212, 1251, 214, 1065, 408, 1250, 17, 853, 211, 1038, 213, 853, 407, 1038, 229, 1474, 414, 1668, 229, 1277, 414, 1471, 15, 1474, 200, 1668, 16, 1278, 201, 1472, 442, 1464, 636, 1649, 638, 1464, 832, 1649, 441, 1252, 635, 1437, 638, 1252, 832, 1437, 846, 1134, 1040, 1319, 1042, 1133, 1236, 1318, 845, 922, 1039, 1107, 1042, 922, 1236, 1107, 1060, 1542, 1245, 1736, 1060, 1346, 1245, 1540, 848, 1543, 1033, 1737, 847, 1348, 1032, 1542, 1469, 1336, 1663, 1521, 1272, 1337, 1466, 1522, 1472, 1550, 1666, 1735, 1274, 1549, 1468, 1734, 1460, 1117, 1645, 1311, 1460, 921, 1645, 1115, 1247, 1117, 1432, 1311, 1248, 921, 1433, 1115, 1130, 718, 1315, 912, 1130, 522, 1315, 716, 918, 719, 1103, 913, 917, 523, 1102, 717, 1543, 514, 1737, 699, 1346, 513, 1540, 698, 1543, 725, 1737, 910, 1346, 725, 1540, 910, 1338, 94, 1523, 288, 1339, 290, 1524, 484, 1552, 93, 1737, 287, 1551, 289, 1736, 483, 1115, 114, 1309, 299, 918, 113, 1112, 298, 1115, 326, 1309, 511, 918, 326, 1112, 511".replace(" ",""),
#                     "effective_tile_boundaries=719, 434, 913, 619, 521, 434, 715, 619, 719, 648, 913, 833, 521, 648, 715, 833, 515, 13, 700, 207, 512, 210, 697, 404, 727, 14, 912, 208, 726, 211, 911, 405, 90, 225, 284, 410, 288, 226, 482, 411, 89, 13, 283, 198, 288, 13, 482, 198, 109, 439, 294, 633, 109, 637, 294, 831, 322, 439, 507, 633, 322, 637, 507, 831, 438, 847, 623, 1041, 440, 1045, 625, 1239, 652, 847, 837, 1041, 653, 1045, 838, 1239, 21, 1064, 215, 1249, 217, 1064, 411, 1249, 20, 850, 214, 1035, 217, 851, 411, 1036, 232, 1474, 417, 1668, 232, 1277, 417, 1471, 17, 1474, 202, 1668, 18, 1277, 203, 1471, 446, 1464, 640, 1649, 642, 1465, 836, 1650, 445, 1252, 639, 1437, 643, 1252, 837, 1437, 851, 1134, 1045, 1319, 1048, 1134, 1242, 1319, 851, 921, 1045, 1106, 1048, 922, 1242, 1107, 1065, 1543, 1250, 1737, 1066, 1347, 1251, 1541, 853, 1544, 1038, 1738, 852, 1349, 1037, 1543, 1475, 1338, 1669, 1523, 1278, 1339, 1472, 1524, 1477, 1550, 1671, 1735, 1279, 1550, 1473, 1735, 1466, 1118, 1651, 1312, 1466, 921, 1651, 1115, 1253, 1118, 1438, 1312, 1255, 921, 1440, 1115, 1137, 717, 1322, 911, 1137, 520, 1322, 714, 924, 718, 1109, 912, 923, 521, 1108, 715, 1550, 513, 1744, 698, 1353, 512, 1547, 697, 1550, 725, 1744, 910, 1353, 725, 1547, 910, 1346, 91, 1531, 285, 1347, 288, 1532, 482, 1560, 92, 1745, 286, 1558, 288, 1743, 482, 1122, 111, 1316, 296, 925, 109, 1119, 294, 1122, 324, 1316, 509, 925, 323, 1119, 508".replace(" ",""),

  ])
  """mark0: no minimization; just evaluate tiles like cxi.plot_cv
            reports fewer spots than cxi.plot_cv because it ignores frames with < 10 spots not
            producing a spot output.
     mark1: lsq fitting of translation of spot predictions, to bring them on top of
             the fictitious tile montage containing spotfinder spots.
     mark2: mark1 plus tile rotations following the translation step; spinning around tile center.
            Controls: Same-sensor tile rotations should be equal; found rmsd difference 0.02 degrees
     mark3: same as mark2, but predictions come from the model instead of the log file.
            The refined tile positions include all sub-pixel corrections (replacing tp038), instead of using the
            log file predictions that implicitly include the tp038 model.
            Controls:
              Agreement of labelit-refined direct beam (frame table) with CV logfile (spotfinder table)
              Agreement of integrated spot position(observation table) w/ CV obs (spotfinder table)
                 ... in both approximate position as well as Miller index tag, when the tp038 model
                     is included.
              Same-sensor ASIC displacement is 2.55 pixels with sigma=0.26 pixels.  It was thought
              that this sigma should be much closer to zero, and this suggests that the current
              sensor placement is is not as accurate as hoped.
     mark4: Independently refine the beam position of each shot.  Before
              refinement, the rmsd displacement is 0.65 pixels.  After, 0.00 pixels.
              Total spot position rmsd is reduced from 1.36 pixels to 1.23 pixels.
              Displacement, per frame:  beam x (sigma=0.037 mm) beam y (sigma=0.027 mm)
                            radial, rms=0.42 pixels, azimuthal, rms=0.50 pixels
              Brings global spot rmsd down to 1.23 pixels.
     mark5: Independent refine distance of each shot.  Set penetration to 0.0 so spots have
              an exactly linear dependence on distance, easier to calculate gradient & curvature.
     mark6: 1) independent refinement of rotz for each shot.
            2) Redesign so that any set of parameters can be chosen for refinement.
              --self.x set defined dynamically; no hardcode of parameter order or semantics. DONE
              --Python/C++ translation is automatic; set defined in Python, executed in C++. DONE
              --Object defines its own semantics "refine on frames, refine distance" or
                 "refine on tiles, refine distance". HARD CODED
              --All C++ code should be centralized in one module. IN PROGRESS
              --C++ code can test if refinement is enabled on a parameter; so control with Python constructor DONE
              --Reorganize so it is very easy to define new parameters. IN PROGRESS
              --Also, the reporting code should be more transparent. DONE
              --The refinement class should be refactored ao it is easy to
                test out finite differences on any given paramaeter. IN PROGRESS
              --Facility to add in chain rule easily; reducing # of parameters.
              --Ability to corefine PSII and thermolysin data.
              --Ability to switch between finite diff, analytical, side-by-side,
                and mix&match production analytical with selected finite diff.  IN PROGRESS
              --for same sensor, compute lateral and transverse asic separation. DONE
     mark7: same as mark6, except the bandpass model is re-expressed so the predicted spot position is
            smooth (twice differentiable) w.r.t. bandpass and mosaicity; use Gaussians instead of
            top hat functions.
     mark9: introduce macrocycles.  Each cycle consists of a) refine tile xy, tile rotation,
            frame xy, frame distance, frame rotz by LBFGS using curvatures. b) each frame:  refine
            rotx and roty by Gauss-Newton non-linear least squares. However, due to the difficulty
            in deriving and then implementing derivatives, let the code estimate the Jacobian
            by finite differences.  It turns out that the scipy
            implementation of NL-LSQ has extension modules incompatible with boost python
            (segmentation fault).  Option 1) use Luc Bourhis' code: either Gauss-Newton
            or Levenburg Marqhardt.  Option 2) (used) import scipy optimize before anything else;
            seems to work.
     mark10:Continue using scipy NL-LSQ with global/per-frame macrocycle.
            Refine mean energy factor.  Can't refine bandwidth factor, mosaicity factor without
            a radial spotwidth target (future).  Instead, refine metrical_matrix g-factors for the
            independent parameters (2 for hexagonal thermolysin).  Abandon mean energy refinement;
            too strongly correlated with g-factor for non-gradient method. Achieved 0.89 pixel rmsd.
            Write out refined parameters for reintroduction into data integration script.
  """
"""
   gradients available from the python layer
"""
