from __future__ import absolute_import, division, print_function
from six.moves import range
import math,copy
import scitbx.math
from scitbx.array_family import flex
from xfel.metrology.mark0 import get_phil,correction_vectors
from scitbx import matrix
from xfel import get_radial_tangential_vectors
from xfel.metrology.mark1 import fit_translation
from xfel.metrology import mark2_iteration,mark3_collect_data
from scitbx.lbfgs.tst_curvatures import lbfgs_with_curvatures_mix_in
from rstbx.bandpass import use_case_bp3,parameters_bp3
from scitbx.matrix import col
from cctbx.crystal import symmetry
from math import pi

class fit_translation2(fit_translation):
  def parameter_based_model(self,params):
    PIXEL_SZ = 0.11 # mm/pixel
    all_model = mark3_collect_data(self.frame_id, self.HKL)

    for iframe in range(len(self.FRAMES["frame_id"])):
      frame_id = self.FRAMES["frame_id"][iframe]
      if frame_id not in self.bandpass_models:

        reserve_orientation = self.FRAMES["orientation"][iframe]
        effective_orientation = reserve_orientation

        #Not necessary to apply the 3 offset rotations; they have apparently
        #  been applied already.\
        #  .rotate_thru((1,0,0),self.FRAMES["rotation100_rad"][iframe]
        # ).rotate_thru((0,1,0),self.FRAMES["rotation010_rad"][iframe]
        # ).rotate_thru((0,0,1),self.FRAMES["rotation001_rad"][iframe])

        detector_origin = col((-self.FRAMES["beam_x"][iframe],
                               -self.FRAMES["beam_y"][iframe], 0.))
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
        ucbp3 = use_case_bp3(parameters=parameters)

        ucbp3.set_active_areas( self.tiles ) #params.effective_tile_boundaries
        integration_signal_penetration=0.5

        ucbp3.set_sensor_model( thickness_mm = 0.5, mu_rho = 8.36644, # CS_PAD detector at 1.3 Angstrom
          signal_penetration = integration_signal_penetration)

        half_mosaicity_rad = self.FRAMES["half_mosaicity_deg"][iframe] * pi/180.
        ucbp3.set_mosaicity(half_mosaicity_rad)
        ucbp3.set_bandpass(self.FRAMES["wave_HE_ang"][iframe],self.FRAMES["wave_LE_ang"][iframe])
        ucbp3.set_orientation(effective_orientation)
        ucbp3.set_domain_size(self.FRAMES["domain_size_ang"][iframe])

        ucbp3.picture_fast_slow_force()
        self.bandpass_models[frame_id]=ucbp3

      all_model.collect(self.bandpass_models[frame_id].hi_E_limit,
                        self.bandpass_models[frame_id].lo_E_limit,
                        self.bandpass_models[frame_id].observed_flag,
                        frame_id);

    sq_displacements = ((all_model.cx - self.spotcx)*(all_model.cx - self.spotcx) +
                        (all_model.cy - self.spotcy)*(all_model.cy - self.spotcy))
    selected_sq_displacements = sq_displacements.select( all_model.flags == True )
    print("Root Mean squared displacement all spots      %8.3f"%math.sqrt(
      flex.sum(selected_sq_displacements)/len(selected_sq_displacements)))
    return all_model.cx, all_model.cy, all_model.flags

  def __init__(self,params):
    self.bandpass_models = {}
    correction_vectors.__init__(self)
    #self.nominal_tile_centers(params.effective_tile_boundaries)
    self.read_data(params)
    self.params = params


    self.x = flex.double([0.]*(len(self.tiles) // 2)+[0.]*(len(self.tiles) // 4)) # x & y displacements for each of 64 tiles

    lbfgs_with_curvatures_mix_in.__init__(self,
      min_iterations=0,
      max_iterations=1000,
      use_curvatures=True)

  def curvatures(self):
    return self.c_curvatures

  def nominal_tile_centers(self,corners):
    self.To_x =flex.double(64)
    self.To_y = flex.double(64)
    for x in range(64):
      self.To_x[x] = (corners[4*x] + corners[4*x+2])/2.
      self.To_y[x] = (corners[4*x+1] + corners[4*x+3])/2.

  def compute_functional_and_gradients(self):

    self.cx, self.cy, self.flags = self.parameter_based_model(self.params)
    engine = mark2_iteration(values = self.x, tox = self.To_x,
                                              toy = self.To_y,
                                              spotcx = self.cx.deep_copy(),
                                              spotcy = self.cy.deep_copy(),
#                                              spotcx = self.spotcx,
#                                              spotcy = self.spotcy,
                                              spotfx = self.spotfx,
                                              spotfy = self.spotfy,
                                              master_tiles = self.master_tiles)

    print("Functional ",math.sqrt(engine.f()/self.cx.size()))
    self.c_curvatures = engine.curvatures()
    self.model_calcx = engine.model_calcx
    self.model_calcy = engine.model_calcy


    ### HATTNE ADDITION FOR TILE 14 ###
    #selection = self.selections[14]
    #print "HATTNE additional 1", self.To_x[14], self.To_y[14], len(self.To_x), len(self.To_y)
    #print "HATTNE additional 2", list(self.model_calcx.select(selection))
    #print "HATTNE additional 3", self.c_curvatures[2 * 14], self.c_curvatures[2 * 14 + 1]
    #import sys
    #sys.exit(0)
    ### HATTNE ADDITION FOR TILE 14 ###


    return engine.f(),engine.gradients()

  def post_min_recalc(self):
    fit_translation.post_min_recalc(self)
    self.frame_delx = {}
    self.frame_dely = {}
    for x in range(len(self.frame_id)):
      frame_id = self.frame_id[x]
      if frame_id not in self.frame_delx:
        self.frame_delx[frame_id] = flex.double()
        self.frame_dely[frame_id] = flex.double()
      self.frame_delx[frame_id].append(self.correction_vector_x[x])
      self.frame_dely[frame_id].append(self.correction_vector_y[x])

  def same_sensor_table(self,verbose=True):
    radii = flex.double() # from-instrument-center distance in pixels
    delrot= flex.double() # delta rotation in degrees
    weight= flex.double() #
    displacement = [] # vector between two same-sensor ASICS in pixels
    for x in range(len(self.tiles) // 8):
      delrot.append(self.x[len(self.tiles) // 2 +2*x] - self.x[len(self.tiles) // 2 + 1 +2*x])
      radii.append((self.radii[2*x]+self.radii[2*x+1])/2)
      weight.append(min([self.tilecounts[2*x],self.tilecounts[2*x+1]]))
      displacement.append(   col((self.To_x[2*x+1], self.To_y[2*x+1]))
                            -col((self.x[2*(2*x+1)], self.x[2*(2*x+1)+1]))
                            -col((self.To_x[2*x], self.To_y[2*x]))
                            +col((self.x[2*(2*x)], self.x[2*(2*x)+1]))  )
    order = flex.sort_permutation(radii)
    if verbose:
      for x in order:
        print("%02d %02d %5.0f"%(2*x,2*x+1,weight[x]), end=' ')
        print("%6.1f"%radii[x], end=' ')
        print("%5.2f"%(delrot[x]), end=' ')
        print("%6.3f"%(displacement[x].length()-194.)) # ASIC is 194; just print gap
    stats = flex.mean_and_variance(flex.double([t.length()-194. for t in displacement]),weight)
    print("sensor gap is %7.3f px +/- %7.3f"%(stats.mean(), stats.gsl_stats_wsd()))

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
        delrot = "%5.2f"%(self.x[len(self.tiles) // 2 +x] - self.x[len(self.tiles) // 2 + 1 +x])
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
        format_value("%5.2f", self.x[2*x]),
        format_value("%5.2f", self.x[2*x+1]),
        copy.copy(delrot),
        format_value("%5.2f", self.x[len(self.tiles) // 2 +x])
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
               (min([self.tilecounts[2*isen],self.tilecounts[2*isen+1]])) * (self.x[len(self.tiles) // 2 +2*isen] - self.x[len(self.tiles) // 2 + 1 +2*isen])**2
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

def run(args):

  work_params = get_phil(args)
  C = fit_translation2(work_params)
  C.post_min_recalc()
  C.print_table()
  sum_sq = 0.
  for key in C.frame_delx.keys():
    mn_x = flex.mean(C.frame_delx[key])
    mn_y = flex.mean(C.frame_dely[key])
    print("frame %d count %4d delx %7.2f  dely %7.2f"%(key,
      len(C.frame_delx[key]),
      mn_x,
      mn_y ))
    sum_sq += mn_x*mn_x + mn_y*mn_y
  displacement = math.sqrt(sum_sq / len(C.frame_delx.keys()))
  print("rms displacement of frames %7.2f"%displacement)
  C.same_sensor_table()
  return None

if (__name__ == "__main__"):

  result = run(args=["mysql.runtag=for_may060corner","mysql.passwd=sql789",
                     "mysql.user=nick","mysql.database=xfelnks",
                     #"show_plots=True",
                     #"show_consistency=True",
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
  """
