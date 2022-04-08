# -*- coding: utf-8 -*-
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.jungfrau_metrics
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
#
from __future__ import absolute_import, division, print_function

import math
from dials.array_family import flex
from libtbx.phil import parse
import glob
import numpy as np
from skimage.transform import SimilarityTransform
from skimage.measure import ransac

help_message = '''

This program is used to analyze the agreement between observed and model
Bragg spots on the Jungfrau 16M detector. The purpose is to assess the
quality of the metrology and the spot prediction model using proximal metrics:
1) Does the overall detector position seem to change over time (by run) ?
2) Does the sensor metrology make sense (<MODEL-OBS> = 0) ?
3) Are the spots predicted correctly, in terms of RMSD(MODEL_OBS) ?
3a) Especially radial vs. transverse RMSD ?
The program should challenge dials.refine (or whatever method is used to
determine the metrology) to produce superior proximal metrics.

Example:

  cctbx.xfel.jungfrau_metrics experiment_glob.expt reflections_glob.refl

Assumptions:  Everything is in the reflections files.  Predicted and observed are
taken as given in refls["xyzcal.mm"] and refls["xyzobs.mm.value"].
No reprediction is done.
Only a single *.expt file is read, to get the detector model.
It is assumed that all the detector models are identical (not verified).
Also assume that the files adhere to a specific naming convention so as to
make the regular expressions work.  This could be generalized.
Only supports the Jungfrau 16M with 8 panels x 32 sensors

Root assumptions:  All the reflections are strong spots that have been indexed.
There are no outliers (although the RANSAC fit will remove some outliers).

Requirements:  Python 3, scikit-image
'''

# Create the phil parameters
phil_scope = parse('''
panel_numbers = True
  .type = bool
  .help = Whether to show panel numbers on each panel
verbose = False
  .type = bool
  .help = Whether to print extra statistics
input {
  refl_glob = None
    .type = str
    .optional = False
  expt_glob = None
    .type = str
    .optional = False
  predictions_glob = None
    .type = str
    .help = instead of refl_glob, give separate globs for predictions refls["xyzcalc.mm"] and observations refls["xyzobs.mm.value"]
  observations_glob = None
    .type = str
    .help = instead of refl_glob, give separate globs for predictions refls["xyzcalc.mm"] and observations refls["xyzobs.mm.value"]
}
dispatch {
  by_run = True
    .type = bool
    .help = if True, print out a per-run metric analysis. \
            Assumes the run number is encoded in the file name as done in the regex RUN below.
  by_panel_mcd_filter = True
    .type = bool
    .help = finally, repeat the sensor analysis but with panel-level 3 feature MCD filter
}
residuals {
  mcd_filter {
    enable = False
      .type = bool
      .help = on the DeltaPsi plot, apply a 3-feature MCD filter \
              FIXME refactor this so it is not repeated from detector_residuals
    mahalanobis_distance = 7
      .type = float(value_min=0.)
      .help = cutoff level expressed as a multivariate Gaussian std dev
    keep = *inliers outliers
      .type = choice
      .help = this should be obvious.  Either keep the good ones or the bad ones.
  }
}
output {
  observations_filtered_replaces = None
    .type = str
    .help = After MCD, write out the filtered observations, to a filename patterned \
            on the observations file, but replacing this string with 'filtered'.

}
include scope xfel.command_line.cspad_detector_congruence.phil_scope
''', process_includes=True)

import re
RUN = re.compile(".*run_000([0-9][0-9][0-9])")
EVENT = re.compile(".*master_([0-9][0-9][0-9][0-9][0-9])")

# Next steps

# demonstrate that re-prediction based on expt model does not change the result.
#    that way I can validate on the fly any new code to write out a modified expt file
#    that way also I can include group B

#    have an option for sorting the output table by sensor radius from beam position

# remove the need to create a "groupA" file ahead of time
#   add the ability to get this analysis with group A+B
# add radial and transverse RMSD
# based on RANSAC, output a modified expt file [then plot the new metrology XY and Z positions]

class Script(object):
  ''' Class to parse the command line options. '''

  def __init__(self):
    ''' Set the expected options. '''
    from dials.util.options import ArgumentParser
    import libtbx.load_env

    # Create the option parser
    usage = "usage: %s [options] /path/to/refined/json/file" % libtbx.env.dispatcher_name
    self.parser = ArgumentParser(
      usage=usage,
      sort_options=True,
      phil=phil_scope,
      read_experiments=True,
      read_reflections=True,
      check_format=False,
      epilog=help_message)

  def build_statistics(self,refl_gen):

    n_total = 0
    n_file = 0
    sum_delta_x = 0
    sum_delta_y = 0
    deltax = {}
    deltay = {}
    delRoverR = {}
    panel_deltax = {}
    panel_deltay = {}
    panel_deltapsi = {}
    cumPANNO = {}
    cumOBS = {}
    cumCALC = {}
    run_nos = {}

    print ("Reading %d refl files, printing the first 10"%(len(refl_gen)))
    for item in refl_gen:
      strong_refls = item["strong_refls"]
      nrefls = len(strong_refls)
      if self.params.dispatch.by_run==True:
        run_match = RUN.match(item["strfile"])
        run_token = int(run_match.group(1)); run_nos[run_token]=run_nos.get(run_token,0); run_nos[run_token]+=1
        deltax[run_token] = deltax.get(run_token, flex.double())
        deltay[run_token] = deltay.get(run_token, flex.double())
        delRoverR[run_token] = delRoverR.get(run_token, flex.double())
      OBS = flex.vec2_double()
      CALC = flex.vec2_double()

      if n_file < 10: print (item["strfile"], nrefls)
      n_total += nrefls
      n_file += 1
      fcalc = strong_refls["xyzcal.mm"]
      fobs = strong_refls["xyzobs.mm.value"]
      delpsi = strong_refls["delpsical.rad"]
      panno = strong_refls["panel"]
      if n_file==1:
        for ipanel in range(256):
          panel_deltax[ipanel] = panel_deltax.get(ipanel, flex.double())
          panel_deltay[ipanel] = panel_deltay.get(ipanel, flex.double())
          panel_deltapsi[ipanel] = panel_deltapsi.get(ipanel, flex.double())
        for isensor in range(32):
          cumOBS[isensor] = cumOBS.get(isensor, flex.vec3_double())
          cumCALC[isensor] = cumCALC.get(isensor, flex.vec3_double())
          cumPANNO[isensor] = cumPANNO.get(isensor, flex.int())
        expt_files = glob.glob(self.params.input.expt_glob)
        expt_file = expt_files[0]
        from dxtbx.model.experiment_list import ExperimentList
        int_expt = ExperimentList.from_file(expt_file, check_format=False)[0]
        D = int_expt.detector
        B = int_expt.beam
        Beam = D.hierarchy().get_beam_centre_lab((0,0,-1)) # hard code the sample to source vector for now
        distance = D.hierarchy().get_distance()
      for irefl in range(len(strong_refls)):
        sum_delta_x += (fcalc[irefl][0] - fobs[irefl][0])
        sum_delta_y += (fcalc[irefl][1] - fobs[irefl][1])
        panel = D[panno[irefl]]
        panel_deltax[panno[irefl]].append( (fcalc[irefl][0] - fobs[irefl][0]) )
        panel_deltay[panno[irefl]].append( (fcalc[irefl][1] - fobs[irefl][1]) )
        panel_deltapsi[panno[irefl]].append( delpsi[irefl] )
        CALC.append(panel.get_lab_coord(fcalc[irefl][0:2])[0:2])  # take only the xy coords; assume z along beam
        OBS.append(panel.get_lab_coord(fobs[irefl][0:2])[0:2])
        cumCALC[panno[irefl]//8].append(panel.get_lab_coord(fcalc[irefl][0:2]))
        cumOBS[panno[irefl]//8].append(panel.get_lab_coord(fobs[irefl][0:2]))
        cumPANNO[panno[irefl]//8].append(panno[irefl])
      if self.params.dispatch.by_run==True:
        for irefl in range(len(strong_refls)):
          deltax[run_token].append( (fcalc[irefl][0] - fobs[irefl][0]) )
          deltay[run_token].append( (fcalc[irefl][1] - fobs[irefl][1]) )

      R_sq = CALC.dot(CALC)
      DIFF = CALC-OBS
      Dr_over_r = (DIFF.dot(CALC))/R_sq
      if self.params.dispatch.by_run==True:
        for irefl in range(len(strong_refls)):
          delRoverR[run_token].append( Dr_over_r[irefl] )

    from libtbx import adopt_init_args
    obj = dict(n_total=n_total, n_file=n_file, sum_delta_x=sum_delta_x,
            sum_delta_y=sum_delta_y, deltax=deltax, deltay=deltay, distance=distance,
            delRoverR=delRoverR, panel_deltax=panel_deltax, panel_deltay=panel_deltay, panel_deltapsi=panel_deltapsi,
            cumOBS=cumOBS, cumCALC=cumCALC, cumPANNO=cumPANNO, O="ok")
    adopt_init_args(self, obj)
    self.ordered_run_nos = sorted(list(run_nos.keys()))
    self.run_nos = run_nos

  def build_statistics_rejecting_outliers(self,refl_gen):
    # XXX duplicates behavior of above function; should refactor to cover both cases with same code.
    n_total = 0
    n_file = 0
    sum_delta_x = 0
    sum_delta_y = 0
    deltax = {}
    deltay = {}
    delRoverR = {}
    panel_deltax = {}
    panel_deltay = {}
    panel_deltapsi = {}
    cumPANNO = {}
    cumOBS = {}
    cumCALC = {}
    run_nos = {}

    print ("Reading %d refl files, WITH REJECTION printing the first 10"%(len(refl_gen)))
    for item in refl_gen:
      strong_refls = item["strong_refls"]
      nrefls = len(strong_refls)
      keep_refls = flex.bool(nrefls, True)
      if self.params.dispatch.by_run==True:
        run_match = RUN.match(item["strfile"])
        run_token = int(run_match.group(1)); run_nos[run_token]=run_nos.get(run_token,0); run_nos[run_token]+=1
        deltax[run_token] = deltax.get(run_token, flex.double())
        deltay[run_token] = deltay.get(run_token, flex.double())
        delRoverR[run_token] = delRoverR.get(run_token, flex.double())
      OBS = flex.vec2_double()
      CALC = flex.vec2_double()

      if n_file < 1000: print (item["strfile"], nrefls)
      n_total += nrefls
      n_file += 1
      fcalc = strong_refls["xyzcal.mm"]
      fobs = strong_refls["xyzobs.mm.value"]
      delpsi = strong_refls["delpsical.rad"]
      panno = strong_refls["panel"]
      if n_file==1:
        for ipanel in range(256):
          panel_deltax[ipanel] = panel_deltax.get(ipanel, flex.double())
          panel_deltay[ipanel] = panel_deltay.get(ipanel, flex.double())
          panel_deltapsi[ipanel] = panel_deltapsi.get(ipanel, flex.double())
        for isensor in range(32):
          cumOBS[isensor] = cumOBS.get(isensor, flex.vec3_double())
          cumCALC[isensor] = cumCALC.get(isensor, flex.vec3_double())
          cumPANNO[isensor] = cumPANNO.get(isensor, flex.int())
        expt_files = glob.glob(self.params.input.expt_glob)
        expt_file = expt_files[0]
        from dxtbx.model.experiment_list import ExperimentList
        int_expt = ExperimentList.from_file(expt_file, check_format=False)[0]
        D = int_expt.detector
        B = int_expt.beam
        Beam = D.hierarchy().get_beam_centre_lab((0,0,-1)) # hard code the sample to source vector for now
        distance = D.hierarchy().get_distance()
      for irefl in range(len(strong_refls)):
        trial_delta_x = (fcalc[irefl][0] - fobs[irefl][0])
        trial_delta_y = (fcalc[irefl][1] - fobs[irefl][1])
        #if irefl==26 and nrefls==400: import ipdb;ipdb.set_trace()
        try:
          maha_distance = self.record_MCD_models[panno[irefl]].rob_cov.mahalanobis(X=[ (trial_delta_x,trial_delta_y,delpsi[irefl]) ])[0]
        except KeyError:
          maha_distance = None
        if (
            maha_distance is None or
            maha_distance >= (self.params.residuals.mcd_filter.mahalanobis_distance)**2
        ):
          keep_refls[irefl]=False
          continue
        sum_delta_x += trial_delta_x
        sum_delta_y += trial_delta_y
        panel = D[panno[irefl]]
        panel_deltax[panno[irefl]].append( trial_delta_x )
        panel_deltay[panno[irefl]].append( trial_delta_y )
        panel_deltapsi[panno[irefl]].append( delpsi[irefl] )
        CALC.append(panel.get_lab_coord(fcalc[irefl][0:2])[0:2])  # take only the xy coords; assume z along beam
        OBS.append(panel.get_lab_coord(fobs[irefl][0:2])[0:2])
        cumCALC[panno[irefl]//8].append(panel.get_lab_coord(fcalc[irefl][0:2]))
        cumOBS[panno[irefl]//8].append(panel.get_lab_coord(fobs[irefl][0:2]))
        cumPANNO[panno[irefl]//8].append(panno[irefl])
      if self.params.dispatch.by_run==True:
        for irefl in range(len(strong_refls)):
          deltax[run_token].append( (fcalc[irefl][0] - fobs[irefl][0]) )
          deltay[run_token].append( (fcalc[irefl][1] - fobs[irefl][1]) )

      R_sq = CALC.dot(CALC)
      DIFF = CALC-OBS
      Dr_over_r = (DIFF.dot(CALC))/R_sq
      if self.params.dispatch.by_run==True:
        for irefl in range(len(Dr_over_r)):
          delRoverR[run_token].append( Dr_over_r[irefl] )

      if self.params.output.observations_filtered_replaces is not None:
        filtered_file = item["obsfile"].replace(self.params.output.observations_filtered_replaces,"filtered")
        item["strong_refls"].select(keep_refls).as_file(filtered_file)

    from libtbx import adopt_init_args
    obj = dict(n_total=n_total, n_file=n_file, sum_delta_x=sum_delta_x,
            sum_delta_y=sum_delta_y, deltax=deltax, deltay=deltay, distance=distance,
            delRoverR=delRoverR, panel_deltax=panel_deltax, panel_deltay=panel_deltay, panel_deltapsi=panel_deltapsi,
            cumOBS=cumOBS, cumCALC=cumCALC, cumPANNO=cumPANNO, O="ok")
    adopt_init_args(self, obj)
    self.ordered_run_nos = sorted(list(run_nos.keys()))
    self.run_nos = run_nos

  def per_run_analysis(self):
    print ()
    print ("Total number of files and reflections for each run")
    print ("Average Delta-x,Deltay in microns, along with standard error of the mean, and RMSD(model-obs)")
    print ()
    print ("Run  n_file  n_refl  <Δx>(μm)     <Δy>(μm)  RMSD Δx(μm)  RMSD Δy(μm)  <distance*Δr/r>(μm)")

    for irun in self.ordered_run_nos:
      Sx = flex.mean_and_variance(1000.*self.deltax[irun])
      Sy = flex.mean_and_variance(1000.*self.deltay[irun])
      print(irun,
        "%6d %8d"%(self.run_nos[irun], len(self.deltax[irun])),
        "%6.2f±%.2f %6.2f±%.2f"%(1000.*flex.mean(self.deltax[irun]),Sx.unweighted_standard_error_of_mean(),
                                 1000.*flex.mean(self.deltay[irun]),Sy.unweighted_standard_error_of_mean()),
        "    %5.1f    %5.1f"%(Sx.unweighted_sample_standard_deviation(),Sy.unweighted_sample_standard_deviation()),
        "       %5.1f"%(1000. * self.distance * flex.mean(self.delRoverR[irun]))
        )
    print("All",
        "%6d %8d"%(self.n_file, self.n_total),
        "%6.2f      %6.2f     "%(1000*self.sum_delta_x/self.n_total,
                                 1000*self.sum_delta_y/self.n_total,)
        )
    print ()

  def per_sensor_analysis(self): # hardcoded Jungfrau 16M geometry
    for isensor in range(32):
      print ("Panel Sensor  <Δx>(μm)     <Δy>(μm)      Nrefl  RMS Δx(μm)  RMS Δy(μm) ")

      if len(self.cumCALC[isensor]) <= 3: continue

      for ipanel in range(8*isensor, 8*(1+isensor)):
        if len(self.panel_deltax[ipanel])<2: continue
        Sx = flex.mean_and_variance(1000.*self.panel_deltax[ipanel])
        Sy = flex.mean_and_variance(1000.*self.panel_deltay[ipanel])
        RMSDx = 1000.*math.sqrt(flex.mean(self.panel_deltax[ipanel]*self.panel_deltax[ipanel]))
        RMSDy = 1000.*math.sqrt(flex.mean(self.panel_deltay[ipanel]*self.panel_deltay[ipanel]))
        print("%3d  %3d"%(ipanel,ipanel//8),"%7.2f±%6.2f %7.2f±%6.2f %6d"%(Sx.mean(),Sx.unweighted_standard_error_of_mean(),
                                                 Sy.mean(),Sy.unweighted_standard_error_of_mean(), len(self.panel_deltax[ipanel])),
            "    %5.1f   %5.1f"%(RMSDx,RMSDy),
        )
      print("")
      cumD = (self.cumCALC[isensor]-self.cumOBS[isensor]).parts()
      print ( "All  %3d %7.2f        %7.2f        %6d"%(isensor,1000.*flex.mean(cumD[0]), 1000.*flex.mean(cumD[1]), len(cumD[0])))
      print("")

  # Now we'll do a linear least squares refinement over sensors:
  #Method 1. Simple rectilinear translation.
      if self.params.verbose:
        veclength = len(self.cumCALC[isensor])
        correction = flex.vec3_double( veclength, (flex.mean(cumD[0]), flex.mean(cumD[1]), flex.mean(cumD[2])) )

        new_delta = (self.cumCALC[isensor]-correction ) -self.cumOBS[isensor]
        for ipanel in range(8*isensor, 8*(1+isensor)):
          panel_delta = new_delta.select(self.cumPANNO[isensor]==ipanel)
          if len(panel_delta)<2: continue
          deltax_part, deltay_part = panel_delta.parts()[0:2]
          RMSDx = 1000.*math.sqrt( flex.mean(deltax_part * deltax_part) )
          RMSDy = 1000.*math.sqrt( flex.mean(deltay_part * deltay_part) )
          Sx = flex.mean_and_variance(1000.*deltax_part)
          Sy = flex.mean_and_variance(1000.*deltay_part)
          print("%3d  %3d"%(ipanel,ipanel//8),"%7.2f±%6.2f %7.2f±%6.2f %6d"%(Sx.mean(),Sx.unweighted_standard_error_of_mean(),
                                                 Sy.mean(),Sy.unweighted_standard_error_of_mean(), len(deltax_part)),
          "    %5.1f   %5.1f"%(RMSDx,RMSDy),
          )
        print()
  # Method 2. Translation + rotation.
      src = []
      dst = []
      for icoord in range(len(self.cumCALC[isensor])):
        src.append(self.cumCALC[isensor][icoord][0:2])
        dst.append(self.cumOBS[isensor][icoord][0:2])
      src = np.array(src)
      dst = np.array(dst)

      # estimate affine transform model using all coordinates
      model = SimilarityTransform()
      model.estimate(src, dst)

      # robustly estimate affine transform model with RANSAC
      model_robust, inliers = ransac((src, dst), SimilarityTransform, min_samples=3,
                               residual_threshold=2, max_trials=10)
      outliers = flex.bool(inliers == False)

      # compare "true" and estimated transform parameters
      if self.params.verbose:
        print("Similarity transform:")
        print("%2d"%isensor, "Scale: %.5f,"%(model.scale),
        "Translation(μm): (%7.2f,"%(1000.*model.translation[0]),
        "%7.2f),"%(1000.*model.translation[1]),
        "Rotation (°): %7.4f"%((180./math.pi)*model.rotation))
      print("RANSAC:")
      print("%2d"%isensor, "Scale: %.5f,"%(model_robust.scale),
      "Translation(μm): (%7.2f,"%(1000.*model_robust.translation[0]),
      "%7.2f),"%(1000.*model_robust.translation[1]),
      "Rotation (°): %7.4f,"%((180./math.pi)*model_robust.rotation),
      "Outliers:",outliers.count(True)
      )
      """from documentation:
      X = a0 * x - b0 * y + a1 = s * x * cos(rotation) - s * y * sin(rotation) + a1
      Y = b0 * x + a0 * y + b1 = s * x * sin(rotation) + s * y * cos(rotation) + b1"""

      oldCALC = self.cumCALC[isensor].parts()

      ransacCALC = flex.vec3_double(
               (float(model_robust.scale) * oldCALC[0] * math.cos(model_robust.rotation) -
               float(model_robust.scale) * oldCALC[1] * math.sin(model_robust.rotation) +
               float(model_robust.translation[0])),
               (float(model_robust.scale) * oldCALC[0] * math.sin(model_robust.rotation) +
               float(model_robust.scale) * oldCALC[1] * math.cos(model_robust.rotation) +
               float(model_robust.translation[1])),
               oldCALC[2]
               )
      new_delta = ransacCALC - self.cumOBS[isensor]
      inlier_delta = new_delta.select(~outliers)
      inlier_panno = self.cumPANNO[isensor].select(~outliers)

      for ipanel in range(8*isensor, 8*(1+isensor)):
        panel_delta = inlier_delta.select(inlier_panno==ipanel)
        if len(panel_delta)<2: continue
        deltax_part, deltay_part = panel_delta.parts()[0:2]
        RMSDx = 1000.*math.sqrt( flex.mean(deltax_part * deltax_part) )
        RMSDy = 1000.*math.sqrt( flex.mean(deltay_part * deltay_part) )
        Sx = flex.mean_and_variance(1000.*deltax_part)
        Sy = flex.mean_and_variance(1000.*deltay_part)
        print("%3d  %3d"%(ipanel,ipanel//8),"%7.2f±%6.2f %7.2f±%6.2f %6d"%(Sx.mean(),Sx.unweighted_standard_error_of_mean(),
                                                 Sy.mean(),Sy.unweighted_standard_error_of_mean(), len(deltax_part)),
        "    %5.1f   %5.1f"%(RMSDx,RMSDy),
        )

      if self.params.verbose:
        print("")
        cumD = (inlier_delta).parts()
        print ( "     %3d %7.2f        %7.2f        %6d\n"%(isensor,1000.*flex.mean(cumD[0]), 1000.*flex.mean(cumD[1]), len(cumD[0])))
      print("----\n")

  def per_panel_mcd_filter(self):
    from xfel.metrology.panel_fitting import Panel_MCD_Filter
    self.record_MCD_models = {}
    for isensor in range(32):
      if len(self.cumCALC[isensor]) <= 3: continue
      for ipanel in range(8*isensor, 8*(1+isensor)):
        if len(self.panel_deltax[ipanel])<2: continue
        MCD = Panel_MCD_Filter(lab_coords_x = self.panel_deltax[ipanel],
                               lab_coords_y = self.panel_deltay[ipanel],
                               data = self.panel_deltapsi[ipanel],
                               i_panel=ipanel, delta_scalar=1., params=self.params, verbose=True)
        sX,sY,sPsi = MCD.scatter_coords()
        self.record_MCD_models[ipanel] = MCD

  def build_reflections_generator(self):
    help = """the application accepts two mutually exclusive input formats
    1) indexed strong spots and corresponding predictions from input.refl_glob
    2) strong from input.observations_glob, predictions from input.predictions_glob"""
    assert (self.params.input.refl_glob is not None).__xor__(
           (self.params.input.observations_glob is not None) and
            self.params.input.predictions_glob is not None
           ), help
    if self.params.input.refl_glob is not None:
      class all_refl_list:
        def __init__(O): O._all_file_list = glob.glob(self.params.input.refl_glob)
        def __len__(O): return len(O._all_file_list)
        def __iter__(O):
          for item in O._all_file_list:
            yield dict(strfile = item, obsfile = item, strong_refls=flex.reflection_table.from_file(item))
      return all_refl_list()

    else:
      class all_refl_list:
        def __init__(O):
          O.pred_list = glob.glob(self.params.input.predictions_glob)
          O.obs_list = glob.glob(self.params.input.observations_glob)
          O.obs_reverse_lookup = {}
          for obs in O.obs_list:
              run_match = RUN.match(obs)
              run_token = int(run_match.group(1))
              event_match = EVENT.match(obs)
              event_token = int(event_match.group(1))
              O.obs_reverse_lookup[(run_token,event_token)] = obs
        def __len__(O): return len(O.pred_list)
        def __iter__(O):
          from cctbx.miller import match_indices
          n_total = 0
          n_file = 0
          n_reindex = 0
          n_matches = 0
          n_strong_no_integration = 0
          n_integrated = 0
          n_common = 0
          n_weak = 0

          for item in O.pred_list:
            run_match = RUN.match(item)
            run_token = int(run_match.group(1))
            event_match = EVENT.match(item)
            event_token = int(event_match.group(1))
            obs = O.obs_reverse_lookup[(run_token,event_token)]
            strong_refls = flex.reflection_table.from_file(obs)
            nrefls = len(strong_refls)
            ri = reindex_miller = strong_refls["miller_index"]
            select_indexed = reindex_miller != (0,0,0)
            reindex = (select_indexed).count(True)
            strong_and_indexed = strong_refls.select(select_indexed)
            #print (obs, nrefls, reindex)
            n_total += nrefls
            n_file += 1
            n_reindex += reindex

            int_refls = flex.reflection_table.from_file(item)
            ii = integration_miller = int_refls["miller_index"]
            MM = match_indices(strong_and_indexed["miller_index"], int_refls["miller_index"])
            #print("  Strong+indexed",reindex, "integrated",len(ii), "in common",len(MM.pairs()),
            #"indexed but no integration", len(MM.singles(0)), "integrated weak", len(MM.singles(1)))
            n_integrated += len(ii)
            n_common += len(MM.pairs())
            n_strong_no_integration += len(MM.singles(0))
            n_weak += len(MM.singles(1))

            P = MM.pairs()
            A = P.column(0)
            B = P.column(1)
            strong_and_indexed = strong_and_indexed.select(A)
            int_refls = int_refls.select(B)
            # transfer over the calculated positions from integrate2 to strong refls
            strong_and_indexed["xyzcal.mm"] = int_refls["xyzcal.mm"]
            strong_and_indexed["xyzcal.px"] = int_refls["xyzcal.px"]
            strong_and_indexed["delpsical.rad"] = int_refls["delpsical.rad"]

            yield dict(strfile = item, obsfile=obs, strong_refls=strong_and_indexed)
          print ("Grand total is %d from %d files of which %d reindexed"%(n_total,n_file, n_reindex))
          print ("TOT Strong+indexed",n_reindex, "integrated",n_integrated, "in common",
          n_common, "indexed but no integration", n_strong_no_integration, "integrated weak", n_weak)
      return all_refl_list()

  def run(self):
    ''' Parse the options. '''
    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    assert params.input.expt_glob is not None, "Must give a input.expt_glob= string"
    self.params = params

    generate_refls = self.build_reflections_generator()
    assert generate_refls is not None

    self.build_statistics(refl_gen = generate_refls)
    print("Finished reading total refl files",self.n_file)

    if self.params.dispatch.by_run==True:
      self.per_run_analysis()
    self.per_sensor_analysis()

    if self.params.dispatch.by_panel_mcd_filter:
      self.per_panel_mcd_filter()
      self.build_statistics_rejecting_outliers(refl_gen = self.build_reflections_generator())
      if self.params.dispatch.by_run==True:
        self.per_run_analysis()
      self.per_sensor_analysis()

if __name__ == '__main__':
  script = Script()
  script.run()

