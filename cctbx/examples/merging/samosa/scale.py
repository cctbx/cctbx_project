from __future__ import division
import os,sys
import libtbx
import iotbx.phil
from libtbx.development.timers import Timer
from cctbx.array_family import flex
from cctbx.examples.merging.task4 import prepare_observations_for_scaling
from cctbx.examples.merging.data_utilities import I_and_G_base_estimate, plot_it, show_correlation
from cctbx.examples.merging.data_subset import mapper_factory
from cctbx import miller

from cctbx.examples.merging.test_levenberg_sparse import xscale6e

class execute_case(object):
 def __init__(self,datadir,work_params,plot=False,esd_plot=False,half_data_flag=0):
  casetag = work_params.output.prefix
  # read the ground truth values back in
  import cPickle as pickle
  reference_millers = pickle.load(open(os.path.join(datadir,casetag+"_miller.pickle"),"rb"))

  obs = pickle.load(open(os.path.join(datadir,casetag+"_observation.pickle"),"rb"))
  print "Read in %d observations"%(len(obs["observed_intensity"]))
  reference_millers.show_summary(prefix="Miller index file ")

  print len(obs["frame_lookup"]),len(obs["observed_intensity"]), flex.max(obs['miller_lookup']),flex.max(obs['frame_lookup'])
  max_frameno = flex.max(obs["frame_lookup"])

  from iotbx import mtz
  mtz_object = mtz.object(file_name="1jb0_Fright.mtz")
  #for array in mtz_object.as_miller_arrays():
  #  this_label = array.info().label_string()
  #  print this_label, array.observation_type()
  I_sim = mtz_object.as_miller_arrays()[0].as_intensity_array()
  I_sim.show_summary()
  MODEL_REINDEX_OP = work_params.model_reindex_op
  I_sim = I_sim.change_basis(MODEL_REINDEX_OP).map_to_asu()

  #match up isomorphous (the simulated fake F's) with experimental unique set
  matches = miller.match_multi_indices(
      miller_indices_unique=reference_millers.indices(),
      miller_indices=I_sim.indices())

  print "original unique",len(reference_millers.indices())
  print "isomorphous set",len(I_sim.indices())
  print "pairs",len(matches.pairs())
  iso_data = flex.double(len(reference_millers.indices()))

  for pair in matches.pairs():
    iso_data[pair[0]] = I_sim.data()[pair[1]]

  reference_data = miller.array(miller_set = reference_millers,
                                data = iso_data)
  reference_data.set_observation_type_xray_intensity()

  FOBS = prepare_observations_for_scaling(work_params,obs=obs,
                                          reference_intensities=reference_data,
                                          half_data_flag=half_data_flag)

  I,I_visited,G,G_visited = I_and_G_base_estimate(FOBS)
  print "I length",len(I), "G length",len(G)
  assert len(reference_data.data()) == len(I)

  #presumably these assertions fail when half data are taken for CC1/2 or d_min is cut
  model_I = reference_data.data()[0:len(I)]

  T = Timer("%d frames"%(len(G), ))

  mapper = mapper_factory(xscale6e)
  minimizer = mapper(I,G,I_visited,G_visited,FOBS)

  del T
  minimizer.show_summary()

  Fit_I, Fit_G, Fit_B = minimizer.e_unpack()
  Gstats=flex.mean_and_variance(Fit_G.select(G_visited==1))
  Bstats=flex.mean_and_variance(Fit_B.select(G_visited==1))
  print "G mean and standard deviation:",Gstats.mean(),Gstats.unweighted_sample_standard_deviation()
  print "B mean and standard deviation:",Bstats.mean(),Bstats.unweighted_sample_standard_deviation()
  show_correlation(Fit_I,model_I,I_visited,"Correlation of I:")
  Fit_I_stddev, Fit_G_stddev, Fit_B_stddev = minimizer.e_unpack_stddev()

  if plot:
    plot_it(Fit_I, model_I, mode="I")
  print

  if esd_plot:
    minimizer.esd_plot()

  from cctbx.examples.merging.show_results import show_overall_observations
  table1,self.n_bins,self.d_min = show_overall_observations(
           Fit_I,Fit_I_stddev,I_visited,
           reference_data,FOBS,title="Statistics for all reflections",
           work_params = work_params)

  self.FSIM=FOBS
  self.ordered_intensities=reference_data
  self.reference_millers=reference_millers
  self.Fit_I=Fit_I
  self.Fit_I_stddev=Fit_I_stddev
  self.I_visited=I_visited

def run(show_plots,args):
  from xfel.command_line.cxi_merge import master_phil
  phil = iotbx.phil.process_command_line(args=args, master_string=master_phil).show()
  work_params = phil.work.extract()
  from xfel.merging.phil_validation import application,samosa
  application(work_params)
  samosa(work_params)
  if ("--help" in args) :
    libtbx.phil.parse(master_phil.show())
    return

  datadir = "."
  written_files = []
  if work_params.levmar.compute_cc_half:

    for half_data_flag in [1,2,0]:
      case = execute_case(datadir, work_params, plot=show_plots, half_data_flag=half_data_flag)
      assert len(case.Fit_I)==len(case.ordered_intensities.indices())==len(case.reference_millers.indices())
      model_subset = case.reference_millers[0:len(case.Fit_I)]
      fitted_miller_array = miller.array (miller_set = model_subset,
                                        data = case.Fit_I, sigmas = case.Fit_I_stddev)
      fitted_miller_array.set_observation_type_xray_intensity()
      output_result = fitted_miller_array.select(case.I_visited==1)
      outfile = "%s_s%1d_levmar.mtz"%(work_params.output.prefix,half_data_flag)
      output_result.show_summary(prefix="%s: "%outfile)
      mtz_out = output_result.as_mtz_dataset(column_root_label="Iobs",title=outfile,wavelength=None)
      mtz_obj = mtz_out.mtz_object()
      mtz_obj.write(outfile)
      written_files.append(outfile)
      print "OK s%1d"%half_data_flag
      #raw_input("OK?")

    work_params.scaling.algorithm="levmar"
    print case.d_min,work_params.d_min
    # XXX FIXME. Put a filter on input data such that the execute_case() function respects work_params.d_min
    # XXX FIXME. Regardless of the tabular output, the levmar analysis is done on the unfiltered set.
    #assert case.d_min == work_params.d_min

    from xfel.cxi.cxi_cc import run_cc
    run_cc(work_params,work_params.model_reindex_op,sys.stdout)

  else:
    execute_case(datadir, work_params, plot=show_plots)

