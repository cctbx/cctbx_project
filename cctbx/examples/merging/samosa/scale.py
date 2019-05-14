from __future__ import absolute_import, division, print_function
from six.moves import range
import os,sys
import libtbx
import iotbx.phil
from libtbx.development.timers import Timer
from cctbx.array_family import flex
from cctbx.examples.merging.task4 import prepare_observations_for_scaling
from cctbx.examples.merging.data_utilities import I_and_G_base_estimate, plot_it, show_correlation
from cctbx.examples.merging.data_utilities import show_histogram
from cctbx.examples.merging.data_subset import mapper_factory
from cctbx import miller
from xfel.merging.database.merging_database_flex import read_experiments
from cctbx.examples.merging.test_levenberg_sparse import xscale6e

class execute_case(object):
 def __init__(self,datadir,work_params,plot=False,esd_plot=False,half_data_flag=0):
  casetag = work_params.output.prefix
  # read the ground truth values back in
  from six.moves import cPickle as pickle
  # it is assumed (for now) that the reference millers contain a complete asymmetric unit
  # of indices, within the (d_max,d_min) region of interest and possibly outside the region.
  reference_millers = pickle.load(open(os.path.join(datadir,casetag+"_miller.pickle"),"rb"))
  experiment_manager = read_experiments(work_params)

  obs = pickle.load(open(os.path.join(datadir,casetag+"_observation.pickle"),"rb"))
  print("Read in %d observations"%(len(obs["observed_intensity"])))
  reference_millers.show_summary(prefix="Miller index file ")

  print(len(obs["frame_lookup"]),len(obs["observed_intensity"]), flex.max(obs['miller_lookup']),flex.max(obs['frame_lookup']))
  max_frameno = flex.max(obs["frame_lookup"])

  from iotbx import mtz
  mtz_object = mtz.object(file_name=work_params.scaling.mtz_file)
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

  print("original unique",len(reference_millers.indices()))
  print("isomorphous set",len(I_sim.indices()))
  print("pairs",len(matches.pairs()))
  iso_data = flex.double(len(reference_millers.indices()))

  for pair in matches.pairs():
    iso_data[pair[0]] = I_sim.data()[pair[1]]

  reference_data = miller.array(miller_set = reference_millers,
                                data = iso_data)
  reference_data.set_observation_type_xray_intensity()

  FOBS = prepare_observations_for_scaling(work_params,obs=obs,
                                          reference_intensities=reference_data,
                                          files = experiment_manager.get_files(),
                                          half_data_flag=half_data_flag)

  I,I_visited,G,G_visited = I_and_G_base_estimate(FOBS,params=work_params)
  print("I length",len(I), "G length",len(G), "(Reference set; entire asymmetric unit)")
  assert len(reference_data.data()) == len(I)

  #presumably these assertions fail when half data are taken for CC1/2 or d_min is cut
  model_I = reference_data.data()[0:len(I)]

  T = Timer("%d frames"%(len(G), ))

  mapper = mapper_factory(xscale6e)
  minimizer = mapper(I,G,I_visited,G_visited,FOBS,params=work_params,
                     experiments=experiment_manager.get_experiments())

  del T
  minimizer.show_summary()

  Fit = minimizer.e_unpack()
  Gstats=flex.mean_and_variance(Fit["G"].select(G_visited==1))
  print("G mean and standard deviation:",Gstats.mean(),Gstats.unweighted_sample_standard_deviation())
  if "Bfactor" in work_params.levmar.parameter_flags:
    Bstats=flex.mean_and_variance(Fit["B"].select(G_visited==1))
    print("B mean and standard deviation:",Bstats.mean(),Bstats.unweighted_sample_standard_deviation())
  show_correlation(Fit["I"],model_I,I_visited,"Correlation of I:")
  Fit_stddev = minimizer.e_unpack_stddev()

  # XXX FIXME known bug:  the length of Fit["G"] could be smaller than the length of experiment_manager.get_files()
  # Not sure if this has any operational drawbacks.  It's a result of half-dataset selection.

  if plot:
    plot_it(Fit["I"], model_I, mode="I")
    if "Rxy" in work_params.levmar.parameter_flags:
      show_histogram(Fit["Ax"],"Histogram of x rotation (degrees)")
      show_histogram(Fit["Ay"],"Histogram of y rotation (degrees)")
  print()

  if esd_plot:
    minimizer.esd_plot()

  from cctbx.examples.merging.show_results import show_overall_observations
  table1,self.n_bins,self.d_min = show_overall_observations(
           Fit["I"],Fit_stddev["I"],I_visited,
           reference_data,FOBS,title="Statistics for all reflections",
           work_params = work_params)

  self.FSIM=FOBS
  self.ordered_intensities=reference_data
  self.reference_millers=reference_millers
  self.Fit_I=Fit["I"]
  self.Fit_I_stddev=Fit_stddev["I"]
  self.I_visited=I_visited
  self.Fit = Fit
  self.experiments = experiment_manager

def run(show_plots,args):
  from xfel.command_line.cxi_merge import master_phil
  phil = iotbx.phil.process_command_line(args=args, master_string=master_phil).show()
  work_params = phil.work.extract()
  from xfel.merging.phil_validation import application,samosa
  application(work_params)
  samosa(work_params)
  if ("--help" in args):
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
      print("OK s%1d"%half_data_flag)
      #raw_input("OK?")

    """Guest code to retrieve the modified orientations after rotational fitting is done"""
    if "Rxy" in work_params.levmar.parameter_flags:
      all_A = [e.crystal.get_A() for e in case.experiments.get_experiments()]
      all_files = case.experiments.get_files()
      all_x = case.Fit["Ax"]
      all_y = case.Fit["Ay"]

      from scitbx import matrix
      x_axis = matrix.col((1.,0.,0.))
      y_axis = matrix.col((0.,1.,0.))
      out = open("aaaaa","w")
      for x in range(len(all_A)):
        Rx = x_axis.axis_and_angle_as_r3_rotation_matrix(angle=all_x[x], deg=True)
        Ry = y_axis.axis_and_angle_as_r3_rotation_matrix(angle=all_y[x], deg=True)
        modified_A = Rx * Ry * all_A[x]
        filename = all_files[x]
        print(filename, " ".join([str(a) for a in modified_A.elems]), file=out)



    work_params.scaling.algorithm="levmar"
    from xfel.cxi.cxi_cc import run_cc
    run_cc(work_params,work_params.model_reindex_op,sys.stdout)

  else:
    execute_case(datadir, work_params, plot=show_plots)
