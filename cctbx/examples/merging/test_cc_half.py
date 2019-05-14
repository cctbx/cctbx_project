from __future__ import absolute_import, division, print_function
import os,sys
from cctbx.examples.merging import test_levenberg_sparse as test
import libtbx.load_env
import math

# test script assumes that you get the data files directly from the author (NKS) and
# install them in the directory "xscale_reserve" at the same dir-level as cctbx_project

def get_model_intensities(work_params,ordered_intensities):
  ordered_intensities.set_observation_type_xray_intensity()
  ordered_intensities.show_summary()
  mtz_out = ordered_intensities.as_mtz_dataset(column_root_label="Iobs",title="PSI simulation",wavelength=None)
  mtz_obj = mtz_out.mtz_object()
  mtz_obj.write(work_params.scaling.mtz_file)

if __name__=="__main__":
  modules_dist = os.path.abspath(os.path.join(libtbx.env.dist_path("cctbx"),"../.."))
  datadir = os.path.join(modules_dist,"xscale_reserve") # Get files directly from author, NKS
  plot_flag=False
  esd_plot_flag=False
  written_files = []

  #for N in [25, 200, 300, 400, 500, 800, 1000, 2000, 5000]:
  for N in [400]:
    #for trans in [1.0, 0.1, 0.01, 0.001, 0.0001, 0.00001]:
    for trans in [0.001,]:
      tagbase = "PSIfulls"
      tag = "%s_N%04d_tE%1d"%(tagbase,N,-math.log10(trans))
      for half_data_flag in [1,2,0]:
        case = test.execute_case(datadir, n_frame=N, transmittance=trans, apply_noise=True,
               plot=plot_flag, esd_plot = esd_plot_flag,half_data_flag = half_data_flag)
        model_subset = case.ordered_intensities[0:len(case.Fit_I)]
        from cctbx.xray import observation_types
        fitted_miller_array = model_subset.customized_copy(data = case.Fit_I, sigmas = case.Fit_I_stddev,
          observation_type=observation_types.intensity())
        output_result = fitted_miller_array.select(case.I_visited==1)
        outfile = "%s_s%1d_levmar.mtz"%(tag,half_data_flag)
        output_result.show_summary(prefix="%s: "%outfile)
        mtz_out = output_result.as_mtz_dataset(column_root_label="Iobs",title=outfile,wavelength=None)
        mtz_obj = mtz_out.mtz_object()
        mtz_obj.write(outfile)
        written_files.append(outfile)
        print("OK")
        #raw_input("OK")
        sys.stdout.flush()

      from xfel.command_line.cxi_merge import master_phil
      import iotbx.phil
      args = ["output.prefix=%s"%tag,"scaling.algorithm=levmar",
              "d_min=%f"%(case.d_min),"output.n_bins=%d"%(case.n_bins),"model=%s"%(os.path.join(datadir,"not_used.pdb")),
              "mtz_file=%s"%("rigged_filename.mtz"),
              "mtz_column_F=iobs",
              "merge_anomalous=True","scaling.show_plots=False","log_cutoff=-20"]
      phil = iotbx.phil.process_command_line(args=args, master_string=master_phil)# .show()
      work_params = phil.work.extract()
      from xfel.merging.phil_validation import application
      application(work_params)
      get_model_intensities(work_params,case.ordered_intensities)
      from xfel.cxi.cxi_cc import run_cc
      run_cc(work_params,"h,k,l",sys.stdout)
