# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME dev.cxi.merge
#
# $Id$
from __future__ import absolute_import, division, print_function

from xfel.command_line.cxi_merge import master_phil
from libtbx.utils import Usage, multi_out
import sys,os
from xfel.command_line.cxi_merge import get_observations, scaling_manager, show_overall_observations, scaling_result
from libtbx import easy_pickle

help_message = '''
Script for merging xfel data
'''

class Script(object):
  '''A class for running the script.'''

  def __init__(self, scaler_class):
    # The script usage
    import libtbx.load_env
    self.usage = "usage: %s [options] [param.phil] " % libtbx.env.dispatcher_name
    self.parser = None
    self.scaler_class = scaler_class

  def initialize(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from iotbx.phil import parse
    phil_scope = parse(master_phil)
    # Create the parser
    self.parser = OptionParser(
      usage=self.usage,
      phil=phil_scope,
      epilog=help_message)
    self.parser.add_option(
        '--plots',
        action='store_true',
        default=False,
        dest='show_plots',
        help='Show some plots.')

    # Parse the command line. quick_parse is required for MPI compatibility
    params, options = self.parser.parse_args(show_diff_phil=True,quick_parse=True)
    self.params = params
    self.options = options

  def validate(self):
    from xfel.merging.phil_validation import application
    application(self.params)
    if ((self.params.d_min is None) or
        (self.params.data is None) or
        ( (self.params.model is None) and self.params.scaling.algorithm != "mark1") ) :
      command_name = os.environ["LIBTBX_DISPATCHER_NAME"]
      raise Usage(command_name + " "
                  "d_min=4.0 "
                  "data=~/scratch/r0220/006/strong/ "
                  "model=3bz1_3bz2_core.pdb")
    if ((self.params.rescale_with_average_cell) and
        (not self.params.set_average_unit_cell)) :
      raise Usage("If rescale_with_average_cell=True, you must also specify "+
        "set_average_unit_cell=True.")
    if [self.params.raw_data.sdfac_auto, self.params.raw_data.sdfac_refine, self.params.raw_data.errors_from_sample_residuals].count(True) > 1:
      raise Usage("Specify only one of sdfac_auto, sdfac_refine or errors_from_sample_residuals.")

  def read_models(self):
    # Read Nat's reference model from an MTZ file.  XXX The observation
    # type is given as F, not I--should they be squared?  Check with Nat!
    log = open("%s.log" % self.params.output.prefix, "w")
    out = multi_out()
    out.register("log", log, atexit_send_to=None)
    out.register("stdout", sys.stdout)
    print("I model", file=out)
    if self.params.model is not None:
      from xfel.merging.general_fcalc import run as run_fmodel
      i_model = run_fmodel(self.params)
      self.params.target_unit_cell = i_model.unit_cell()
      self.params.target_space_group = i_model.space_group_info()
      i_model.show_summary()
    else:
      i_model = None

    print("Target unit cell and space group:", file=out)
    print("  ", self.params.target_unit_cell, file=out)
    print("  ", self.params.target_space_group, file=out)
    from xfel.command_line.cxi_merge import consistent_set_and_model
    self.miller_set, self.i_model = consistent_set_and_model(self.params,i_model)
    self.frame_files = get_observations(self.params)
    self.out = out

  def scale_all(self):
    scaler = self.scaler_class(
      miller_set=self.miller_set,
      i_model=self.i_model,
      params=self.params,
      log=self.out)
    scaler.scale_all(self.frame_files)
    return scaler

  def finalize(self, scaler):
    scaler.show_unit_cell_histograms()
    if (self.params.rescale_with_average_cell) :
      average_cell_abc = scaler.uc_values.get_average_cell_dimensions()
      average_cell = uctbx.unit_cell(list(average_cell_abc) +
        list(self.params.target_unit_cell.parameters()[3:]))
      self.params.target_unit_cell = average_cell
      print("", file=out)
      print("#" * 80, file=out)
      print("RESCALING WITH NEW TARGET CELL", file=out)
      print("  average cell: %g %g %g %g %g %g" % \
        self.params.target_unit_cell.parameters(), file=out)
      print("", file=out)
      scaler.reset()
      scaler.scale_all(frame_files)
      scaler.show_unit_cell_histograms()
    if False : #(self.params.output.show_plots) :
      try :
        plot_overall_completeness(completeness)
      except Exception as e :
        print("ERROR: can't show plots")
        print("  %s" % str(e))
    print("\n", file=self.out)

    sum_I, sum_I_SIGI = scaler.sum_intensities()

    miller_set_avg = self.miller_set.customized_copy(
      unit_cell=self.params.target_unit_cell)
    table1 = show_overall_observations(
      obs=miller_set_avg,
      redundancy=scaler.completeness,
      redundancy_to_edge=scaler.completeness_predictions,
      summed_wt_I=scaler.summed_wt_I,
      summed_weight=scaler.summed_weight,
      ISIGI=scaler.ISIGI,
      n_bins=self.params.output.n_bins,
      title="Statistics for all reflections",
      out=self.out,
      work_params=self.params)
    print("", file=self.out)
    if self.params.model is not None:
      n_refl, corr = scaler.get_overall_correlation(sum_I)
    else:
      n_refl, corr = ((scaler.completeness > 0).count(True), 0)
    print("\n", file=self.out)
    table2 = show_overall_observations(
      obs=miller_set_avg,
      redundancy=scaler.summed_N,
      redundancy_to_edge=scaler.completeness_predictions,
      summed_wt_I=scaler.summed_wt_I,
      summed_weight=scaler.summed_weight,
      ISIGI=scaler.ISIGI,
      n_bins=self.params.output.n_bins,
      title="Statistics for reflections where I > 0",
      out=self.out,
      work_params=self.params)
    #from libtbx import easy_pickle
    #easy_pickle.dump(file_name="stats.pickle", obj=stats)
    #stats.report(plot=self.params.plot)
    #miller_counts = miller_set_p1.array(data=stats.counts.as_double()).select(
    #  stats.counts != 0)
    #miller_counts.as_mtz_dataset(column_root_label="NOBS").mtz_object().write(
    #  file_name="nobs.mtz")
    if self.params.data_subsubsets.subsubset is not None and self.params.data_subsubsets.subsubset_total is not None:
      easy_pickle.dump("scaler_%d.pickle"%self.params.data_subsubsets.subsubset, scaler)
    explanation = """
  Explanation:
  Completeness       = # unique Miller indices present in data / # Miller indices theoretical in asymmetric unit
  Asu. Multiplicity  = # measurements / # Miller indices theoretical in asymmetric unit
  Obs. Multiplicity  = # measurements / # unique Miller indices present in data
  Pred. Multiplicity = # predictions on all accepted images / # Miller indices theoretical in asymmetric unit"""
    print(explanation, file=self.out)
    mtz_file, miller_array = scaler.finalize_and_save_data()
    #table_pickle_file = "%s_graphs.pkl" % self.params.output.prefix
    #easy_pickle.dump(table_pickle_file, [table1, table2])
    loggraph_file = os.path.abspath("%s_graphs.log" % self.params.output.prefix)
    f = open(loggraph_file, "w")
    f.write(table1.format_loggraph())
    f.write("\n")
    f.write(table2.format_loggraph())
    f.close()
    result = scaling_result(
      miller_array=miller_array,
      plots=scaler.get_plot_statistics(),
      mtz_file=mtz_file,
      loggraph_file=loggraph_file,
      obs_table=table1,
      all_obs_table=table2,
      n_reflections=n_refl,
      overall_correlation=corr)
    easy_pickle.dump("%s.pkl" % self.params.output.prefix, result)
    return result

  def show_plot(self,result):
    if result is not None:
      if (self.options.show_plots) :
        try :
          result.plots.show_all_pyplot()
          from wxtbx.command_line import loggraph
          loggraph.run([result.loggraph_file])
        except Exception as e :
          print("Can't display plots")
          print("You should be able to view them by running this command:")
          print("  wxtbx.loggraph %s" % result.loggraph_file)
          raise e

  def run(self):
    '''Execute the script.'''
    self.initialize()
    self.validate()
    self.read_models()
    scaler = self.scale_all()
    result = self.finalize(scaler)
    return result

if __name__ == '__main__':
  script = Script(scaling_manager)
  result = script.run()
  script.show_plot(result)
  print("DONE")

