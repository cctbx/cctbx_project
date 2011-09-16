from mmtbx import density_modification
import mmtbx.utils
from cctbx.array_family import flex
from cctbx import crystal
from iotbx.option_parser import option_parser
import iotbx.phil
from iotbx.reflection_file_reader import any_reflection_file
from libtbx.utils import show_times_at_exit, Sorry
from libtbx import runtime_utils
from libtbx import easy_pickle
from libtbx import adopt_init_args
import mmtbx.maps
from mmtbx.ncs import ncs
from scitbx.math import nearest_phase, phase_error
import os, sys

master_params_including_IO_str = """\
density_modification {
  input
    .style = noauto
  {
    reflection_data {
      %s
    }
    experimental_phases {
      %s
    }
    map_coefficients
      .optional=True
      .help = Optional starting map coefficients
    {
      %s
    }
    ncs_file_name = None
      .type = path
      .optional = True
      .short_caption = NCS file
      .help = .ncs_spec file produced by phenix.find_ncs or phenix.find_ncs_from_density
    unit_cell = None
      .type = unit_cell
      .optional = False
      .style = bold noauto
    space_group = None
      .type = space_group
      .optional = False
      .style = bold noauto
  }
  output
    .style = noauto
  {
    map {
      file_name = None
        .type = path
        .help = The file name for the final density-modified map
        .short_caption = Output map file
        .style = bold noauto new_file file_type:mtz
      format = xplor *ccp4
        .type = choice
        .short_caption = Map format
      scale = *sigma volume
        .type = choice(multi=False)
        .short_caption = Map scaling
        .expert_level=2
    }
    mtz {
      file_name = None
        .type = path
        .help = The file name for the coefficients of the final density-modified map
        .short_caption = Output map coefficients
        .style = bold noauto new_file
      output_hendrickson_lattman_coefficients = False
        .type = bool
        .help = Output density modified phase probability distributions
    }
  }
  job_title = None
    .type = str
    .input_size = 400
    .help = Job title in PHENIX GUI, not used on command line
    .style = noauto
%s
}
""" %(mmtbx.utils.data_and_flags_str,
      mmtbx.utils.experimental_phases_params_str,
      mmtbx.utils.map_coefficents_params_str,
      density_modification.master_params_str)

def defaults(log):
  parsed = iotbx.phil.parse(
    master_params_including_IO_str, process_includes=True)
  print >> log
  return parsed

def run(args, log = sys.stdout, as_gui_program=False):
  if(len(args)==0):
    parsed = defaults(log=log)
    parsed.show(prefix="  ", out=log)
    return
  command_line = (option_parser()
                  .enable_symmetry_comprehensive()
                  .option("-q", "--quiet",
                          action="store_true",
                          default=False,
                          help="suppress output")
                  .option("--output_plots",
                          action="store_true",
                          default=False)
                  ).process(args=args)
  parsed = defaults(log=log)
  processed_args = mmtbx.utils.process_command_line_args(
    args=command_line.args,
    cmd_cs=command_line.symmetry,
    master_params=parsed,
    log=sys.stdout,
    suppress_symmetry_related_errors=True)
  processed_args.params.show()
  params = processed_args.params.extract().density_modification
  output_plots = command_line.options.output_plots

  crystal_symmetry = crystal.symmetry(
    unit_cell=params.input.unit_cell,
    space_group_info=params.input.space_group)
  reflection_files = {}
  for rfn in (params.input.reflection_data.file_name,
              params.input.experimental_phases.file_name,
              params.input.map_coefficients.file_name):
    if os.path.isfile(str(rfn)) and rfn not in reflection_files:
      reflection_files.setdefault(
        rfn, iotbx.reflection_file_reader.any_reflection_file(
          file_name=rfn, ensure_read_access=False))
  server = iotbx.reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    reflection_files=reflection_files.values())
  fo = mmtbx.utils.determine_data_and_flags(
    server,
    parameters=params.input.reflection_data,
    extract_r_free_flags=False).f_obs
  hl_coeffs = mmtbx.utils.determine_experimental_phases(
    server,
    params.input.experimental_phases,
    log=sys.stdout,
    parameter_scope="",
    working_point_group=None,
    symmetry_safety_check=True,
    ignore_all_zeros=True)
  if params.input.map_coefficients.file_name is not None:
    map_coeffs = server.get_phases_deg(
      file_name=params.input.map_coefficients.file_name,
      labels=params.input.map_coefficients.labels,
      convert_to_phases_if_necessary=False,
      original_phase_units=None,
      parameter_scope="",
      parameter_name="labels").map_to_asu()
  else:
    map_coeffs = None
  ncs_object = None
  if params.input.ncs_file_name is not None:
    ncs_object = ncs.ncs()
    ncs_object.read_ncs(params.input.ncs_file_name)
    ncs_object.display_all(log=log)

  fo = fo.map_to_asu()
  hl_coeffs = hl_coeffs.map_to_asu()

  fo = fo.eliminate_sys_absent().average_bijvoet_mates()
  hl_coeffs = hl_coeffs.eliminate_sys_absent().average_bijvoet_mates()

  model_map = None
  model_map_coeffs = None
  if len(processed_args.pdb_file_names):
    pdb_file = mmtbx.utils.pdb_file(
      pdb_file_names=processed_args.pdb_file_names)
    xs = pdb_file.pdb_inp.xray_structure_simple()
    fo_, hl_ = fo, hl_coeffs
    if params.change_basis_to_niggli_cell:
      change_of_basis_op = xs.change_of_basis_op_to_niggli_cell()
      xs = xs.change_basis(change_of_basis_op)
      fo_ = fo_.change_basis(change_of_basis_op).map_to_asu()
      hl_ = hl_.change_basis(change_of_basis_op).map_to_asu()
    #fo_, hl_ = fo_.common_sets(hl_)
    fmodel_refined = mmtbx.utils.fmodel_simple(
      f_obs=fo_,
      scattering_table="wk1995",#XXX pva: 1) neutrons? 2) move up as a parameter.
      xray_structures=[xs],
      bulk_solvent_correction=True,
      anisotropic_scaling=True,
      r_free_flags=fo_.array(data=flex.bool(fo_.size(), False)))
    fmodel_refined.update(abcd=hl_)

    master_phil = mmtbx.maps.map_and_map_coeff_master_params()
    map_params = master_phil.fetch(iotbx.phil.parse("""\
map_coefficients {
  map_type = 2mFo-DFc
  isotropize = True
}
""")).extract().map_coefficients[0]
    model_map_coeffs = mmtbx.maps.map_coefficients_from_fmodel(
      fmodel_refined, map_params)
    model_map = model_map_coeffs.fft_map(
      resolution_factor=params.grid_resolution_factor).real_map_unpadded()

  import time

  t0 = time.time()
  dm = density_modify(
    params,
    fo,
    hl_coeffs,
    ncs_object=ncs_object,
    map_coeffs=map_coeffs,
    model_map_coeffs=model_map_coeffs,
    log=log,
    as_gui_program=as_gui_program)
  time_dm = time.time()-t0
  print >> log, "Time taken for density modification: %.2fs" %time_dm

  # run cns
  if 0:
    from cctbx.development import cns_density_modification
    cns_result = cns_density_modification.run(params, fo, hl_coeffs)
    print cns_result.modified_map.all()
    print dm.map.all()
    dm_map_coeffs = dm.map_coeffs_in_original_setting
    from cctbx import maptbx, miller
    crystal_gridding = maptbx.crystal_gridding(
      dm_map_coeffs.unit_cell(),
      space_group_info=dm_map_coeffs.space_group().info(),
      pre_determined_n_real=cns_result.modified_map.all())
    dm_map = miller.fft_map(crystal_gridding, dm_map_coeffs).apply_sigma_scaling()
    corr = flex.linear_correlation(cns_result.modified_map.as_1d(), dm_map.real_map_unpadded().as_1d())
    print "CNS dm/mmtbx dm correlation:"
    corr.show_summary()
    if dm.model_map_coeffs is not None:
      model_map = miller.fft_map(
        crystal_gridding,
        dm.miller_array_in_original_setting(dm.model_map_coeffs)).apply_sigma_scaling()
      corr = flex.linear_correlation(cns_result.modified_map.as_1d(), model_map.real_map_unpadded().as_1d())
      print "CNS dm/model correlation:"
      corr.show_summary()

  if output_plots:
    plots_to_make = (
      "fom", "skewness",
      "r1_factor", "r1_factor_fom", "mean_solvent_density", "mean_protein_density",
      "f000_over_v", "k_flip", "rms_solvent_density", "rms_protein_density",
      "standard_deviation_local_rms", "mean_delta_phi", "mean_delta_phi_initial",
      )
    from matplotlib.backends.backend_pdf import PdfPages
    from libtbx import pyplot

    stats = dm.get_stats()
    pdf = PdfPages("density_modification.pdf")

    if len(dm.correlation_coeffs) > 1:
      if 0:
        start_coeffs, model_coeffs = dm.map_coeffs_start.common_sets(model_map_coeffs)
        model_phases = model_coeffs.phases(deg=True).data()
        exptl_phases = nearest_phase(
          model_phases, start_coeffs.phases(deg=True).data(), deg=True)
        corr = flex.linear_correlation(exptl_phases, model_phases)
        corr.show_summary()
        fig = pyplot.figure()
        ax = fig.add_subplot(1,1,1)
        ax.set_title("phases start")
        ax.set_xlabel("Experimental phases")
        ax.set_ylabel("Phases from refined model")
        ax.scatter(exptl_phases,
                   model_phases,
                   marker="x", s=10)
        pdf.savefig(fig)
        #
        dm_coeffs, model_coeffs = dm.map_coeffs.common_sets(model_map_coeffs)
        model_phases = model_coeffs.phases(deg=True).data()
        dm_phases = nearest_phase(
          model_phases, dm_coeffs.phases(deg=True).data(), deg=True)
        corr = flex.linear_correlation(dm_phases, model_phases)
        corr.show_summary()
        fig = pyplot.figure()
        ax = fig.add_subplot(1,1,1)
        ax.set_title("phases dm")
        ax.set_xlabel("Phases from density modification")
        ax.set_ylabel("Phases from refined model")
        ax.scatter(dm_phases,
                   model_phases,
                   marker="x", s=10)
        pdf.savefig(fig)
      #
      data = dm.correlation_coeffs
      fig = pyplot.figure()
      ax = fig.add_subplot(1,1,1)
      ax.set_title("correlation coefficient")
      ax.plot(range(1, dm.i_cycle+2), data)
      pdf.savefig(fig)
      #
      data = dm.mean_phase_errors
      fig = pyplot.figure()
      ax = fig.add_subplot(1,1,1)
      ax.set_title("Mean effective phase errors")
      ax.plot(range(1, dm.i_cycle+2), data)
      pdf.savefig(fig)

    for plot in plots_to_make:
      data = [getattr(stats.get_cycle_stats(i), plot) for i in range(1, dm.i_cycle+2)]
      fig = pyplot.figure()
      ax = fig.add_subplot(1,1,1)
      ax.set_title(plot.replace("_", " "))
      ax.plot(range(1, dm.i_cycle+2), data)
      pdf.savefig(fig)

    data = [stats.get_cycle_stats(i).rms_solvent_density/
            stats.get_cycle_stats(i).rms_protein_density
            for i in range(1, dm.i_cycle+2)]
    fig = pyplot.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_title("RMS solvent/protein density ratio")
    ax.plot(range(1, dm.i_cycle+2), data)
    pdf.savefig(fig)

    pdf.close()

  dm_map_coeffs = dm.map_coeffs_in_original_setting
  dm_hl_coeffs = dm.hl_coeffs_in_original_setting

  # output map if requested
  map_params = params.output.map
  if map_params.file_name is not None:
    fft_map = dm_map_coeffs.fft_map(resolution_factor=params.grid_resolution_factor)
    if map_params.scale == "sigma":
      fft_map.apply_sigma_scaling()
    else:
      fft_map.apply_volume_scaling()
    gridding_first = gridding_last = None
    title_lines = []
    if map_params.format == "xplor":
      fft_map.as_xplor_map(
        file_name      = map_params.file_name,
        title_lines    = title_lines,
        gridding_first = gridding_first,
        gridding_last  = gridding_last)
    else :
      fft_map.as_ccp4_map(
        file_name      = map_params.file_name,
        gridding_first = gridding_first,
        gridding_last  = gridding_last,
        labels=title_lines)

  # output map coefficients if requested
  mtz_params = params.output.mtz
  if mtz_params.file_name is not None:
    label_decorator=iotbx.mtz.ccp4_label_decorator()
    fo = dm.miller_array_in_original_setting(dm.f_obs_complete).common_set(dm_map_coeffs)
    mtz_dataset = fo.as_mtz_dataset(
      column_root_label="F",
      label_decorator=label_decorator)
    mtz_dataset.add_miller_array(
      dm_map_coeffs,
      column_root_label="FWT",
      label_decorator=label_decorator)
    phase_source = dm.miller_array_in_original_setting(dm.phase_source).common_set(dm_map_coeffs)
    mtz_dataset.add_miller_array(
      phase_source.array(data=flex.abs(phase_source.data())),
      column_root_label="FOM",
      label_decorator=label_decorator)
    mtz_dataset.add_miller_array(
      phase_source.array(data=phase_source.phases(deg=True).data()),
      column_root_label="PHIB",
      label_decorator=None)
    if mtz_params.output_hendrickson_lattman_coefficients:
      mtz_dataset.add_miller_array(
        dm_hl_coeffs,
        column_root_label="P",
        label_decorator=label_decorator)
    mtz_dataset.mtz_object().write(mtz_params.file_name)

  return result(
    map_file=map_params.file_name,
    mtz_file=mtz_params.file_name,
    stats=dm.get_stats())

# just for development purposes, compare the correlation of the
# density-modified map with map calculated from the model at each cycle
class density_modify(density_modification.density_modification):

  def __init__(self, params,
                     fo,
                     hl_coeffs,
                     ncs_object=None,
                     map_coeffs=None,
                     model_map_coeffs=None,
                     log=None,
                     as_gui_program=False):
    self.model_map_coeffs = model_map_coeffs
    self.correlation_coeffs = flex.double()
    self.mean_phase_errors = flex.double()
    density_modification.density_modification.__init__(
      self, params, fo, hl_coeffs,
      ncs_object=ncs_object,
      map_coeffs=map_coeffs,
      as_gui_program=as_gui_program)
    if len(self.correlation_coeffs) > 1:
      model_coeffs, start_coeffs = self.model_map_coeffs.common_sets(self.map_coeffs_start)
      model_fft_map = model_coeffs.fft_map(
        resolution_factor=self.params.grid_resolution_factor).apply_sigma_scaling()
      fft_map = start_coeffs.fft_map(
        resolution_factor=self.params.grid_resolution_factor
      ).apply_sigma_scaling()
      corr = flex.linear_correlation(
        model_fft_map.real_map_unpadded().as_1d(), fft_map.real_map_unpadded().as_1d())
      print "Starting dm/model correlation: %.6f" %corr.coefficient()
      print "Final dm/model correlation:    %.6f" %self.correlation_coeffs[-1]
      fft_map.as_ccp4_map(file_name="starting.map", labels=[])

  def compute_map(self):
    density_modification.density_modification.compute_map(self)
    if self.model_map_coeffs is not None:
      model_coeffs, dm_coeffs = self.model_map_coeffs.common_sets(self.map_coeffs)
      fft_map = model_coeffs.fft_map(
        resolution_factor=self.params.grid_resolution_factor).apply_sigma_scaling()
      dm_map = dm_coeffs.fft_map(
        resolution_factor=self.params.grid_resolution_factor).apply_sigma_scaling()
      print
      corr = flex.linear_correlation(
        fft_map.real_map_unpadded().as_1d(), dm_map.real_map_unpadded().as_1d())
      print "dm/model correlation:"
      corr.show_summary()
      self.correlation_coeffs.append(corr.coefficient())
      self.mean_phase_errors.append(flex.mean(phase_error(
        flex.arg(model_coeffs.data()),
        flex.arg(dm_coeffs.data())))/density_modification.pi_180)

def validate_params (params) :
  params_ = params.density_modification
  if (params_.input.reflection_data.file_name is None) :
    raise Sorry("No reflection data provided.")
  if (params_.input.reflection_data.labels is None) :
    raise Sorry("Data labels not specified.")
  if (params_.input.experimental_phases.file_name is None) :
    raise Sorry("Experimental phases (Hendrickson-Lattman coefficients " +
                "not specified.")
  if (params_.input.experimental_phases.labels is None) :
    raise Sorry("Experimental phase labels not specified.")
  if ((params_.output.map.file_name is None) and
      (params_.output.mtz.file_name is None)) :
    raise Sorry("No output requested!")
  if (params_.solvent_fraction is None) :
    raise Sorry("Please specify the solvent fraction!")

class launcher (runtime_utils.simple_target) :
  def __call__ (self) :
    result = run(args=list(self.args),
                 log=sys.stdout,
                 as_gui_program=True)
    config_file = self.args[0]
    pkl_file = config_file[:-4] + ".pkl"
    easy_pickle.dump(pkl_file, result)
    return result

class result (object) :
  def __init__ (self, map_file, mtz_file, stats) :
    adopt_init_args(self, locals())

  def extract_loggraph (self) :
    return self.stats.extract_loggraph()

  def get_final_job_statistics (self) :
    stats = [
      ("FOM", self.stats.get_cycle_stats(-1).fom),
      ("Skewness", self.stats.get_cycle_stats(-1).skewness)
    ]
    return stats

  def finish_job (self) :
    output_files = []
    if (self.mtz_file is not None) :
      output_files.append((self.mtz_file, "Map coefficients"))
    if (self.map_file is not None) :
      output_files.append((self.map_file, "Real-space map"))
    stats = self.get_final_job_statistics()
    return (output_files, stats)

if __name__ == '__main__':
  show_times_at_exit()
  run(sys.argv[1:])
