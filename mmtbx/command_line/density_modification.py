from cctbx.array_family import flex
from cctbx import crystal
from iotbx.option_parser import option_parser
import iotbx.phil
from iotbx.reflection_file_reader import any_reflection_file
import mmtbx.utils
from libtbx.utils import show_times_at_exit

from mmtbx import density_modification
import mmtbx.utils

import os, sys

master_params_including_IO_str = """\
density_modification {
  input {
    reflection_data {
      %s
    }
    experimental_phases {
      %s
    }
    unit_cell = None
      .type = unit_cell
      .optional = False
      .style = bold noauto
    space_group = None
      .type = space_group
      .optional = False
      .style = bold noauto
  }
  output {
    map {
      file_name = None
        .type = path
        .help = The file name for the final density-modified map
        .short_caption = Output map file
        .style = bold noauto
      format = xplor *ccp4
        .type = choice
        .short_caption = Map format
      scale = *sigma volume
        .type = choice(multi=False)
        .short_caption = Map scaling
        .expert_level=2
    }
    map_coefficients {
      file_name = none
        .type = path
        .help = The file name for the coefficients of the final density-modified map
        .short_caption = Output map coefficients
      format = *mtz cns
        .optional=True
        .type = choice
    }
  }
%s
}
""" %(mmtbx.utils.data_and_flags_str,
      mmtbx.utils.experimental_phases_params_str,
      density_modification.master_params_str)

def defaults(log):
  parsed = iotbx.phil.parse(
    master_params_including_IO_str, process_includes=True)
  print >> log
  return parsed

def run(args, log = sys.stdout):
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
  crystal_symmetry = crystal.symmetry(
    unit_cell=params.input.unit_cell,
    space_group_info=params.input.space_group)
  reflection_files = []
  for rfn in (params.input.reflection_data.file_name,
              params.input.experimental_phases.file_name):
    if os.path.isfile(str(rfn)):
      reflection_files.append(iotbx.reflection_file_reader.any_reflection_file(
        file_name=rfn, ensure_read_access=False))
  server = iotbx.reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    reflection_files=reflection_files)
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

  fo = fo.map_to_asu()
  hl_coeffs = hl_coeffs.map_to_asu()

  fo = fo.eliminate_sys_absent().average_bijvoet_mates()
  hl_coeffs = hl_coeffs.eliminate_sys_absent().average_bijvoet_mates()
  fo, hl_coeffs = fo.common_sets(hl_coeffs)

  model_map = None
  if len(processed_args.pdb_file_names):
    pdb_file = mmtbx.utils.pdb_file(
      pdb_file_names=processed_args.pdb_file_names)
    xs = pdb_file.pdb_inp.xray_structure_simple()
    if params.change_basis_to_niggli_cell:
      change_of_basis_op = xs.change_of_basis_op_to_niggli_cell()
      xs = xs.change_basis(change_of_basis_op)
    fmodel = mmtbx.f_model.manager(f_obs=fo.change_basis(change_of_basis_op),
                                   abcd=hl_coeffs.change_basis(change_of_basis_op),
                                   xray_structure=xs)
    true_phases = fmodel.f_model().phases()
    model_map = fmodel.electron_density_map().fft_map(
      resolution_factor=params.grid_resolution_factor,
      map_type="2mFo-DFc").real_map_unpadded()

  dm = density_modify(fo, hl_coeffs, params, model_map=model_map)

  map_coeffs = dm.map_coeffs_in_original_setting

  # output map if requested
  map_params = params.output.map
  if map_params.file_name is not None:
    fft_map = map_coeffs.fft_map(resolution_factor=params.grid_resolution_factor)
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
  map_coeff_params = params.output.map_coefficients
  if map_coeff_params.file_name is not None:
    # XXX TO DO: write out map coeffs!
    pass


# just for development purposes, compare the correlation of the
# density-modified map with map calculated from the model at each cycle
class density_modify(density_modification.density_modification):

  def __init__(self, fo, hl_coeffs, params, model_map=None):
    self.model_map = model_map
    self.correlation_coeffs = flex.double()
    density_modification.density_modification.__init__(
      self, fo, hl_coeffs, params)
    if len(self.correlation_coeffs) > 1:
      fft_map = self.map_coeffs_start.fft_map(
        resolution_factor=self.params.grid_resolution_factor
      ).apply_sigma_scaling()
      corr = flex.linear_correlation(
        self.model_map.as_1d(), fft_map.real_map_unpadded().as_1d())
      print "Starting dm/model correlation: %.6f" %corr.coefficient()
      print "Final dm/model correlation:    %.6f" %self.correlation_coeffs[-1]
      fft_map.as_ccp4_map(file_name="starting.map", labels=[])

  def compute_map(self):
    density_modification.density_modification.compute_map(self)
    if self.model_map is not None:
      print
      corr = flex.linear_correlation(self.model_map.as_1d(), self.map.as_1d())
      print "dm/model correlation:"
      corr.show_summary()
      self.correlation_coeffs.append(corr.coefficient())

if __name__ == '__main__':
  show_times_at_exit()
  run(sys.argv[1:])
