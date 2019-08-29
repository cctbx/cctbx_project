
from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry
import sys

def run(args, out=sys.stdout):
  from mmtbx.refinement import anomalous_scatterer_groups
  import mmtbx.command_line
  master_phil = mmtbx.command_line.generate_master_phil_with_inputs(
    phil_string="""
map_type = *anom_residual llg
  .type = choice
exclude_waters = False
  .type = bool
exclude_non_water_light_elements = True
  .type = bool
n_cycles_max=None
  .type = int
map_sigma_min = 3.0
  .type = float
wavelength = None
  .type = float
refine = *f_prime *f_double_prime
  .type = choice(multi=True)
reset_water_u_iso = True
  .type = bool
""",
    enable_automatic_twin_detection=True)
  usage_string = """\
mmtbx.refine_anomalous_substructure model.pdb data.mtz [options]

Iterative identification of anomalously scattering atoms in the anomalous
residual map (simple or Phaser LLG), followed by refinement of the anomalous
scattering coefficients.  Intended as a diagnostic/development tool only!
"""
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=master_phil,
    out=out,
    process_pdb_file=False,
    prefer_anomalous=True,
    usage_string=usage_string)
  fmodel = cmdline.fmodel
  if (not fmodel.f_obs().anomalous_flag()):
    raise Sorry("Anomalous data required.")
  pdb_hierarchy = cmdline.pdb_hierarchy
  params = cmdline.params
  return anomalous_scatterer_groups.refine_anomalous_substructure(
    fmodel=fmodel,
    pdb_hierarchy=pdb_hierarchy,
    wavelength=params.wavelength,
    map_type=params.map_type,
    exclude_waters=params.exclude_waters,
    exclude_non_water_light_elements=params.exclude_non_water_light_elements,
    n_cycles_max=params.n_cycles_max,
    map_sigma_min=params.map_sigma_min,
    refine=params.refine,
    reset_water_u_iso=params.reset_water_u_iso,
    out=out)

if (__name__ == "__main__"):
  run(sys.argv[1:])
