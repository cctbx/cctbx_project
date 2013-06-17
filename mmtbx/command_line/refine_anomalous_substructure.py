
from __future__ import division
from libtbx.utils import Usage, Sorry
import libtbx.phil
import sys

def run (args, out=sys.stdout) :
  from mmtbx.refinement import anomalous_scatterer_groups
  import mmtbx.utils
  master_phil = libtbx.phil.parse("""
include scope mmtbx.utils.cmdline_input_phil_str
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
""", process_includes=True)
  if (len(args) == 0) or ("--help" in args) :
    raise Usage("""\
mmtbx.refine_anomalous_substructure model.pdb data.mtz [options]

Iterative identification of anomalously scattering atoms in the anomalous
residual map (simple or Phaser LLG), followed by refinement of the anomalous
scattering coefficients.  Intended as a diagnostic/development tool only!

Full options:
%s""" % master_phil.as_str(prefix="  "))
  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
    update_f_part1_for="refinement",
    args=args,
    master_phil=master_phil,
    out=out,
    process_pdb_file=False,
    prefer_anomalous=True)
  fmodel = cmdline.fmodel
  if (not fmodel.f_obs().anomalous_flag()) :
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

if (__name__ == "__main__") :
  run(sys.argv[1:])
