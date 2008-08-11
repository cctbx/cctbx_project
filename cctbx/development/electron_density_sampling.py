from cctbx import xray
from cctbx import miller
from cctbx.command_line.structure_factor_timings import read_structure
import sys

def timings(structure, wing_cutoff=1.e-3):
  print "wing_cutoff for following fft calculations: %3.1e"%wing_cutoff
  for calc_type,exp_table_one_over_step_size in (("exp function:",0),
                                                 ("exp table:",-100)):
    print calc_type
    for d_min in [4,3,2,1]:
      structure_ng = structure.deep_copy_scatterers()
      structure_ng.scattering_type_registry(d_min=d_min, table="n_gaussian")
      structure_4g = structure.deep_copy_scatterers()
      structure_4g.scattering_type_registry(table="it1992")
      miller_set = miller.build_set(
        crystal_symmetry=structure,
        d_min=d_min,
        anomalous_flag=False)
      miller_set.show_summary()
      times = []
      for structure in (structure_ng, structure_4g):
        structure.scattering_type_registry().show_summary()
        f_calc_object = xray.structure_factors.from_scatterers(
        miller_set=miller_set,
        wing_cutoff=wing_cutoff,
        exp_table_one_over_step_size=exp_table_one_over_step_size)(
          xray_structure=structure,
          miller_set=miller_set,
          algorithm="fft")
        times.append(f_calc_object.manager().estimate_time_fft.time_sampling)
        print "  %.2f seconds," % times[-1]
      print "d_min=%d: %.2f s / %.2f s" % (d_min, times[0], times[1]),
      if (times[1] != 0):
        print "= %.2f" % (times[0] / times[1]),
      print
      sys.stdout.flush()
  print

def run(args):
  assert len(args) == 1
  structure = read_structure(args[0])
  structure.show_summary()
  print
  timings(structure=structure)

if (__name__ == "__main__"):
  run(sys.argv[1:])
