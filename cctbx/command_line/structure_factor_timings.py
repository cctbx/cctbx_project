from __future__ import absolute_import, division, print_function
from cctbx import xray
from cctbx import eltbx
from cctbx import miller
from cctbx.array_family import flex
import iotbx.pdb
from libtbx.utils import user_plus_sys_time
from libtbx import easy_pickle
import sys

def timings(structure, d_min, fft_only=False,
            wing_cutoff_reference=1.e-6,
            wing_cutoff_others=1.e-3):
  structure_ng = structure.deep_copy_scatterers()
  structure_5g = structure.deep_copy_scatterers()
  structure_4g = structure.deep_copy_scatterers()
  structure_2g = structure.deep_copy_scatterers()
  structure_1g = structure.deep_copy_scatterers()
  structure_ng.scattering_type_registry(d_min=d_min, table="n_gaussian")
  structure_5g.scattering_type_registry(table="wk1995")
  structure_4g.scattering_type_registry(table="it1992")
  custom_dict = {}
  custom_dict.update(eltbx.xray_scattering.two_gaussian_agarwal_isaacs.table)
  structure_2g.scattering_type_registry(custom_dict=custom_dict)
  custom_dict.update(eltbx.xray_scattering.one_gaussian_agarwal_1978.table)
  structure_1g.scattering_type_registry(custom_dict=custom_dict)
  miller_set = miller.build_set(
    crystal_symmetry=structure,
    d_min=d_min,
    anomalous_flag=False)
  miller_set.show_summary()
  print("d_min:", d_min)
  if (fft_only):
    timer = user_plus_sys_time()
    f_calc_reference = xray.structure_factors.from_scatterers(
      miller_set=miller_set,
      wing_cutoff=wing_cutoff_reference,
      exp_table_one_over_step_size=0)(
        xray_structure=structure,
        miller_set=miller_set,
        algorithm="fft").f_calc().data()
    print("fft exp function wing_cutoff=%3.1e: %.2f seconds" % (
      wing_cutoff_reference, timer.elapsed()))
  else:
    timer = user_plus_sys_time()
    f_calc_reference = xray.structure_factors_simple(
      structure_5g.unit_cell(),
      structure_5g.space_group(),
      miller_set.indices(),
      structure_5g.scatterers(),
      structure_5g.scattering_type_registry()).f_calc()
    print("direct simple: %.2f seconds" % timer.elapsed())
  f_calc_reference = flex.abs(f_calc_reference)
  print("wing_cutoff for following fft calculations: %3.1e"%wing_cutoff_others)
  for structure in (structure_ng, structure_5g, structure_4g,
                    structure_2g, structure_1g):
    structure.scattering_type_registry().show_summary()
    if (not fft_only):
      for calc_type,cos_sin_flag in (("direct cos+sin function:",False),
                                     ("direct cos+sin table:",True)):
        timer = user_plus_sys_time()
        f_calc = miller_set.structure_factors_from_scatterers(
          xray_structure=structure,
          algorithm="direct",
          cos_sin_table=cos_sin_flag).f_calc()
        print("  %-24s %.2f seconds," % (calc_type, timer.elapsed()), end=' ')
        ls = xray.targets_least_squares_residual(
          f_calc_reference, f_calc.data(), False, 1)
        print("r-factor: %.6f" % ls.target())
    for calc_type,exp_table_one_over_step_size in (("fft exp function:",0),
                                                   ("fft exp table:",-100)):
      timer = user_plus_sys_time()
      f_calc = xray.structure_factors.from_scatterers(
        miller_set=miller_set,
        wing_cutoff=wing_cutoff_others,
        exp_table_one_over_step_size=exp_table_one_over_step_size)(
          xray_structure=structure,
          miller_set=miller_set,
          algorithm="fft").f_calc()
      print("  %-24s %.2f seconds," % (calc_type, timer.elapsed()), end=' ')
      ls = xray.targets_least_squares_residual(
        f_calc_reference, f_calc.data(), False, 1)
      print("r-factor: %.6f" % ls.target())
  print()

def read_structure(file_name):
  try: return easy_pickle.load(file_name)
  except KeyboardInterrupt: raise
  except Exception: pass
  try:
    if file_name.endswith('.res') or file_name.endswith('.ins'):
      return xray.structure.from_shelx(filename=file_name)
    else:
      return iotbx.pdb.input(file_name=file_name).xray_structure_simple()
  except KeyboardInterrupt: raise
  except Exception: pass
  raise RuntimeError("Unknown file format: %s" % file_name)

def run(args):
  "usage: cctbx.structure_factor_timings" \
  + " [--fft_only] coordinate_file [d-spacing ...]"
  try: i = args.index("--fft_only")
  except Exception: fft_only = False
  else:
    fft_only = True
    args = args[:]
    del args[i]
  if (len(args) == 0):
    print(run.__doc__, file=sys.stderr)
    return
  file_name = args[0]
  try:
    d_spacings = [float(arg) for arg in args[1:]]
  except Exception:
    print(run.__doc__, file=sys.stderr)
    return
  if (d_spacings == []):
    d_spacings = [4,3,2,1]
  structure = read_structure(file_name)
  structure.show_summary()
  print()
  for d_min in d_spacings:
    timings(structure=structure, d_min=d_min, fft_only=fft_only)

if (__name__ == "__main__"):
  run(sys.argv[1:])
