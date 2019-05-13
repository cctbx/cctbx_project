from __future__ import division, absolute_import
from __future__ import print_function
import clipper
from iotbx import reflection_file_converter
from cctbx.array_family import flex
import scitbx.math
import math
import sys

def exercise_Compute_phifom_from_abcd_interface(args):
  miller_arrays = reflection_file_converter.run(
    args=args, simply_return_all_miller_arrays=True)
  complex_input = None
  for miller_array in miller_arrays:
    if (miller_array.is_complex_array()):
      complex_input = miller_array
      break
  if (complex_input is not None):
    print("complex_input.info():", complex_input.info())
  fom_input = None
  for miller_array in miller_arrays:
    if (miller_array.is_real_array()
        and miller_array.info().lower().find("fom")):
      fom_input = miller_array
      break
  if (fom_input is not None):
    print("fom_input.info():", fom_input.info())
  for miller_array in miller_arrays:
    if (isinstance(miller_array.data(), flex.hendrickson_lattman)):
      print("Hendrickson-Lattman coefficients:", miller_array.info())
      n_bad_figures_of_merit = 0
      n_bad_phases = 0
      if (miller_array.anomalous_flag()):
        print("Clipper cannot currently handle anomalous arrays" \
            + " with Hendrickson-Lattman coefficients.")
        continue
      phi_fom = clipper.Compute_phifom_from_abcd_interface(
        unit_cell=miller_array.unit_cell(),
        space_group=miller_array.space_group(),
        miller_indices=miller_array.indices(),
        phase_probabilities=miller_array.data())
      phi_clipper = phi_fom.centroid_phases()*180/math.pi
      fom_clipper = phi_fom.figures_of_merit()
      fom_deltas = fom_input.data() - fom_clipper
      if (fom_input is not None):
        perm = flex.sort_permutation(flex.abs(fom_deltas), reverse=True)
        fom_deltas_sorted = fom_deltas.select(perm)
        fom_clipper_sorted = fom_clipper.select(perm)
        fom_input_sorted = fom_input.data().select(perm)
        pp_sorted = miller_array.data().select(perm)
        indices_sorted = miller_array.indices().select(perm)
        centric_flags_sorted = fom_input.centric_flags().data().select(perm)
        print("FOM comparison")
        for f,i,c,d,p,ix in zip(centric_flags_sorted,
                                fom_input_sorted,
                                fom_clipper_sorted,
                                fom_deltas_sorted,
                                pp_sorted,
                                indices_sorted):
          print(f, "%.3f %.3f %.3f" % (i,c,d), end=' ')
          if (abs(d) > 0.01):
            print("LOOK", "%.2f %.2f %.2f %.2f" % p, ix, end=' ')
            n_bad_figures_of_merit += 1
          print()
        print()
      if (complex_input is not None):
        phi_input = complex_input.phases(deg=True).data()
        phi_deltas = flex.double()
        for phi1,phi2 in zip(phi_input, phi_clipper):
          phi_deltas.append(
            scitbx.math.signed_phase_error(phi1, phi2, deg=True))
        perm = flex.sort_permutation(flex.abs(phi_deltas), reverse=True)
        phi_deltas_sorted = phi_deltas.select(perm)
        phi_clipper_sorted = phi_clipper.select(perm)
        phi_input_sorted = phi_input.select(perm)
        pp_sorted = miller_array.data().select(perm)
        indices_sorted = miller_array.indices().select(perm)
        centric_flags_sorted =complex_input.centric_flags().data().select(perm)
        fom_clipper_sorted = fom_clipper.select(perm)
        print("PHI comparison")
        for f,i,c,d,m,p,ix in zip(centric_flags_sorted,
                                phi_input_sorted,
                                phi_clipper_sorted,
                                phi_deltas_sorted,
                                fom_clipper_sorted,
                                pp_sorted,
                                indices_sorted):
          if (m < 0.01): continue
          print(f, "%.3f %.3f %.3f" % (i,c,d), end=' ')
          if (abs(d) > 3 and (max(p) < 1000 or abs(d) > 10)):
            print("LOOK", "%.2f %.2f %.2f %.2f" % p, "fom=%.3f" % m, ix, end=' ')
            n_bad_phases += 1
          print()
        print()
      assert n_bad_figures_of_merit == 0
      assert n_bad_phases == 0

def run():
  exercise_Compute_phifom_from_abcd_interface(sys.argv[1:])
  print("OK")

if (__name__ == "__main__"):
  run()
