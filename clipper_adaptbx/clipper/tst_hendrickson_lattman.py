import clipper
from iotbx import reflection_file_converter
from cctbx import utils
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from stdlib import math
import sys

def exercise_Compute_phifom_from_abcd_interface(args):
  miller_arrays = reflection_file_converter.run(
    args=args, simply_return_all_miller_arrays=0001)
  complex_array = None
  for miller_array in miller_arrays:
    if (miller_array.is_complex()):
      complex_array = miller_array
      break
  if (complex_array is not None):
    print "complex_array.info():", complex_array.info()
  fom_array = None
  for miller_array in miller_arrays:
    if (miller_array.is_real() and miller_array.info().lower().find("fom")):
      fom_array = miller_array
      break
  if (fom_array is not None):
    print "fom_array.info():", fom_array.info()
  for miller_array in miller_arrays:
    if (isinstance(miller_array.data(), flex.hendrickson_lattman)):
      print "Hendrickson-Lattman coefficients:", miller_array.info()
      if (miller_array.anomalous_flag()):
        print "Clipper cannot currently handle anomalous arrays" \
            + " with Hendrickson-Lattman coefficients."
        continue
      phi_fom = clipper.Compute_phifom_from_abcd_interface(
        unit_cell=miller_array.unit_cell(),
        space_group=miller_array.space_group(),
        miller_indices=miller_array.indices(),
        phase_probabilities=miller_array.data())
      phi_clipper = phi_fom.centroid_phases()*180/math.pi
      fom_clipper = phi_fom.figures_of_merit()
      fom_deltas = fom_array.data() - fom_clipper
      if (fom_array is not None):
        perm = flex.sort_permutation(flex.abs(fom_deltas), reverse=0001)
        fom_deltas_sorted = fom_deltas.select(perm)
        fom_clipper_sorted = fom_clipper.select(perm)
        fom_input_sorted = fom_array.data().select(perm)
        pp_sorted = miller_array.data().select(perm)
        indices_sorted = miller_array.indices().select(perm)
        centric_flags_sorted = fom_array.centric_flags().data().select(perm)
        print "FOM comparison"
        for f,i,c,d,p,ix in zip(centric_flags_sorted,
                                fom_input_sorted,
                                fom_clipper_sorted,
                                fom_deltas_sorted,
                                pp_sorted,
                                indices_sorted):
          print f, "%.3f %.3f %.3f" % (i,c,d),
          if (abs(d) > 0.01): print "LOOK", "%.2f %.2f %.2f %.2f" % p, ix,
          print
        print
      if (complex_array is not None):
        phi_input = complex_array.phases(deg=0001).data()
        phi_deltas = flex.double()
        for phi1,phi2 in zip(phi_input, phi_clipper):
          phi_deltas.append(utils.signed_phase_error(phi1, phi2, deg=0001))
        perm = flex.sort_permutation(flex.abs(phi_deltas), reverse=0001)
        phi_deltas_sorted = phi_deltas.select(perm)
        phi_clipper_sorted = phi_clipper.select(perm)
        phi_input_sorted = phi_input.select(perm)
        pp_sorted = miller_array.data().select(perm)
        indices_sorted = miller_array.indices().select(perm)
        centric_flags_sorted =complex_array.centric_flags().data().select(perm)
        print "PHI comparison"
        for f,i,c,d,p,ix in zip(centric_flags_sorted,
                                phi_input_sorted,
                                phi_clipper_sorted,
                                phi_deltas_sorted,
                                pp_sorted,
                                indices_sorted):
          if (approx_equal(p, [0,0,0,0])): continue
          print f, "%.3f %.3f %.3f" % (i,c,d),
          if (abs(d) > 3): print "LOOK", "%.2f %.2f %.2f %.2f" % p, ix,
          print
        print

def run():
  exercise_Compute_phifom_from_abcd_interface(sys.argv[1:])
  print "OK"

if (__name__ == "__main__"):
  run()
