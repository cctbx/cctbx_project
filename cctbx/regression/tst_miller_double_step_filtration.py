from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from cctbx import miller
import iotbx.pdb

def map_cc(map_coeffs_1, map_coeffs_2):
  fft_map_1 = map_coeffs_1.fft_map(resolution_factor=0.25)
  map_1 = fft_map_1.real_map_unpadded()
  fft_map_2 = miller.fft_map(
    crystal_gridding = fft_map_1,
    fourier_coefficients = map_coeffs_2)
  map_2 = fft_map_2.real_map_unpadded()
  assert map_1.size() == map_2.size()
  m1 = flex.double(list(map_1))
  m2 = flex.double(list(map_2))
  return flex.linear_correlation(x = m1, y = m2).coefficient()

def run():
  import iotbx.pdb
  xrs_str = """
CRYST1    8.000    8.000    8.000  90.00  90.00  90.00 P 1
HETATM  115  O   HOH A  18       4.000   4.000   4.000  1.00 10.00           O
TER
END
  """
  pdb_inp = iotbx.pdb.input(source_info=None, lines=xrs_str)
  xrs = pdb_inp.xray_structure_simple()
  xrs.scattering_type_registry(table = "wk1995")
  cntr = 0
  for remove_fraction in [0.2, 0.5]:
    print
    print "remove_fraction:", remove_fraction
    for option in ["random", "smallest", "highest"]:
      print
      print "data incompleteness:",option,"-"*30
      #
      f_exact = xrs.structure_factors(d_min = 1.0).f_calc()
      #
      if(option=="highest"):
        s = flex.sort_permutation(abs(f_exact).data(), reverse=True)
        f_exact = f_exact.select(s)
        n_remove = int(s.size()*remove_fraction)
        f_poor = f_exact.customized_copy(
          data    = f_exact.data()[n_remove:],
          indices = f_exact.indices()[n_remove:])
      elif(option == "smallest"):
        s = flex.sort_permutation(abs(f_exact).data(), reverse=True)
        f_exact = f_exact.select(s)
        n_remove = int(s.size()*remove_fraction)
        sz = f_exact.data().size()
        f_poor = f_exact.customized_copy(
          data    = f_exact.data()[:sz-n_remove],
          indices = f_exact.indices()[:sz-n_remove])
      elif(option == "random"):
        s = flex.random_bool(f_exact.data().size(), 1.-remove_fraction)
        f_poor = f_exact.select(s)
      else: assert 0
      #
      print "number of all data:", f_exact.data().size()
      print "number of incomplete data:", f_poor.data().size()
      cc1 = map_cc(map_coeffs_1=f_exact, map_coeffs_2=f_poor)
      print "start CC(exact_map, poor_map): ", cc1
      #
      f_dsf = f_poor.double_step_filtration(
        vol_cutoff_plus_percent =0.1,
        vol_cutoff_minus_percent=0.1,
        complete_set=f_exact)
      f_new = f_poor.complete_with(other = f_dsf)
      cc2 = map_cc(map_coeffs_1=f_exact, map_coeffs_2=f_new)
      print "start CC(exact_map, filled_map)1: ", cc2
      #
      f_dsf = f_poor.double_step_filtration(
        vol_cutoff_plus_percent =0.1,
        vol_cutoff_minus_percent=0.1,
        complete_set=f_exact,
        scale_to=f_exact)
      f_new = f_poor.complete_with(other = f_dsf)
      cc3 = map_cc(map_coeffs_1=f_exact, map_coeffs_2=f_new)
      print "start CC(exact_map, filled_map)2: ", cc3
      #
      if(option=="highest"):
        if(remove_fraction==0.2):
          assert cc1<0.92
          assert cc2>0.99 and cc3>0.99
          cntr += 1
        elif(remove_fraction==0.5):
          assert cc1<0.16
          assert cc2>0.76
          assert cc3>0.99
          cntr += 1
  #
  assert cntr == 2 # make sure it's gone through all if statements

if (__name__ == "__main__"):
  run()
