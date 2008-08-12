from mmtbx.monomer_library import rna_sugar_pucker_analysis
from libtbx.test_utils import approx_equal
import sys

def exercise(args):
  assert len(args) == 0
  params = rna_sugar_pucker_analysis.master_phil.fetch().extract()
  analysis = rna_sugar_pucker_analysis.evaluate(
    params=params,
    sites_cart_1=[
      (21.661,25.659,31.824), (21.101,30.048,33.325), (22.546,31.003,35.039),
      (22.184,28.704,36.539), (21.527,28.782,35.281), (22.379,28.287,34.119),
      (22.676,26.795,34.032), (19.628,30.157,32.953)],
    sites_cart_2=[
      (33.058,33.903,32.613), (36.334,31.063,34.876), (37.147,29.459,33.274),
      (35.933,30.835,31.292), (35.553,31.218,32.607), (36.283,32.484,33.032),
      (35.685,33.793,32.540), (35.480,30.849,36.050)])
  assert analysis.epsilon is None
  assert analysis.delta is None
  assert analysis.p_distance_c1_n_line is None
  assert analysis.is_2p_epsilon is None
  assert analysis.is_2p_delta is None
  assert analysis.is_2p_p_distance_c1_n_line is None
  assert analysis.is_2p is None
  analysis = rna_sugar_pucker_analysis.evaluate(
    params=params,
    sites_cart_1=[
      (53.886,12.435,32.393), (50.473,9.942,36.204), (48.109,9.956,36.120),
      (48.952,10.236,33.644), (49.474,11.213,34.497), (50.978,11.174,34.278),
      (51.710,12.501,34.095), (50.926,10.084,37.569)],
    sites_cart_2=[
      (48.114,10.641,32.347), (43.341,8.882,31.694), (44.062,6.631,32.479),
      (43.515,8.286,34.518), (44.572,8.714,33.657), (44.396,10.191,33.275),
      (45.664,11.004,33.132), (43.317,8.783,30.235)])
  assert approx_equal(analysis.epsilon, 250.715932662)
  assert approx_equal(analysis.delta, 133.811229901)
  assert approx_equal(analysis.p_distance_c1_n_line, 1.52374220064)
  assert not analysis.is_2p_epsilon
  assert analysis.is_2p_delta
  assert analysis.is_2p_p_distance_c1_n_line
  assert analysis.is_2p
  analysis = rna_sugar_pucker_analysis.evaluate(
    params=params,
    sites_cart_1=[
      (44.576,14.676,36.147), (42.736,12.998,31.661), (44.142,14.508,30.426),
      (42.922,16.529,31.706), (42.428,15.271,31.984), (43.132,14.710,33.185),
      (42.545,15.146,34.511), (41.433,12.332,31.466)],
    sites_cart_2=[
      (41.856,17.581,31.191), (40.301,14.847,27.279), (40.175,16.098,25.250),
      (40.845,18.224,26.422), (40.676,17.174,27.296), (41.974,16.449,27.569),
      (42.684,17.049,28.754), (39.245,14.368,28.207)])
  assert not analysis.is_2p_epsilon
  assert not analysis.is_2p_delta
  assert not analysis.is_2p_p_distance_c1_n_line
  assert not analysis.is_2p
  #
  params.epsilon_range_not_2p_min = None
  params.epsilon_range_not_2p_max = None
  params.delta_range_2p_min = None
  params.delta_range_2p_max = None
  params.p_distance_c1_n_line_2p_max = None
  analysis = rna_sugar_pucker_analysis.evaluate(
    params=params,
    sites_cart_1=[
      (53.886,12.435,32.393), (50.473,9.942,36.204), (48.109,9.956,36.120),
      (48.952,10.236,33.644), (49.474,11.213,34.497), (50.978,11.174,34.278),
      (51.710,12.501,34.095), (50.926,10.084,37.569)],
    sites_cart_2=[
      (48.114,10.641,32.347), (43.341,8.882,31.694), (44.062,6.631,32.479),
      (43.515,8.286,34.518), (44.572,8.714,33.657), (44.396,10.191,33.275),
      (45.664,11.004,33.132), (43.317,8.783,30.235)])
  assert approx_equal(analysis.epsilon, 250.715932662)
  assert approx_equal(analysis.delta, 133.811229901)
  assert approx_equal(analysis.p_distance_c1_n_line, 1.52374220064)
  assert analysis.is_2p_epsilon is None
  assert analysis.is_2p_delta is None
  assert analysis.is_2p_p_distance_c1_n_line is None
  assert analysis.is_2p is None
  #
  print "OK"

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
