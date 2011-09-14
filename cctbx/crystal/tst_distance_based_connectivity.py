def run(args):
  assert len(args) == 0
  from cctbx.crystal.distance_based_connectivity import \
    build_simple_two_way_bond_sets
  from scitbx.array_family import flex
  sites_cart = flex.vec3_double([
    (25.655, 43.266, 42.630),
    (24.038, 43.853, 43.337),
    (23.048, 44.525, 43.290),
    (21.223, 44.207, 41.475),
    (24.951, 47.170, 37.585),
    (19.298, 46.942, 51.808)])
  elements = flex.std_string(["S", "C", "N", "CU", "ZN", "CA"])
  bond_list = build_simple_two_way_bond_sets(
    sites_cart=sites_cart, elements=elements)
  assert [sorted(b) for b in bond_list] == [[1], [0,2], [1,3], [2], [], []]
  #
  # caffeine
  sites_cart = flex.vec3_double([
    (-2.986, 0.015, 1.643),
    (-1.545, 0.015, 1.643),
    (-0.733, 0.015, 2.801),
    (0.592, 0.015, 2.395),
    (0.618, 0.015, 1.034),
    (1.758, 0.015, 0.102),
    (3.092, -0.06, 0.694),
    (1.525, 0.015, -1.360),
    (2.489, -0.024, -2.139),
    (0.158, 0.015, -1.888),
    (-0.025, 0.024, -3.330),
    (-0.986, 0.015, -0.959),
    (-2.155, 0.008, -1.408),
    (-0.733, 0.015, 0.565),
    (-3.346, 0.016, 2.662),
    (-3.347, 0.896, 1.133),
    (-3.347, -0.868, 1.136),
    (-1.083, 0.02, 3.822),
    (3.184, -0.975, 1.26),
    (3.245, 0.785, 1.348),
    (3.835, -0.047, -0.09),
    (0.508, 0.861, -3.756),
    (-1.076, 0.113, -3.560),
    (0.358, -0.896, -3.748)
  ])
  elements = flex.std_string([
    ' C', ' N', ' C', ' N', ' C', ' N', ' C', ' C', ' O', ' N', ' C', ' C',
    ' O', ' C', ' H', ' H', ' H', ' H', ' H', ' H', ' H', ' H', ' H', ' H'])
  bonds = build_simple_two_way_bond_sets(
    sites_cart=sites_cart, elements=elements)
  assert bonds.size() == sites_cart.size()
  assert list(bonds[0]) == [1, 14, 15, 16]
  #
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
