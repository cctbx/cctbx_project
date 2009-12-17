def run(args):
  assert len(args) == 0
  from cctbx.eltbx.distance_based_connectivity import build_edge_list
  from scitbx.array_family import flex
  sites_cart = flex.vec3_double([
    (25.655, 43.266, 42.630),
    (24.038, 43.853, 43.337),
    (23.048, 44.525, 43.290),
    (21.223, 44.207, 41.475),
    (24.951, 47.170, 37.585),
    (19.298, 46.942, 51.808)])
  elements = ["S", "C", "N", "CU", "ZN", "CA"]
  edge_list = build_edge_list(sites_cart=sites_cart, elements=elements)
  assert edge_list == [(0,1), (1,2)]
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
