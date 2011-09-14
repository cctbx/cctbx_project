def exercise_cubicles_max_memory():
  import scitbx.cubicle_neighbors as cn
  assert cn.cubicles_max_memory_allocation_get() != 0
  mm = cn.cubicles_max_memory_allocation_get()
  cn.cubicles_max_memory_allocation_set(number_of_bytes=10)
  assert cn.cubicles_max_memory_allocation_get() == 10
  cn.cubicles_max_memory_allocation_set(number_of_bytes=0)
  assert cn.cubicles_max_memory_allocation_get() == 0
  cn.cubicles_max_memory_allocation_set(number_of_bytes=mm)
  assert cn.cubicles_max_memory_allocation_get() == mm
  # more tests in cctbx/crystal/tst_ext.py, exercise_cubicles_max_memory()

def neighbors_simple(main_sites_cart, other_sites_cart, distance_cutoff_sq):
  from scitbx.matrix import col
  result = {}
  for j,sj in enumerate(other_sites_cart):
    vj = col(sj)
    for i,si in enumerate(main_sites_cart):
      vi = col(si)
      if ((vj-vi).length_sq() <= distance_cutoff_sq):
        result.setdefault(j, []).append(i)
  return result

def run(args):
  assert len(args) == 0
  exercise_cubicles_max_memory()
  from scitbx.cubicle_neighbors import cubicle_neighbors
  from scitbx.array_family import flex
  main_sites_cart = flex.vec3_double()
  cn = cubicle_neighbors(main_sites_cart=main_sites_cart, cubicle_edge=5)
  nb = cn.neighbors_of(other_sites_cart=main_sites_cart, distance_cutoff_sq=1)
  assert nb.size() == 0
  for xyz in [(0,0,0), (0.1, 0.2, -0.3)]:
    main_sites_cart = flex.vec3_double([xyz])
    cn = cubicle_neighbors(main_sites_cart=main_sites_cart, cubicle_edge=5)
    nb = cn.neighbors_of(other_sites_cart=main_sites_cart, distance_cutoff_sq=1)
    assert nb.size() == 1
    assert nb.keys() == [0]
    assert list(nb[0]) == [0]
    nb = cn.neighbors_of(
      other_sites_cart=flex.vec3_double([(2,2,2)]), distance_cutoff_sq=1)
    assert nb.size() == 0
    nb = cn.neighbors_of(
      other_sites_cart=flex.vec3_double([(2,2,2)]), distance_cutoff_sq=25)
    assert nb.size() == 1
  mt = flex.mersenne_twister(seed=0)
  for nm in [3,5,8]:
    for no in [1,7,9]:
      main_sites_cart = flex.vec3_double(zip(
        mt.random_double(size=nm)*2-1,
        mt.random_double(size=nm)*2-1,
        mt.random_double(size=nm)*2-1))
      other_sites_cart = flex.vec3_double(zip(
        mt.random_double(size=no)*2-1,
        mt.random_double(size=no)*2-1,
        mt.random_double(size=no)*2-1))
      for distance_cutoff in [0.5, 1]:
        distance_cutoff_sq = distance_cutoff**2
        cn = cubicle_neighbors(main_sites_cart=main_sites_cart, cubicle_edge=1)
        nb = cn.neighbors_of(
          other_sites_cart=other_sites_cart,
          distance_cutoff_sq=distance_cutoff_sq)
        nb_simple = neighbors_simple(
          main_sites_cart=main_sites_cart,
          other_sites_cart=other_sites_cart,
          distance_cutoff_sq=distance_cutoff_sq)
        assert sorted(nb.keys()) == sorted(nb_simple.keys())
        for j_seq,i_seqs_simple in nb_simple.items():
          i_seqs = nb[j_seq]
          assert sorted(i_seqs) == sorted(i_seqs_simple)
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
