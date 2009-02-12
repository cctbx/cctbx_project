from scitbx.graph.tardy_tree import find_paths
from scitbx.graph import rigidity
from scitbx.graph import utils
from libtbx.utils import null_out
import sys

def exercise_simple_loops(out, loop_size_max=8):
  for n_vertices in xrange(3, loop_size_max+1):
    edge_list = [tuple(sorted((i,(i+1)%n_vertices)))
      for i in xrange(n_vertices)]
    edge_sets = utils.construct_edge_sets(
      n_vertices=n_vertices, edge_list=edge_list)
    jps = find_paths(edge_sets=edge_sets, out=out)
    if (n_vertices <= 6):
      assert len(jps[0]) == 2
      for jp in jps[0]:
        assert len(jp[1][0]) == n_vertices-2
    else:
      assert len(jps[0]) == 0

def three_archs_grow_edge_list(edge_list, offs, size):
  result = list(edge_list)
  i = 0
  for j in xrange(offs, offs+size):
    result.append((i,j))
    i = j
  result.append((1,i))
  return result

def exericse_three_archs(out, arch_size_max=8):
  for arch_size_1 in xrange(1, arch_size_max):
    edge_list_1 = three_archs_grow_edge_list(
      [], 2, arch_size_1)
    for arch_size_2 in xrange(1, arch_size_max):
      edge_list_12 = three_archs_grow_edge_list(
        edge_list_1, 2+arch_size_1, arch_size_2)
      for arch_size_3 in xrange(1, arch_size_max):
        n_vertices = 2 + arch_size_1 + arch_size_2 + arch_size_3
        edge_list_123 = three_archs_grow_edge_list(
          edge_list_12, 2+arch_size_1+arch_size_2, arch_size_3)
        es = utils.construct_edge_sets(
          n_vertices=n_vertices, edge_list=edge_list_123)
        bbes = utils.bond_bending_edge_sets(edge_sets=es)
        bbel = utils.extract_edge_list(edge_sets=bbes)
        dofs = [rigidity.determine_degrees_of_freedom(
          n_dim=3, n_vertices=n_vertices, edge_list=bbel, method=method)
            for method in ["float", "integer"]]
        assert dofs[0] == dofs[1]
        print >> out, " ", arch_size_1, arch_size_2, arch_size_3, dofs[0]-6
        print >> out, "edge_list:", edge_list_123
        expected = max(
          6,
          max(arch_size_1, arch_size_2, arch_size_3) + 1,
          arch_size_1 + arch_size_2 + arch_size_3 - 3)
        assert expected == dofs[0]
        is_rigid = dofs[0] == 6
        inferred_is_rigid = (
              arch_size_1 < 6
          and arch_size_2 < 6
          and arch_size_3 < 6
          and arch_size_1 + arch_size_2 + arch_size_3 < 10)
        assert inferred_is_rigid == is_rigid
        s = (arch_size_1 < 6) \
          + (arch_size_2 < 6) \
          + (arch_size_3 < 6)
        jps = find_paths(edge_sets=es, out=out)
        assert len(jps[1]) == s
        if (len(jps[1]) == 3):
          inferred_is_rigid = sum([len(jp[1][0]) for jp in jps[1]]) < 7
          assert inferred_is_rigid == is_rigid

def run(args):
  assert args in [[], ["--verbose"]]
  if (len(args) == 0):
    out = null_out()
  else:
    out = sys.stdout
  exercise_simple_loops(out=out)
  exericse_three_archs(out=out)
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
