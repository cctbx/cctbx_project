from scitbx.graph.tardy_tree import find_paths
from scitbx.graph import rigidity
from scitbx.graph import utils
import sys

def exercise_minimal():
  edge_sets = utils.construct_edge_sets(n_vertices=1, edge_list=[])
  assert find_paths(edge_sets=edge_sets, iv=0) == {0: {}}
  edge_sets = utils.construct_edge_sets(n_vertices=2, edge_list=[(0,1)])
  for iv in [0,1]:
    assert find_paths(edge_sets=edge_sets, iv=iv) == {0: {1: []}, 1: {0: []}}

def exercise_simple_loops(loop_size_max=8):
  for n_vertices in xrange(3, loop_size_max+1):
    edge_list = [tuple(sorted((i,(i+1)%n_vertices)))
      for i in xrange(n_vertices)]
    edge_sets = utils.construct_edge_sets(
      n_vertices=n_vertices, edge_list=edge_list)
    jv_kv_paths = find_paths(edge_sets=edge_sets, iv=0)
    assert len(jv_kv_paths[0]) == 2
    if (n_vertices <= 6):
      jv_kv_paths[0][1] == n_vertices-2
      jv_kv_paths[0][n_vertices-1] == n_vertices-2
    else:
      jv_kv_paths[0][1] == 0
      jv_kv_paths[0][n_vertices-1] == 0
    if (n_vertices == 3):
      assert jv_kv_paths == {
        0: {1: [[2]], 2: [[1]]},
        1: {0: [], 2: []},
        2: {0: [], 1: []}}
    elif (n_vertices == 4):
      assert jv_kv_paths == {
        0: {1: [[3, 2]], 3: [[1, 2]]},
        1: {0: [], 2: [[3]]},
        2: {1: [], 3: []},
        3: {0: [], 2: [[1]]}}

def exercise_knot():
  edge_sets = utils.construct_edge_sets(
    n_vertices=4,
    edge_list=[(0,1), (1,2), (2,3), (1,3)])
  expected_jv_kb_paths = [
    {0: {1: []},
     1: {0: [], 2: [[1, 3]], 3: [[1, 2]]},
     2: {1: [], 3: [[1]]},
     3: {1: [], 2: [[1]]}},
    {0: {1: []},
     1: {0: [], 2: [[3]], 3: [[2]]},
     2: {1: [], 3: []},
     3: {1: [], 2: []}},
    {0: {1: [[3]]},
     1: {2: [], 3: []},
     2: {1: [[3]], 3: [[1]]},
     3: {1: [], 2: []}},
    {0: {1: [[2]]},
     1: {2: [], 3: []},
     2: {1: [], 3: []},
     3: {1: [[2]], 2: [[1]]}}]
  for iv in xrange(4):
    jv_kv_paths = find_paths(edge_sets=edge_sets, iv=iv)
    assert jv_kv_paths == expected_jv_kb_paths[iv]

def three_archs_grow_edge_list(edge_list, offs, size):
  result = list(edge_list)
  i = 0
  for j in xrange(offs, offs+size):
    result.append((i,j))
    i = j
  result.append((1,i))
  return result

def exericse_three_archs(arch_size_max=8):
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
        for iv in [0,1]:
          jv_kv_paths = find_paths(edge_sets=es, iv=iv)
          kv_paths = jv_kv_paths.get(1-iv, [])
          assert len(kv_paths) == s
          if (len(kv_paths) == 3):
            sum_len = 0
            for paths in kv_paths.values():
              assert len(paths) in [0,1]
              if (len(paths) == 1):
                sum_len += len(paths[0])
            inferred_is_rigid = sum_len < 7
            assert inferred_is_rigid == is_rigid
          else:
            assert not is_rigid

def run(args):
  assert len(args) == 0
  exercise_minimal()
  exercise_simple_loops()
  exercise_knot()
  exericse_three_archs()
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
