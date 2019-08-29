from __future__ import absolute_import, division, print_function
from scitbx.graph.tardy_tree import find_paths, construct
from scitbx.graph import rigidity
from scitbx.graph import utils
import sys
from six.moves import range

def exercise_minimal():
  edge_sets = utils.construct_edge_sets(n_vertices=1, edge_list=[])
  assert find_paths(edge_sets=edge_sets).search_from(iv=0) == ({}, {})
  edge_sets = utils.construct_edge_sets(n_vertices=2, edge_list=[(0,1)])
  for iv in [0,1]:
    assert find_paths(edge_sets=edge_sets).search_from(iv=iv) == ({}, {})

def exercise_simple_loops(loop_size_max=10):
  for n_vertices in range(3, loop_size_max+1):
    edge_list = [tuple(sorted((i,(i+1)%n_vertices)))
      for i in range(n_vertices)]
    edge_sets = utils.construct_edge_sets(
      n_vertices=n_vertices, edge_list=edge_list)
    loops, dendrites = find_paths(edge_sets=edge_sets).search_from(iv=0)
    if (n_vertices <= 7):
      assert len(loops) == 2
      assert sorted(dendrites.keys()) == list(range(1,n_vertices))
    else:
      assert len(loops) == 0
      assert sorted(dendrites.keys()) == list(range(2,n_vertices-1))
    if (n_vertices == 3):
      assert loops == {1: [[2]], 2: [[1]]}
      assert dendrites == {1: [set([2])], 2: [set([1])]}
    elif (n_vertices == 4):
      assert loops == {1: [[2, 3]], 3: [[2, 1]]}
      assert dendrites == {
        1: [set([2, 3])], 2: [set([1]), set([3])], 3: [set([1, 2])]}

def exercise_knot():
  edge_sets = utils.construct_edge_sets(
    n_vertices=4,
    edge_list=[(0,1), (1,2), (2,3), (1,3)])
  expected_loops_dendrites = [
    ({}, {2: [set([1]), set([1, 3])], 3: [set([1, 2]), set([1])]}),
    ({2: [[3]], 3: [[2]]}, {2: [set([3])], 3: [set([2])]}),
    ({1: [[3]], 3: [[1]]}, {0: [set([1])], 1: [set([3])], 3: [set([1])]}),
    ({1: [[2]], 2: [[1]]}, {0: [set([1])], 1: [set([2])], 2: [set([1])]})]
  for iv in range(4):
    loop_dendrites = find_paths(edge_sets=edge_sets).search_from(iv=iv)
    assert loop_dendrites == expected_loops_dendrites[iv]

def archs_grow_edge_list(edge_list, offs, size, av=0, bv=1):
  result = list(edge_list)
  i = av
  for j in range(offs, offs+size):
    result.append((i,j))
    i = j
  result.append((bv,i))
  return result

def arch_dof(n_vertices, edge_list):
  es = utils.construct_edge_sets(n_vertices=n_vertices, edge_list=edge_list)
  bbes = utils.bond_bending_edge_sets(edge_sets=es)
  bbel = utils.extract_edge_list(edge_sets=bbes)
  dofs = [rigidity.determine_degrees_of_freedom(
    n_dim=3, n_vertices=n_vertices, edge_list=bbel, method=method)
      for method in ["float", "integer"]]
  assert dofs[0] == dofs[1]
  return es, dofs[0]

def exercise_fused_loops(arch_size_max=8):
  for arch_size_1 in range(1, arch_size_max+1):
    edge_list_1 = archs_grow_edge_list(
      [(0,1)], 2, arch_size_1)
    for arch_size_2 in range(1, arch_size_max+1):
      n_vertices = 2 + arch_size_1 + arch_size_2
      edge_list_12 = archs_grow_edge_list(
        edge_list_1, 2+arch_size_1, arch_size_2)
      es, dof = arch_dof(n_vertices=n_vertices, edge_list=edge_list_12)
      is_rigid = (dof == 6)
      inferred_is_rigid = (
            arch_size_1 < 6
        and arch_size_2 < 6
        and arch_size_1 + arch_size_2 < 10)
      assert inferred_is_rigid == is_rigid
      #
      tt = construct(
        n_vertices=n_vertices, edge_list=edge_list_12).build_tree()
      inferred_is_rigid = (len(tt.cluster_manager.clusters) == 1)
      assert inferred_is_rigid == is_rigid

def exercise_three_archs(arch_size_max=8):
  for arch_size_1 in range(1, arch_size_max+1):
    edge_list_1 = archs_grow_edge_list(
      [], 2, arch_size_1)
    for arch_size_2 in range(1, arch_size_max+1):
      edge_list_12 = archs_grow_edge_list(
        edge_list_1, 2+arch_size_1, arch_size_2)
      for arch_size_3 in range(1, arch_size_max+1):
        n_vertices = 2 + arch_size_1 + arch_size_2 + arch_size_3
        edge_list_123 = archs_grow_edge_list(
          edge_list_12, 2+arch_size_1+arch_size_2, arch_size_3)
        es, dof = arch_dof(n_vertices=n_vertices, edge_list=edge_list_123)
        expected = max(
          6,
          max(arch_size_1, arch_size_2, arch_size_3) + 1,
          arch_size_1 + arch_size_2 + arch_size_3 - 3)
        assert expected == dof
        is_rigid = (dof == 6)
        inferred_is_rigid = (
              arch_size_1 < 6
          and arch_size_2 < 6
          and arch_size_3 < 6
          and arch_size_1 + arch_size_2 + arch_size_3 < 10)
        assert inferred_is_rigid == is_rigid
        #
        tt = construct(
          n_vertices=n_vertices, edge_list=edge_list_123).build_tree()
        inferred_is_rigid = (len(tt.cluster_manager.clusters) == 1)

def run(args):
  assert len(args) == 0
  exercise_minimal()
  exercise_simple_loops()
  exercise_knot()
  exercise_fused_loops()
  exercise_three_archs()
  print("OK")

if (__name__ == "__main__"):
  run(sys.argv[1:])
