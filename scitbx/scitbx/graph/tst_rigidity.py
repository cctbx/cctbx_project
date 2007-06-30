from scitbx.graph.rigidity import gcd, determine_degrees_of_freedom
import sys

def exercise_gcd():
  assert gcd(8,0) == 8
  assert gcd(-4,0) == 4
  assert gcd(1,1) == 1
  assert gcd(10,15) == 5
  assert gcd(-18,42) == 6

def exercise_double_banana():
  n_vertices = 8
  edge_list = [
    (0,1), (0,2), (1,2),
    (0,3), (1,3), (2,3),
    (0,4), (1,4), (2,4),
    (5,6), (5,7), (6,7),
    (5,3), (6,3), (7,3),
    (5,4), (6,4), (7,4)]
  assert len(edge_list) == 18
  assert determine_degrees_of_freedom(
    n_dim=2, n_vertices=n_vertices, edge_list=edge_list) == 3
  assert determine_degrees_of_freedom(
    n_dim=3, n_vertices=n_vertices, edge_list=edge_list) == 7
  assert determine_degrees_of_freedom(
    n_dim=4, n_vertices=n_vertices, edge_list=edge_list) == 14
  # remove each edge in turn
  for remove in xrange(len(edge_list)):
    assert determine_degrees_of_freedom(
      n_dim=3,
      n_vertices=n_vertices,
      edge_list=edge_list[:remove]+edge_list[remove+1:]) == 7
  # add all possible edges
  dofs = []
  for i in xrange(n_vertices-1):
    for j in xrange(i+1,n_vertices):
      dofs.append(determine_degrees_of_freedom(
        n_dim=3, n_vertices=n_vertices, edge_list=edge_list+[(i,j)]))
  assert dofs == [
    7, 7, 7, 7, 6, 6, 6, 7, 7, 7, 6, 6, 6, 7,
    7, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]

def exercise_mt1996():
  # H. Maehara and N. Tokushige
  # A Spatial Unit-Bar-Framework Which Is Rigid and Triangle-Free
  # Graphs and Combinatorics (1996) 12:341-344
  n_vertices = 26
  edge_list = [
    (0,8), (0,9), (0,12), (0,13), (0,16), (0,17), (0,20), (0,21),
    (0,22), (0,23), (1,8), (1,9), (1,12), (1,13), (1,18), (1,19),
    (1,20), (1,21), (1,24), (1,25), (2,8 ), (2,9), (2,14), (2,15),
    (2,18), (2,19), (2,22), (2,23), (2,24), (2,25), (3,8), (3,9),
    (3,14), (3,15), (3,16), (3,17), (4,10), (4,11), (4,12), (4,13),
    (4,16), (4,17), (4,22), (4,23), (4,24), (4,25), (5,10), (5,11),
    (5,12), (5,13), (5,18), (5,19), (6,10), (6,11), (6,14), (6,15),
    (6,18), (6,19), (6,20), (6,21), (6,22), (6,23), (7,10), (7,11),
    (7,14), (7,15), (7,16), (7,17), (7,20), (7,21), (7,24), (7,25),
    (8,11), (9,10), (12,15), (13,14), (16,19), (17,18)]
  assert len(edge_list) == 78
  assert determine_degrees_of_freedom(
    n_dim=3, n_vertices=n_vertices, edge_list=edge_list) == 6

def exercise_k6_6_minus_six_parallel_edges():
  # "K6,6 minus six parallel edges" (Figure 3.23 of J.E. Graver,
  # Counting on Frameworks, 2001).
  n_vertices = 12
  edge_list = [
    (0,7), (0,8), (0,9), (0,10), (0,11), (1,6), (1,8), (1,9), (1,10),
    (1,11), (2,6), (2,7), (2,9), (2,10), (2,11), (3,6), (3,7), (3,8),
    (3,10), (3,11), (4,6), (4,7), (4,8), (4,9), (4,11), (5,6), (5,7),
    (5,8), (5,9), (5,10)]
  assert len(edge_list) == 30
  assert determine_degrees_of_freedom(
    n_dim=3, n_vertices=n_vertices, edge_list=edge_list) == 6

def exercise_p120():
  # http://www.rwgrayprojects.com/Lynn/Coordinates/coord01.html
  n_vertices = 62
  edge_list = [
    (0,1),(0,3),(0,5),(0,7),(1,2),(1,3),(1,7),(1,8),
    (1,9),(1,10),(1,17),(1,18),(1,19),(2,3),(2,10),
    (2,11),(3,4),(3,5),(3,11),(4,5),(4,11),(4,12),(5,6),
    (5,7),(5,12),(5,13),(5,14),(5,15),(5,22),(6,7),
    (6,15),(6,16),(7,8),(7,16),(8,16),(8,17),(9,10),(9,19),
    (9,26),(10,11),(10,20),(10,26),(11,12),(11,20),
    (11,27),(11,28),(11,21),(11,29),(12,13),(12,21),(12,30),
    (13,22),(13,30),(14,15),(14,22),(14,32),(15,16),
    (15,23),(15,32),(16,17),(16,23),(16,24),(16,33),(16,34),
    (16,35),(17,18),(17,24),(17,36),(18,19),(18,36),
    (19,25),(19,26),(19,36),(20,26),(20,27),(21,29),(21,30),
    (22,30),(22,31),(22,32),(23,32),(23,33),(24,35),
    (24,36),(25,26),(25,36),(25,37),(26,27),(26,37),(26,38),
    (26,43),(26,44),(27,28),(27,38),(27,45),(28,29),
    (28,45),(29,30),(29,39),(29,45),(30,31),(30,39),(30,40),
    (30,46),(30,47),(31,32),(31,40),(32,33),(32,40),
    (32,41),(32,48),(32,49),(33,34),(33,41),(33,50),(34,35),
    (34,50),(35,36),(35,42),(35,50),(36,37),(36,42),
    (36,51),(36,52),(37,43),(37,52),(37,53),(38,44),(38,45),
    (39,45),(39,46),(40,47),(40,48),(40,57),(41,49),
    (41,50),(42,50),(42,51),(43,44),(43,53),(44,54),(44,45),
    (45,46),(45,54),(45,55),(45,56),(46,47),(46,56),
    (46,57),(47,57),(48,49),(48,57),(49,50),(49,57),(49,58),
    (50,51),(50,58),(50,59),(50,60),(51,52),(51,60),
    (51,53),(52,53),(53,54),(53,55),(53,59),(53,60),(53,61),
    (54,55),(55,56),(55,57),(55,61),(56,57),(57,58),
    (57,59),(57,61),(58,59),(59,60),(59,61)]
  assert len(edge_list) == 179
  assert determine_degrees_of_freedom(
    n_dim=3, n_vertices=n_vertices, edge_list=edge_list) == 7
  assert determine_degrees_of_freedom(
    n_dim=3, n_vertices=n_vertices, edge_list=edge_list+[(28,24)]) == 6

def exercise():
  exercise_gcd()
  exercise_double_banana()
  exercise_mt1996()
  exercise_k6_6_minus_six_parallel_edges()
  if ("--all" in sys.argv[1:]):
    exercise_p120()
  print "OK"

if (__name__ == "__main__"):
  exercise()
