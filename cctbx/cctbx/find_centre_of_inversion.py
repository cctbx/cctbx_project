from cctbx import sgtbx
from scitbx.python_utils import list_algebra
from scitbx.python_utils.misc import adopt_init_args

def pair_sort_function(pair_a, pair_b):
  return cmp(pair_a[0], pair_b[0])

class match_refine:

  def __init__(self, equiv_sites, i_pivot1, i_pivot2, tolerance):
    assert len(equiv_sites) > 0
    adopt_init_args(self, locals())
    self.singles = range(len(equiv_sites))
    self.singles.remove(self.i_pivot1)
    if (self.i_pivot1 != self.i_pivot2):
      self.singles.remove(self.i_pivot2)
    self.pairs = [(self.i_pivot1, self.i_pivot2)]
    self.t_ctr = list_algebra.plus(
      self.equiv_sites[i_pivot1].coordinates()[0],
      self.equiv_sites[i_pivot2].coordinates()[0])
    self.unit_cell = equiv_sites[0].unit_cell()
    self.add_pairs()
    self.eliminate_weak_pairs()
    self.pairs.sort(pair_sort_function)

  def add_pairs(self):
    while (len(self.singles)):
      shortest_dist = 2 * self.tolerance
      new_pair = 0
      for is1 in self.singles:
        s1 = self.apply_t_ctr(is1)
        for is2 in self.singles:
          dist = sgtbx.min_sym_equiv_distance_info(
            self.equiv_sites[is2], s1).dist()
          if (dist < shortest_dist):
            shortest_dist = dist
            new_pair = (is1, is2)
      if (new_pair == 0):
        break
      self.pairs.append(new_pair)
      self.singles.remove(new_pair[0])
      if (new_pair[0] != new_pair[1]):
        self.singles.remove(new_pair[1])
      self.refine_t_ctr()

  def eliminate_weak_pairs(self):
    while 1:
      weak_pair = 0
      max_dist = 0
      for pair in self.pairs[1:]:
        dist = self.calculate_pair_info("distance", pair)
        if (dist > max_dist):
          weak_pair = pair
          max_dist = dist
      if (weak_pair == 0): break
      if (max_dist < self.tolerance):
        dist = self.calculate_pair_info("distance", self.pairs[0])
        if (dist < self.tolerance):
          break
      assert len(self.pairs) > 1
      self.pairs.remove(weak_pair)
      self.singles.append(weak_pair[0])
      if (weak_pair[0] != weak_pair[1]):
        self.singles.append(weak_pair[1])
      self.refine_t_ctr()

  def apply_t_ctr(self, i_site):
    return [-x+t for x,t in zip(self.equiv_sites[i_site].coordinates()[0],
                                self.t_ctr)]

  def calculate_pair_info(self, request, pair):
    assert request in ("sites", "distance")
    s1 = self.equiv_sites[pair[0]].coordinates()[0]
    s2 = self.apply_t_ctr(pair[0])
    dist_info = sgtbx.min_sym_equiv_distance_info(
      self.equiv_sites[pair[1]], s2)
    if (request == "sites"):
      return s1, dist_info.sym_op() * s2
    return dist_info.dist()

  def refine_t_ctr(self):
    sum_s1_plus_s2_cart = [0,0,0]
    for pair in self.pairs:
      s1, s2 = self.calculate_pair_info("sites", pair)
      s1_cart = self.unit_cell.orthogonalize(s1)
      s2_cart = self.unit_cell.orthogonalize(s2)
      s1_plus_s2 = list_algebra.plus(s1_cart, s2_cart)
      sum_s1_plus_s2_cart = list_algebra.plus(sum_s1_plus_s2_cart, s1_plus_s2)
    mean_s1_plus_s2_cart = [s / len(self.pairs) for s in sum_s1_plus_s2_cart]
    mean_s1_plus_s2_frac = self.unit_cell.fractionalize(mean_s1_plus_s2_cart)
    self.t_ctr = mean_s1_plus_s2_frac

def match_sort_function(match_a, match_b):
  return -cmp(len(match_a.pairs), len(match_b.pairs))

def weed_refined_matches(refined_matches):
  n_matches = len(refined_matches)
  if (n_matches == 0): return
  best_n_pairs = len(refined_matches[0].pairs)
  is_redundant = [0] * n_matches
  for i in xrange(n_matches-1):
    match_i = refined_matches[i]
    if (is_redundant[i]): continue
    for j in xrange(i+1, n_matches):
      match_j = refined_matches[j]
      if (match_i.pairs == match_j.pairs):
        is_redundant[j] = 1
  for i in xrange(n_matches-1, -1, -1):
    if (is_redundant[i]):
      del refined_matches[i]

def find_matches(xray_structure):
  #    (-1, t_ctr) * site_i_0 = site_j_0
  # => -site_i_0 + t_ctr = site_j_0
  # => t_ctr = site_i_0 + site_j_0
  #    site_j_s = (-site_i_s + t_ctr) + delta_i_j
  # => delta_i_j = site_j_s + site_i_s - t_ctr
  # minimize: sum(delta_i_j**2) as a function of t_ctr
  # => delta_i_j**2 = (site_j_s + site_i_s - t_ctr)**2
  # => d(delta_i_j**2)/d(t_ctr) = -2*(site_j_s + site_i_s - t_ctr)
  # => minimum if sum(-2*(site_j_s + site_i_s - t_ctr)) == 0
  # => t_ctr = sum(site_i_s + site_j_s) / n
  equiv_sites = []
  for scatterer in xray_structure.scatterers():
    site_symmetry = xray_structure.site_symmetry(scatterer.site)
    equiv_sites.append(sgtbx.sym_equiv_sites(site_symmetry))
  refined_matches = []
  n_scatterers = xray_structure.scatterers().size()
  for i_scatterer in xrange(n_scatterers):
    for j_scatterer in xrange(i_scatterer, n_scatterers):
      match = match_refine(equiv_sites, i_scatterer, j_scatterer, tolerance=1)
      refined_matches.append(match)
  refined_matches.sort(match_sort_function)
  weed_refined_matches(refined_matches)
  return refined_matches
