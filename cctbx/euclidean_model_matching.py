from __future__ import absolute_import, division, print_function

import operator
from cctbx import crystal
from cctbx import sgtbx
from cctbx.array_family import flex
from scitbx import matrix
from scitbx.python_utils import dicts
from libtbx.utils import user_plus_sys_time
from libtbx import adopt_init_args
import sys, math

import boost_adaptbx.boost.python as bp
from six.moves import range
from six.moves import zip
ext = bp.import_ext("cctbx_emma_ext")

def sgtbx_rt_mx_as_matrix_rt(s):
  return matrix.rt((s.r().as_double(), s.t().as_double()))

class position(object):

  def __init__(self, label, site):
    adopt_init_args(self, locals())

  def __repr__(self):
    return "%-4s %7.4f %7.4f %7.4f" % ((self.label,) + tuple(self.site))

class model(crystal.special_position_settings):

  def __init__(self, special_position_settings, positions=None):
    crystal.special_position_settings._copy_constructor(
      self, special_position_settings)
    self.reset_cb_op()
    self._positions = []
    if (positions is not None):
      self.add_positions(positions)

  def cb_op(self):
    return self._cb_op

  def reset_cb_op(self):
    self._cb_op = sgtbx.change_of_basis_op()
    return self

  def positions(self):
    return self._positions

  def __len__(self):
    return len(self._positions)

  def size(self):
    return len(self._positions)

  def __getitem__(self, key):
    return self._positions[key]

  def add_position(self, pos):
    self._positions.append(position(
      pos.label,
      self.site_symmetry(pos.site).exact_site()))

  def add_positions(self, positions):
    for pos in positions:
      self.add_position(pos)

  def change_basis(self, cb_op):
    positions = []
    for pos in self._positions:
      positions.append(position(pos.label, cb_op(pos.site)))
    result = model(
      crystal.special_position_settings.change_basis(self, cb_op),
      positions)
    result._cb_op = cb_op.new_denominators(self.cb_op()) * self.cb_op()
    return result

  def transform_to_reference_setting(self):
    cb_op = self.space_group_info().type().cb_op()
    result = self.change_basis(cb_op)
    assert result.space_group_info().is_reference_setting()
    return result

  def change_hand(self):
    ch_op = self.space_group_info().type().change_of_hand_op()
    return self.change_basis(ch_op)

  def expand_to_p1(self):
    new_model = model(
      crystal.special_position_settings(
        crystal.symmetry.cell_equivalent_p1(self)))
    for pos in self._positions:
      site_symmetry = self.site_symmetry(pos.site)
      equiv_sites = sgtbx.sym_equiv_sites(site_symmetry)
      i = 0
      for site in equiv_sites.coordinates():
        i += 1
        new_model.add_position(position(
          label=pos.label+"_%03d"%i,
          site=site))
    return new_model

  def show(self, title, f=None):
    if (f is None): f = sys.stdout
    print(title, file=f)
    crystal.special_position_settings.show_summary(self, f)
    if (not self.cb_op().is_identity_op()):
      print("Change of basis:", file=f)
      print("  c:", self.cb_op().c(), file=f)
      print("  c_inv:", self.cb_op().c_inv(), file=f)
    for pos in self.positions(): print(pos, file=f)
    print(file=f)

  def as_xray_structure(self, scatterer=None):
    from cctbx import xray
    if (scatterer is None):
      scatterer = xray.scatterer(scattering_type="const")
    result = xray.structure(special_position_settings=self)
    for position in self.positions():
      result.add_scatterer(scatterer.customized_copy(
        label=position.label,
        site=position.site))
    return result

  def combine_with_other(self, other_model, tolerance = 1.5,
    models_are_diffraction_index_equivalent = True,
    break_if_match_with_no_singles=True, f=sys.stdout,
    improved_only=True,new_model_number=0):
    # 2013-01-25 tt superpose other on this one and return composite
    match_list=self.best_superpositions_on_other(other_model,
      tolerance=tolerance,models_are_diffraction_index_equivalent=
           models_are_diffraction_index_equivalent,
           break_if_match_with_no_singles=break_if_match_with_no_singles,
           f=f,specifically_test_inverse=True)

    new_model_list=[]
    for x in other_model.component_model_numbers:
      if x in self.component_model_numbers:
        return new_model_list # cannot combine with something already used
    for match in match_list:
      if match is None: continue
      if improved_only and (len(match.singles2) ==0):
         continue

      test_new_model= model(special_position_settings=self)
      new_model= model(special_position_settings=self)
      new_model.model_number=new_model_number
      new_model.component_model_numbers=[new_model_number]+ \
         self.component_model_numbers+other_model.component_model_numbers
      new_model_number+=1
      for pos in self.positions():
        new_model.add_position(position(label=pos.label, site=(pos.site)))
      i=new_model.size()-1
      for s in match.singles2:
        site=match.ref_model2[s].site
        new_model.add_position(position(
           label="ATOM_%03d"%i,
           site=(match.rt*site).elems))
        test_new_model.add_position(position(
           label="ATOM_%03d"%i,
           site=(match.rt*site).elems))
      new_model_list.append(new_model)
    return new_model_list


  def best_superpositions_on_other(self, other_model, tolerance = 1.5,
    models_are_diffraction_index_equivalent = True,
    break_if_match_with_no_singles=True, f=sys.stdout,
    specifically_test_inverse=True):
    # 2013-01-25 tt.Find best match to other_model and return it
    # 2013-01-19 return list of best ones (can be alternatives)

    # if you want the superposed model use:
    #  superposed_model2=match.get_transformed_model2(
    #     template=other_model)

    if not hasattr(self,'match_dict'):
      self.match_dict={}
      match_list=None
    else:
      match_list=self.match_dict.get(other_model,None)

    if match_list is None:   # need to get it
      from cctbx import euclidean_model_matching as emma
      test_list=[other_model]
      match_list=[]
      if specifically_test_inverse:
        from copy import deepcopy
        inv_other_model=other_model.change_hand()
        inv_other_model._cb_op=deepcopy(other_model._cb_op)
        test_list.append(inv_other_model)
      for test_model in test_list:
        matches = emma.model_matches(
          model1 = self,
          model2 = test_model,
          tolerance = tolerance,
          models_are_diffraction_index_equivalent = \
              models_are_diffraction_index_equivalent,
          break_if_match_with_no_singles=break_if_match_with_no_singles,
          )

        if (matches.n_matches() > 0):
          best_number_of_matches=len(matches.refined_matches[0].pairs)
          for match in matches.refined_matches:
            n=len(match.pairs)
            if n > 0 and n >= best_number_of_matches:
              match_list.append(match)

      self.match_dict[other_model]=match_list # save it

    if not match_list: match_list=[None]
    return match_list


def filter_shift(continuous_shift_flags, shift, selector=1):
  filtered_shift = [0,0,0]
  for i in range(3):
    if (continuous_shift_flags[i] == selector):
      filtered_shift[i] = shift[i]
  return filtered_shift

class euclidean_match_symmetry(object):

  def __init__(self, space_group_info, use_k2l, use_l2n):
    adopt_init_args(self, locals())
    search_symmetry = sgtbx.search_symmetry(
      flags=sgtbx.search_symmetry_flags(
        use_space_group_symmetry=False,
        use_space_group_ltr=-1,
        use_seminvariants=True,
        use_normalizer_k2l=use_k2l,
        use_normalizer_l2n=use_l2n),
      space_group_type=space_group_info.type(),
      seminvariant=space_group_info.structure_seminvariants())
    self.rt_mx = search_symmetry.subgroup()
    self.continuous_shifts = search_symmetry.continuous_shifts()
    assert search_symmetry.continuous_shifts_are_principal()
    self.continuous_shift_flags = search_symmetry.continuous_shift_flags()

  def filter_shift(self, shift, selector=1):
    return filter_shift(self.continuous_shift_flags, shift, selector)

  def show(self, title="", f=None):
    if (f is None): f = sys.stdout
    print(("euclidean_match_symmetry: " + title).rstrip(), file=f)
    print(self.rt_mx.type().lookup_symbol(), file=f)
    print(self.continuous_shifts, file=f)

def generate_singles(n, i):
  singles = list(range(n))
  del singles[i]
  return singles

def pair_sort_function(pair_a, pair_b):
  # Deprecated. Do not use
  from libtbx.math_utils import cmp
  return cmp(pair_a[0], pair_b[0])

def inside_zero_one(c):
  new_c=[]
  for x in c:
    new_c.append(math.fmod(x+100.,1.0))
  return matrix.col(new_c)

def match_refine_times():
  return dicts.easy(
    exclude_pairs=0,
    add_pairs=0,
    eliminate_weak_pairs=0,
    refine_adjusted_shift=0)

class match_refine(object):

  def __init__(self, tolerance,
               ref_model1, ref_model2,
               match_symmetry,
               add_pair_ext,
               i_pivot1, i_pivot2,
               eucl_symop,
               initial_shift,
               times=None):
    adopt_init_args(self, locals(), exclude=("initial_shift",))
    self.singles1 = generate_singles(self.ref_model1.size(), self.i_pivot1)
    self.singles2 = generate_singles(self.ref_model2.size(), self.i_pivot2)
    self.pairs = [(self.i_pivot1, self.i_pivot2)]
    self.adjusted_shift = matrix.col(initial_shift)
    if (self.times is None):
      self.times = match_refine_times()
    self.exclude_pairs()
    self.add_pairs()
    self.eliminate_weak_pairs()
    self.ref_eucl_rt = sgtbx_rt_mx_as_matrix_rt(self.eucl_symop) \
                     + self.adjusted_shift
    self.pairs.sort(key=operator.itemgetter(0))
    self.singles1.sort()
    self.singles2.sort()
    self.calculate_rms()

  def exclude_pairs(self):
    # exclude all pairs with dist >= 4 * tolerance
    # dist_allowed is invariant under refine_adjusted_shift:
    #   exclude all pairs with dist_allowed >= tolerance
    # if 0 continuous shifts: dist_allowed == dist:
    #   exclude all pairs with dist >= tolerance
    timer = user_plus_sys_time()
    self.add_pair_ext.next_pivot(
      self.match_symmetry.continuous_shift_flags,
      self.eucl_symop,
      self.adjusted_shift,
      flex.int(self.singles1),
      flex.int(self.singles2))
    self.times.exclude_pairs += timer.delta()

  def add_pairs(self):
    # XXX possible optimizations:
    #   if 0 continuous shifts:
    #     tabulate dist
    #     keep all < tolerance
    #     we do not need eliminate_weak_matches
    timer = user_plus_sys_time()
    while (len(self.singles1) and len(self.singles2)):
      if (not self.add_pair_ext.next_pair(
        self.adjusted_shift,
        flex.int(self.singles1),
        flex.int(self.singles2))):
        break
      new_pair = (self.add_pair_ext.new_pair_1(),
                  self.add_pair_ext.new_pair_2())
      self.pairs.append(new_pair)
      self.singles1.remove(new_pair[0])
      self.singles2.remove(new_pair[1])
      self.refine_adjusted_shift()
    self.times.add_pairs += timer.delta()

  def eliminate_weak_pairs(self):
    timer = user_plus_sys_time()
    while 1:
      weak_pair = 0
      max_dist = 0
      for pair in self.pairs[1:]:
        dist = self.calculate_shortest_dist(pair)
        if (dist > max_dist):
          weak_pair = pair
          max_dist = dist
      if (weak_pair == 0): break
      if (max_dist < self.tolerance):
        dist = self.calculate_shortest_dist(self.pairs[0])
        if (dist < self.tolerance):
          break
      assert len(self.pairs) > 1
      self.pairs.remove(weak_pair)
      self.singles1.append(weak_pair[0])
      self.singles2.append(weak_pair[1])
      self.refine_adjusted_shift()
    self.times.eliminate_weak_pairs += timer.delta()

  def apply_eucl_ops(self, i_model2):
    c2 = matrix.col(self.eucl_symop * self.ref_model2[i_model2].site)
    return c2 + self.adjusted_shift

  def calculate_shortest_diff(self, pair):
    c2 = self.apply_eucl_ops(pair[1])
    return sgtbx.min_sym_equiv_distance_info(
      self.add_pair_ext.equiv1(pair[0]), c2).diff()

  def calculate_shortest_dist(self, pair):
    length = self.ref_model1.unit_cell().length
    return length(self.calculate_shortest_diff(pair))

  def calculate_shortest_diffs(self):
    shortest_diffs = []
    for pair in self.pairs:
      shortest_diffs.append(self.calculate_shortest_diff(pair))
    return shortest_diffs

  def refine_adjusted_shift(self):
    timer = user_plus_sys_time()
    unit_cell = self.ref_model1.unit_cell()
    sum_diff_cart = matrix.col([0.,0.,0.])
    for diff in self.calculate_shortest_diffs():
      diff_allowed = self.match_symmetry.filter_shift(diff, selector=1)
      diff_cart = unit_cell.orthogonalize(diff_allowed)
      sum_diff_cart += matrix.col(diff_cart)
    mean_diff_cart = sum_diff_cart / len(self.pairs)
    mean_diff_frac = matrix.col(unit_cell.fractionalize(mean_diff_cart))
    self.adjusted_shift = matrix.col(self.adjusted_shift) + mean_diff_frac
    self.times.refine_adjusted_shift += timer.delta()

  def calculate_rms(self):
    length = self.ref_model1.unit_cell().length
    self.shortest_distances = flex.double([
      length(d) for d in self.calculate_shortest_diffs() ])
    self.rms = math.sqrt(flex.sum_sq(self.shortest_distances)/len(self.pairs))

  def show(self, f=None, truncate_singles=None, singles_per_line=5):
    if (f is None): f = sys.stdout
    print("Match summary:", file=f)
    print("  Operator:", file=f)
    print("       rotation:", self.rt.r.mathematica_form(format="%.6g"), file=f)
    print("    translation:", \
      self.rt.t.transpose().mathematica_form(format="%.6g")[1:-1], file=f)
    print("  rms coordinate differences: %.2f" % (self.rms,), file=f)
    print("  Pairs:", len(self.pairs), file=f)
    for pair in self.pairs:
      print("   ", self.ref_model1[pair[0]].label, end=' ', file=f)
      print(self.ref_model2[pair[1]].label, end=' ', file=f)
      print("%.3f" % (self.calculate_shortest_dist(pair),), file=f)
    for i_model,ref_model,singles in ((1,self.ref_model1,self.singles1),
                                      (2,self.ref_model2,self.singles2)):
      print("  Singles model %s:" % i_model, len(singles), end=' ', file=f)
      i = 0
      for s in singles:
        if (i == truncate_singles): break
        if (i % singles_per_line == 0):
          print(file=f)
          print(" ", end=' ', file=f)
        print(" ", ref_model[s].label, end=' ', file=f)
        i += 1
      print(file=f)
    print(file=f)

  def get_transformed_model2(self,output_pdb=None,
    scattering_type="SE",f=sys.stdout,
    return_superposed_model2=True,template_pdb_inp=None):
      # tt 2013-01-25; 2016-10-31
      from cctbx import xray
      xray_scatterer = xray.scatterer( scattering_type = scattering_type)
      model2=self.ref_model2.as_xray_structure(xray_scatterer)
      from cctbx.array_family import flex
      new_coords=flex.vec3_double()
      for i_model2 in range(self.ref_model2.size()):
        c2 = matrix.col(self.eucl_symop * self.ref_model2[i_model2].site)
        c2 += self.adjusted_shift
        c2=inside_zero_one(c2)
        new_coords.append(c2)
      model2.set_sites_frac(new_coords)


      if output_pdb is not None:
        assert template_pdb_inp is not None
        # Set up new xrs with these sites and with scattering types, occ, b,
        #   labels from original 2nd model
        xrs=xray.structure(model2.xray_structure())
        assert len(model2.scatterers())==len(template_pdb_inp.atoms())
        b_iso_values=flex.double()
        for scatterer,atom in zip(model2.scatterers(),template_pdb_inp.atoms()):
          b_iso_values.append(atom.b)
          new_scatterer = xray.scatterer(
            scattering_type = atom.element,
            label=atom.name,
            occupancy=atom.occ,
            site=scatterer.site)
          xrs.add_scatterer(new_scatterer)
        xrs.set_b_iso(values = b_iso_values)
        pdb_string=xrs.as_pdb_file()
        ff=open(output_pdb,'w')
        print(pdb_string, file=ff)
        ff.close()
        print("\nWrote model 2 mapped to model 1 to file %s " %(output_pdb), file=f)

      if return_superposed_model2:
        return model2.as_emma_model()

def match_sort_function(match_a, match_b):
  # Deprecated. Do not use
  from libtbx.math_utils import cmp
  i = -cmp(len(match_a.pairs), len(match_b.pairs))
  if (i): return i
  return cmp(match_a.rms, match_b.rms)

def weed_refined_matches(space_group_number, refined_matches,
                         rms_penalty_per_site):
  n_matches = len(refined_matches)
  if (n_matches == 0): return
  best_rms = refined_matches[0].rms
  best_n_pairs = len(refined_matches[0].pairs)
  is_redundant = [0] * n_matches
  for i in range(n_matches-1):
    match_i = refined_matches[i]
    if (is_redundant[i]): continue
    if (match_i.rms < best_rms):
      best_rms = match_i.rms
      best_n_pairs = len(match_i.pairs)
    for j in range(i+1, n_matches):
      match_j = refined_matches[j]
      if (   match_i.pairs == match_j.pairs
          or (    rms_penalty_per_site
              and match_j.rms > best_rms * (1 - rms_penalty_per_site * (
                    best_n_pairs - len(match_j.pairs))))):
        is_redundant[j] = 1
  for i in range(n_matches-1, -1, -1):
    if (is_redundant[i]):
      del refined_matches[i]
  if (space_group_number == 1 and n_matches > 0):
    trivial_matches_only = True
    for match in refined_matches:
      if (len(match.pairs) > 1):
        trivial_matches_only = False
        break
    if (trivial_matches_only):
      while (    len(refined_matches) > 1
             and refined_matches[0].pairs[0] != (0,0)): del refined_matches[0]
      while (len(refined_matches) > 1): del refined_matches[-1]
    else:
      while (len(refined_matches[-1].pairs) == 1): del refined_matches[-1]

def match_rt_from_ref_eucl_rt(model1_cb_op, model2_cb_op, ref_eucl_rt):
  inv_m1 = sgtbx_rt_mx_as_matrix_rt(model1_cb_op.c_inv())
  m2 = sgtbx_rt_mx_as_matrix_rt(model2_cb_op.c())
  # X2_orig -> X2_ref -> X1_ref -> X1_orig
  #          m2   ref_eucl_rt   inv_m1
  # X1 = inv_m1 * ref_eucl_rt * m2 * X2
  return inv_m1 * ref_eucl_rt * m2

def compute_refined_matches(ref_model1, ref_model2,
                            tolerance,
                            models_are_diffraction_index_equivalent,
                            shall_break):
  match_symmetry = euclidean_match_symmetry(
    ref_model1.space_group_info(),
    use_k2l=True, use_l2n=(not models_are_diffraction_index_equivalent))
  ref_model1_sites = flex.vec3_double([pos.site for pos in ref_model1])
  ref_model2_sites = flex.vec3_double([pos.site for pos in ref_model2])
  add_pair_ext = ext.add_pair(
    tolerance,
    ref_model1.unit_cell(),
    ref_model1.space_group(),
    ref_model1.min_distance_sym_equiv(),
    ref_model1_sites,
    ref_model2_sites)
  accumulated_match_refine_times = match_refine_times()
  refined_matches = []
  for i_pivot1 in range(ref_model1.size()):
    for i_pivot2 in range(ref_model2.size()):
      for eucl_symop in match_symmetry.rt_mx:
        c2 = eucl_symop * ref_model2[i_pivot2].site
        dist_info = sgtbx.min_sym_equiv_distance_info(
          add_pair_ext.equiv1(i_pivot1),
          c2,
          match_symmetry.continuous_shift_flags)
        if (dist_info.dist() < tolerance):
          allowed_shift = dist_info.continuous_shifts()
          match = match_refine(tolerance,
                               ref_model1, ref_model2,
                               match_symmetry,
                               add_pair_ext,
                               i_pivot1, i_pivot2,
                               eucl_symop,
                               allowed_shift,
                               accumulated_match_refine_times)
          match.rt = match_rt_from_ref_eucl_rt(
            ref_model1.cb_op(),
            ref_model2.cb_op(),
            match.ref_eucl_rt)
          refined_matches.append(match)
          if shall_break(match):
            return refined_matches
  #print accumulated_match_refine_times
  return refined_matches


class delegating_model_matches(object):
  """
  emma loops over all pairs of sites from the first and second structure.
  If a pair is closer than tolerance (under Euclidean symmetry) it
  initiates a "match and refine" procedure which incrementally adds the
  "next closest" pair, with a distance below 2*tolerance. Each time
  a pair is added the allowed origin shifts (if any) are refined to
  minimize the rmsd of all the pairs matched so far. When there are
  no more pairs to add, emma goes into a "weeding" procedure. This
  procedure sorts the pairs by distance. If there are distances above
  tolerance, the worst pair is removed, followed by the same allowed
  origin shift refinement as before. The weeding is repeated until
  there are no pairs above tolerance.

  emma was written with substructures in mind, which usually have a few
  (by macromolecular standards) sites placed far apart in space.
  """

  def __init__(self, model1, model2,
                     tolerance=1.,
                     models_are_diffraction_index_equivalent=False,
                     shall_break=None,
                     rms_penalty_per_site=0.05):
    if shall_break is None:
      def shall_break(match): return False
    adopt_init_args(self, locals())
    assert model1.cb_op().is_identity_op()
    assert model2.cb_op().is_identity_op()
    ref_model1 = model1.transform_to_reference_setting()
    ref_model2 = model2.transform_to_reference_setting()
    assert ref_model1.unit_cell().is_similar_to(ref_model2.unit_cell())
    if (ref_model1.space_group() != ref_model2.space_group()):
      ref_model2 = ref_model2.change_hand()
    assert ref_model1.unit_cell().is_similar_to(ref_model2.unit_cell())
    assert ref_model1.space_group() == ref_model2.space_group()
    self.refined_matches = compute_refined_matches(
      ref_model1, ref_model2,
      tolerance,
      models_are_diffraction_index_equivalent,
      shall_break)
    self.refined_matches.sort(key=lambda element: (-len(element.pairs), element.rms))
    weed_refined_matches(model1.space_group_info().type().number(),
                         self.refined_matches, rms_penalty_per_site)


  def n_matches(self):
    return len(self.refined_matches)

  def n_pairs_best_match(self):
    if (len(self.refined_matches) == 0): return 0
    return len(self.refined_matches[0].pairs)

  def consensus_model(self, i_model=1, i_refined_matches=0):
    assert i_model in (1,2)
    assert 0 <= i_refined_matches < self.n_matches()
    if (i_model == 1):
      source_model = self.model1
    else:
      source_model = self.model2
    if (self.n_matches() == source_model.size()):
      return source_model
    result = model(special_position_settings=source_model)
    i_model -= 1
    for pair in self.refined_matches[i_refined_matches].pairs:
      result.add_position(source_model.positions()[pair[i_model]])
    return result

  def transform_model(self, i_model, i_refined_matches=0):
    assert i_model in (1,2)
    assert 0 <= i_refined_matches < self.n_matches()
    rt = self.refined_matches[i_refined_matches].rt
    if (i_model == 1):
      result = model(special_position_settings=self.model2)
      source_model = self.model1
      rt = rt.inverse()
    else:
      result = model(special_position_settings=self.model1)
      source_model = self.model2
    for pos in source_model.positions():
      result.add_position(position(label=pos.label, site=(rt*pos.site).elems))
    return result


class model_matches(delegating_model_matches):

  def __init__(self, model1, model2,
                     tolerance=1.,
                     models_are_diffraction_index_equivalent=False,
                     break_if_match_with_no_singles=False,
                     rms_penalty_per_site=0.05):
    if break_if_match_with_no_singles:
      def shall_break(match):
        return len(match.singles1) == 0 or len(match.singles2) == 0
    else:
      shall_break = None
    super(model_matches, self).__init__(
      model1, model2,
      tolerance,
      models_are_diffraction_index_equivalent,
      shall_break,
      rms_penalty_per_site)
