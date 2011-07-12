from cctbx import sgtbx
from scitbx.array_family import flex
import sys

class partition_t(list): pass

class left_decomposition(object):
  def __init__(self, g, h):
    g_lattice_translations = [ g(i,0,0).mod_short() for i in xrange(g.n_ltr())]

    if self.is_subgroup(g,h):
      self.h_name = str( sgtbx.space_group_info( group = h ) )
      self.g_name = str( sgtbx.space_group_info( group = g ) )
      g = [s for s in g] # for speed, convert to plain Python list
      h = [s for s in h]
      assert len(g) % len(h) == 0
      assert h[0].is_unit_mx()
      self.partition_indices = [-1] * len(g)
      self.partitions = []
      for i,gi in enumerate(g):
        if (self.partition_indices[i] != -1): continue
        self.partition_indices[i] = len(self.partitions)
        partition = partition_t([gi])
        for hj in h[1:]:
          gihj = gi.multiply(hj)
          for k in xrange(i+1,len(g)):
            if (self.partition_indices[k] != -1): continue
            gk = g[k]
            tmp_g = gk.inverse().multiply( gihj ).mod_short()
            if tmp_g.as_xyz() in [ltr.as_xyz() for ltr in g_lattice_translations]:
              self.partition_indices[k] = len(self.partitions)
              partition.append(gk)
              break
          else:
            raise RuntimeError("h is not a subgroup of g")
        if (len(partition) != len(h)):
          raise RuntimeError("h is not a subgroup of g")
        self.partitions.append(partition)
      if (len(self.partitions) * len(h) != len(g)):
        raise RuntimeError("h is not a subgroup of g")
      #sort cosets by operator order
      self.sort_cosets()
    else:
      raise RuntimeError("h is not a subgroup of g")

  def is_subgroup(self, g, h):
    tst_group = str( sgtbx.space_group_info(group=g) )
    tst_group = sgtbx.space_group_info( tst_group ).group()
    h = [s for s in h]
    for s in h:
      tst_group.expand_smx( s )
    if str(sgtbx.space_group_info(group=tst_group)) == str( sgtbx.space_group_info(group=g) ):
      return True
    else:
      return False

  def sort_cosets(self):
    new_partitions = []
    for pp in self.partitions:
      orders = []
      for op in pp:
        orders.append( op.r().info().type() )
      orders = flex.sort_permutation( flex.int(orders) )
      tmp = partition_t()
      for ii in orders:
        tmp.append( pp[ii] )
      new_partitions.append( tmp )
    self.partitions = new_partitions

  def show(self,out=None, cb_op=None, format="cosets_form"):
    from cctbx.sgtbx.literal_description import literal_description
    if out is None:
      out = sys.stdout
    count=0
    group_g = self.g_name
    group_h = self.h_name
    if cb_op is not None:
      group_g = str( sgtbx.space_group_info( self.g_name ).change_basis( cb_op ) )
      group_h = str( sgtbx.space_group_info( self.h_name ).change_basis( cb_op ) )

    print >> out, "Left cosets of :"
    print >> out, "  subgroup  H: %s"%( group_h )
    print >> out, "  and group G: %s"%( group_g )
    for part in self.partitions:
      extra_txt="   (all operators from H)"

      tmp_group = sgtbx.space_group_info( self.h_name ).group()
      tmp_group.expand_smx( part[0] )
      if cb_op is None:
        tmp_group = sgtbx.space_group_info( group=tmp_group)
      else:
        tmp_group = sgtbx.space_group_info( group=tmp_group).change_basis(cb_op)

      if count>0:
        extra_txt = "   (H+coset[%i] = %s)"%(count,tmp_group)
      print >> out
      print >> out, "  Coset number : %5s%s"%(count, extra_txt)
      print >> out
      count += 1
      for item in part:
        tmp_item = None
        if cb_op is None:
          tmp_item = item
        else:
          tmp_item = cb_op.apply( item )

        print >> out, literal_description(tmp_item).select(format)


class left_decomposition_point_groups_only(object):

  def __init__(self, g, h):
    def convert_to_plain_list(group): # for speed
      assert group.n_ltr() == 1
      result = []
      for s in group:
        assert s.t().is_zero()
        result.append(s)
      return result
    g = convert_to_plain_list(group=g)
    h = convert_to_plain_list(group=h)
    assert len(g) % len(h) == 0
    assert h[0].is_unit_mx()
    self.partition_indices = [-1] * len(g)
    self.partitions = []
    for i,gi in enumerate(g):
      if (self.partition_indices[i] != -1): continue
      self.partition_indices[i] = len(self.partitions)
      partition = [gi]
      for hj in h[1:]:
        gihj = gi.multiply(hj)
        for k in xrange(i+1,len(g)):
          if (self.partition_indices[k] != -1): continue
          gk = g[k]
          if (gk.r().num() == gihj.r().num()):
            self.partition_indices[k] = len(self.partitions)
            partition.append(gk)
            break
        else:
          raise RuntimeError("h is not a subgroup of g")
      if (len(partition) != len(h)):
        raise RuntimeError("h is not a subgroup of g")
      self.partitions.append(partition)
    if (len(self.partitions) * len(h) != len(g)):
      raise RuntimeError("h is not a subgroup of g")

  def best_partition_representatives(self,
        cb_op=None,
        omit_first_partition=False,
        omit_negative_determinants=False):
    result = []
    for partition in self.partitions:
      best_choice, best_choice_as_hkl = None, None
      for choice in partition:
        if (omit_first_partition and choice.is_unit_mx()):
          assert best_choice is None
          break
        if (omit_negative_determinants and choice.r().determinant() < 0):
          continue
        if (cb_op is not None):
          choice = cb_op.apply(choice)
        choice_as_hkl = choice.r().as_hkl()
        if (best_choice_as_hkl is None
            or sgtbx.compare_cb_op_as_hkl(
                 best_choice_as_hkl, choice_as_hkl) > 0):
          best_choice, best_choice_as_hkl = choice, choice_as_hkl
      if (best_choice is not None):
        result.append(best_choice)
    return result

def double_unique(g, h1, h2):
  """g is the supergroup
     h1 and h2 are subgroups
  """
  # Make lists of symops for all groups
  g = [s for s in g]
  h1 = [s for s in h1]
  h2 = [s for s in h2]
  # this is our final result
  result = []
  # This set keeps track of equivalent symops
  done = set()
  #
  for a in g:
    if (str( a ) in done): continue
    result.append(a)
    for hi in h1:
      for hj in h2:
        b = hi.multiply(a).multiply(hj)
        done.add(str( b ))
  return result

def construct_nice_cb_op(coset,
                         sym_transform_1_to_2,
                         to_niggli_1,
                         to_niggli_2):
  best_choice = None
  best_choice_as_hkl = None
  to_niggli_2 = to_niggli_2.new_denominators( to_niggli_1 )
  sym_transform_1_to_2 = sym_transform_1_to_2.new_denominators( to_niggli_1 )

  for coset_element in coset:
    tmp_coset_element = sgtbx.change_of_basis_op(coset_element)
    tmp_coset_element = tmp_coset_element.new_denominators( to_niggli_1 )
    tmp_op = to_niggli_1.inverse() \
           * (tmp_coset_element * sym_transform_1_to_2) \
           * to_niggli_2
    if (best_choice_as_hkl is None
        or sgtbx.compare_cb_op_as_hkl(
             best_choice_as_hkl, tmp_op.as_hkl()) > 0):
      best_choice = tmp_op
      best_choice_as_hkl =  tmp_op.as_hkl()
  assert best_choice is not None

  tmptmp = sgtbx.change_of_basis_op(coset[0])

  return best_choice

class double_cosets(object):
  def __init__(self,g, h1, h2, enforce_det_ge_1=True):
    """g is the supergroup
       h1 and h2 are subgroups
    """
    # Make lists of symops for all groups
    g = [s for s in g]
    h1 = [s for s in h1]
    h2 = [s for s in h2]

    # a list of lists with our double cosets
    self.double_cosets = []

    #
    for a in g:
      # first we have to check whether or not
      # this symmetry operator is allready in a coset we
      # might have constructured earlier
      if not self.is_in_list_of_cosets( a ):
        # not present, make de double coset please
        tmp_double_coset = []
        tmp_double_coset.append( a )
        # The other members will now be made
        for hi in h1:
          for hj in h2:
            b = hi.multiply(a).multiply(hj)
            #check if this element is allready in this coset please
            if not self.is_in_coset( b, tmp_double_coset ):
              tmp_double_coset.append( b )
        self.double_cosets.append( tmp_double_coset )
    if enforce_det_ge_1:
      self.clear_up_cosets()

  def clear_up_cosets(self):
    temp_cosets = []
    for cs in self.double_cosets:
      if cs[0].r().determinant() > 0:
        temp_cosets.append( cs )
    self.double_cosets = temp_cosets

  def is_in_coset(self, a, coset_list):
    found_it=False
    for hi in coset_list:
      if ( str(hi.mod_positive() ) == str(a.mod_positive() ) ):
        found_it = True
        break
    return found_it

  def is_in_list_of_cosets( self, a ):
    found_it = False
    for cs in self.double_cosets:
      if self.is_in_coset( a, cs ):
        found_it = True
    return found_it

  def have_duplicates(self):
    n_cosets = len(self.double_cosets)
    for ics in xrange(n_cosets):
      tmp_cs = self.double_cosets[ics]
      for jcs in xrange(n_cosets):
        if ics != jcs :
          tmp_cs_2 = self.double_cosets[jcs]
          # now check each element of tmp_cs
          for hi in tmp_cs:
            if (self.is_in_coset(hi, tmp_cs_2)): return True
    return False

  def show(self,out=None):
    if out == None:
      out = sys.stdout
    print >> out, "The double cosets are listed below"
    for cs in self.double_cosets:
      for a in cs:
        print >> out, "("+str(a)+")    ",
      print >> out

def test_double_coset_decomposition():
  from  cctbx.sgtbx import subgroups
  for space_group_number in xrange(17,44):
    parent_group_info = sgtbx.space_group_info(space_group_number)
    subgrs = subgroups.subgroups(parent_group_info).groups_parent_setting()
    g = parent_group_info.group()
    for h1 in subgrs:
      for h2 in subgrs:
        tmp_new = double_cosets(g, h1, h2)
        assert not tmp_new.have_duplicates()

def test_lattice_translation_aware_left_decomposition():
  from cctbx import crystal
  def generate_cases():
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'C 2 2 21' ),
                                      unit_cell = (74.033, 96.747, 109.085, 90, 90, 90) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'C 2 1 1' ),
                                      unit_cell = (74.033, 96.747, 109.085, 89.98, 90, 90) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'F 2 3' ),
                                      unit_cell = (124.059, 124.059, 124.059, 90, 90, 90) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'R 3 :H (-1/2*x+z,1/2*x-1/2*y+z,1/2*y+z)' ),
                                      unit_cell = (124.059, 124.059, 124.059, 90.0003, 90.0003, 90.0003) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 6 2 2' ),
                                      unit_cell = (94.705, 94.705, 57.381, 90, 90, 120) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 6' ),
                                      unit_cell = (94.705, 94.705, 57.381, 90, 90, 120) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'I 2 3' ),
                                      unit_cell = (80.8965, 80.8965, 80.8965, 90, 90, 90) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'R 3 :H (y+1/2*z,-x+1/2*z,x-y+1/2*z)' ),
                                      unit_cell = (80.8965, 80.8965, 80.8965, 89.993, 89.993, 89.993) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 61 2 2' ),
                                      unit_cell = (115.7, 115.7, 247.2, 90, 90, 120) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 31 2 1 (a,b,c-1/6)' ),
                                      unit_cell = (115.7, 115.7, 247.2, 90, 90, 120) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 41 21 2' ),
                                      unit_cell = (64.766, 64.766, 270.651, 90, 90, 90) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 41 (a,b+1/2,c)' ),
                                      unit_cell = (64.766, 64.766, 270.651, 90, 90, 90) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'I 4 2 2' ),
                                      unit_cell = (124.282, 124.282, 117.485, 90, 90, 90) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'I 2 1 1' ),
                                      unit_cell = (124.265, 124.299, 117.485, 89.9919, 90, 90) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 43 21 2' ),
                                      unit_cell = (100.228, 100.228, 294.099, 90, 90, 90) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'C 2 2 21 (x-y,x+y,z)' ),
                                      unit_cell = (100.228, 100.228, 294.099, 90, 90, 90.0036) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 43 3 2' ),
                                      unit_cell = (170.4, 170.4, 170.4, 90, 90, 90) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 21 3' ),
                                      unit_cell = (170.4, 170.4, 170.4, 90, 90, 90) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'R 3 2 :H' ),
                                      unit_cell = (221.819, 221.819, 55.601, 90, 90, 120) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'C 1 2 1 (3/2*a-1/2*b+c,b,c)' ),
                                      unit_cell = (221.835, 221.787, 55.601, 90, 89.9905, 119.993) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 64 2 2' ),
                                      unit_cell = (142.63, 142.63, 108.35, 90, 90, 120) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 31' ),
                                      unit_cell = (142.63, 142.63, 108.35, 90, 90, 120) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'I 41 3 2' ),
                                      unit_cell = (134.242, 134.242, 134.242, 90, 90, 90) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'I 21 3' ),
                                      unit_cell = (134.242, 134.242, 134.242, 90, 90, 90) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 41 21 2' ),
                                      unit_cell = (86.8055, 86.8055, 75.688, 90, 90, 90) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 21 21 21 (a+1/4,b,c-3/8)' ),
                                      unit_cell = (86.782, 86.829, 75.688, 90, 90, 90) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 31 2 1' ),
                                      unit_cell = (97.868, 97.868, 208.953, 90, 90, 120) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'C 1 2 1 (x+y,-x+y,z)' ),
                                      unit_cell = (97.8535, 97.8535, 208.953, 89.9913, 90.0087, 119.971) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'I 4 3 2' ),
                                      unit_cell = (121.805, 121.805, 121.805, 90, 90, 90) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'I 2 3' ),
                                      unit_cell = (121.805, 121.805, 121.805, 90, 90, 90) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 43 2 2' ),
                                      unit_cell = (67.472, 67.472, 174.41, 90, 90, 90) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 2 2 21 (a,b,c-1/4)' ),
                                      unit_cell = (67.479, 67.465, 174.41, 90, 90, 90) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'C 2 2 21' ),
                                      unit_cell = (47.623, 120.166, 57.648, 90, 90, 90) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 1 (a+b,a-b,-c)' ),
                                      unit_cell = (47.623, 120.166, 57.648, 90.0124, 89.86, 89.9641) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'I 2 2 2' ),
                                      unit_cell = (57.816, 75.599, 155.949, 90, 90, 90) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 1 (b+c,a+c,a+b)' ),
                                      unit_cell = (57.816, 75.599, 155.949, 89.9824, 90.0357, 90.02) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 41 21 2' ),
                                      unit_cell = (97.6539, 97.6539, 215.352, 90, 90, 90) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'C 1 2 1 (x-y,x+y,z-1/4)' ),
                                      unit_cell = (97.6539, 97.6539, 215.352, 90, 90, 90.022) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 62 2 2' ),
                                      unit_cell = (58.306, 58.306, 137.882, 90, 90, 120) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 32 1 2 (a,b,c-1/6)' ),
                                      unit_cell = (58.306, 58.306, 137.882, 90, 90, 120) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 3 2 1' ),
                                      unit_cell = (56.2477, 56.2477, 158.706, 90, 90, 120) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'C 1 2 1 (x+y,-x+y,z)' ),
                                      unit_cell = (56.2235, 56.2235, 158.706, 89.9654, 90.0346, 119.915) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'F 2 2 2' ),
                                      unit_cell = (88.492, 90.284, 90.6894, 90, 90, 90) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'C 1 2 1 (a,b,a+2*c)' ),
                                      unit_cell = (88.492, 90.284, 90.6894, 90, 89.9856, 90) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'I 41' ),
                                      unit_cell = (59.2061, 59.2061, 59.48, 90, 90, 90) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 1 (b+c,a+c,a+b)' ),
                                      unit_cell = (59.1874, 59.2247, 59.48, 90.0179, 89.956, 89.9542) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 3 1 2' ),
                                      unit_cell = (56.9, 56.9, 62.77, 90, 90, 120) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 3' ),
                                      unit_cell = (56.9, 56.9, 62.77, 90, 90, 120) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 43' ),
                                      unit_cell = (54.111, 54.111, 71.637, 90, 90, 90) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 1 1 21' ),
                                      unit_cell = (54.09, 54.132, 71.637, 90, 90, 90.04) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 4 3 2' ),
                                      unit_cell = (122.638, 122.638, 122.638, 90, 90, 90) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 2 3' ),
                                      unit_cell = (122.638, 122.638, 122.638, 90, 90, 90) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 62' ),
                                      unit_cell = (115.246, 115.246, 67.375, 90, 90, 120) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 1 1 2' ),
                                      unit_cell = (115.246, 115.242, 67.375, 90, 90, 119.997) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 63' ),
                                      unit_cell = (103.371, 103.371, 91.38, 90, 90, 120) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 1 1 21' ),
                                      unit_cell = (103.35, 103.36, 91.38, 90, 90, 119.97) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'P 65 2 2' ),
                                      unit_cell = (183.226, 183.226, 141.117, 90, 90, 120) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'C 1 2 1 (x-y,2*x,z)' ),
                                      unit_cell = (183.108, 183.285, 141.117, 90.0433, 90, 119.968) )}
    yield { 'group':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'F 41 3 2' ),
                                      unit_cell = (228.763, 228.763, 228.763, 90, 90, 90) ),
         'subgroup':crystal.symmetry( space_group_info = sgtbx.space_group_info( 'I 41 (c,a+b,-a+b)' ),
                                      unit_cell = (228.574, 228.858, 228.858, 90, 90, 90) )}

  for case in generate_cases():
    C = left_decomposition(g=case['group'].space_group(), h = case['subgroup'].space_group())

    # double check that the first partition is identical to the subgroup:
    for element in case['subgroup'].space_group():
      assert element in C.partitions[0]

    # make sure each coset is equal to the product of the coset representative and the subgroup
    for x in xrange(1,len(C.partitions)):
      part = C.partitions[x]
      coset_representative = part[0]
      for element in case['subgroup'].space_group():
        assert coset_representative.multiply(element).mod_short().as_xyz() in [
          i.mod_short().as_xyz() for i in part]

def run():
  test_double_coset_decomposition()
  test_lattice_translation_aware_left_decomposition()
  print "OK"

if (__name__ == "__main__"):
  run()
