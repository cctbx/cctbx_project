from cctbx import sgtbx
import sys

class left_decomposition(object):

  def __init__(self, g, h):
    g = [s for s in g] # for speed, convert to plain Python list
    h = [s for s in h]
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
  # This dictionary keeps track of equivalent symops
  done = {}
  #
  for a in g:
    if (str( a ) in done): continue
    result.append(a)
    for hi in h1:
      for hj in h2:
        b = hi.multiply(a).multiply(hj)
        done[str( b )] = None
  return result

def compare_cb_op_as_hkl(a, b):
  if (len(a) < len(b)): return -1
  if (len(a) > len(b)): return  1
  return cmp(a, b)

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
    tmp_op = to_niggli_1.inverse() * (tmp_coset_element * sym_transform_1_to_2) * to_niggli_2
    if ( (best_choice_as_hkl is None) or
         (compare_cb_op_as_hkl( best_choice_as_hkl, tmp_op.as_hkl() ) >0 ) ):
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
      # might have cnostructured earlier
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

def run():
  test_double_coset_decomposition()
  print "OK"

if (__name__ == "__main__"):
  run()
