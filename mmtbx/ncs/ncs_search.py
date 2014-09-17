from __future__ import division
from scitbx.array_family import flex
from libtbx.utils import Sorry

def simple_alignment(seq_a, seq_b, similarity=0.9):
  """
  Do basic alignment of seq_a and seq_b using dynamic programing

  Give penalty for misalignment and gaps. Do not give any points for
  alignment. (score only change when "bad" things happens)
  Stop if the alignment is less than the "similarity"

  Args:
    seq_a, seq_b (lists of str): list of characters to compare
    similarity (float): minimum boundary on how good should be the alignment

  Returns:
    aligned_sel_a (flex.size_t): the indices of the aligned components of seq_a
    aligned_sel_b (flex.size_t): the indices of the aligned components of seq_b
  """
  a = len(seq_a)
  b = len(seq_b)
  gap_penalty = 1
  misalign_penalty = 2 * max(a,b)
  # limit the number of mis-alignments
  n_max = int((1 - similarity) * min(a,b))
  # Starting score according to the required similarity
  # Fixme: check what the max score should be
  score = n_max + 100
  # build score matrix
  R = build_matrix(row=a,col=b,max_score=score,gap_penalty=gap_penalty)
  # populate score matrix
  for j in xrange(1,b + 1):
    for i in xrange(1,a + 1):
      not_aligned = (seq_a[i-1].lower() != seq_b[j-1].lower())
      s1 = R[i-1][j-1] - misalign_penalty * not_aligned
      s2 = R[i-1][j] - gap_penalty
      s3 = R[i][j-1] - gap_penalty
      s = max(s1,s2,s3)
      R[i][j] = s
  aligned_sel_a, aligned_sel_b = get_matching_indices(
    R=R,row=a,col=b,seq_a=seq_a,seq_b=seq_b,gap_penalty=gap_penalty)
  if aligned_sel_a and aligned_sel_b:
    assert aligned_sel_a.size() == aligned_sel_b.size()
  if (a - aligned_sel_a.size()) > n_max:
    return flex.size_t([]), flex.size_t([])
  return aligned_sel_a, aligned_sel_b

def build_matrix(row, col, max_score=None,gap_penalty=1):
  """
  Build a matrix, a list of lists for use in "simple_alignment"
  The matrix values are initialized according the the gap penalties.
  Since "simple_alignment" halts when two sequences are to different (when
  the score is zero) only the allowed positions are being filled

  Args:
    row (int): number of rows
    col (int): number of columns
    max_score (int): the score assigned for perfect alignment
    gap_penalty (int): penalty for gaps

  Returns:
    R (list of lists of int): the score matrix
  """
  R = []
  for i in xrange(row + 1):
    R.append([0] * (col + 1))
  # populate zero cases
  for i in xrange(row + 1):
    R[i][0] = max_score - (i * gap_penalty)
  for i in xrange(col + 1):
    R[0][i] = max_score - (i * gap_penalty)
  return R

def get_matching_indices(R,row,col,seq_a,seq_b,gap_penalty):
  """
  Trace back best solution and collect the indices of matching pairs
  (Note that mis-alignment is not allowed

  Args:
    R (list of lists of int): the score matrix
    row (int): number of rows
    col (int): number of columns
    seq_a, seq_b (lists of str): list of characters to compare
    gap_penalty (int): penalty for gaps

  Returns:
    sel_a (flex.size_t): matching indices sequence a
    sel_b (flex.size_t): matching indices sequence b
  """
  # index of best score from last row
  j_max = R[row].index(max(R[row]))
  i_max = row
  sel_a = []
  sel_b = []
  assert (len(seq_a) == row) and (len(seq_b) == col)
  # Fixme:  does not work if max score is not properly assign
  if R[i_max][j_max] <= 0:
    # best alignment is not good
    return flex.size_t(sel_a), flex.size_t(sel_b)
  else:
    while (i_max > 0) and (j_max > 0):
      aligned = (seq_a[i_max-1].lower() == seq_b[j_max-1].lower())
      if (R[i_max][j_max] == R[i_max - 1][j_max - 1]) and aligned:
        sel_a.append(i_max - 1)
        sel_b.append(j_max - 1)
        i_max -= 1
        j_max -= 1
      elif R[i_max][j_max] == R[i_max - 1][j_max] - gap_penalty:
        i_max -= 1
      elif R[i_max][j_max] == R[i_max][j_max - 1] - gap_penalty:
        j_max -= 1
      else:
        raise Sorry('Sequence alignment error')
    sel_a.reverse()
    sel_b.reverse()
  return flex.size_t(sel_a), flex.size_t(sel_b)

def get_residue_sequence(chain_ph):
  """
  Get a list of pairs (residue name, residue number) and the chain ID

  Args:
    chain_ph (iotbx_pdb_hierarchy_ext): hierarchy of a single chain

  Returns:
    chain_id (str)
    res_seq (list)
  """
  res_list = [x.resname for x in chain_ph.atom_groups()]
  resid_list = [x.resseq for x in chain_ph.residue_groups()]
  res_seq = [(x,y) for (x,y) in zip(res_list, resid_list) if x.lower() != 'hoh']
  #
  atoms = chain_ph.atoms()
  chain_id = atoms[0].chain().id
  return chain_id, res_seq

def align_residues(hierarchy_a, hierarchy_b, similarity=0.9):
  """
  Get the selections of common atoms atoms and selections of atoms that are
  only present in either hierarchy_a or hierarchy_b.

  Args:
    hierarchy_a (hierarchy object)
    hierarchy_b (hierarchy object)
    similarity (float): minimum boundary on how good should be the alignment

  Returns:
    sel_a (flex.size_t): selection of common atoms in chain a
    sel_b (flex.size_t): selection of common atoms in chain b
    not_sel_a (flex.size_t): selection of atoms no longer included in ncs
    not_sel_b (flex.size_t): selection of atoms no longer included in ncs
    chain_id_a, chain_id_b (str)
  """
  chain_id_a, res_a = get_residue_sequence(hierarchy_a)
  chain_id_b, res_b = get_residue_sequence(hierarchy_b)
  # get residue lists
  seq_a = [x for (x,y) in res_a]
  seq_b = [x for (x,y) in res_b]
  res_sel_a, res_sel_b = simple_alignment(seq_a, seq_b, similarity=similarity)

  res_num_a = [res_a[i][1].strip() for i in res_sel_a]
  res_num_b = [res_b[i][1].strip() for i in res_sel_b]

  # collect all atoms that present in both chains
  atom_cache_a = hierarchy_a.atom_selection_cache().selection
  atom_cache_b = hierarchy_b.atom_selection_cache().selection
  sel_a = flex.bool([False]*hierarchy_a.atoms().size())
  sel_b = flex.bool([False]*hierarchy_b.atoms().size())
  for i,(x,y) in enumerate(zip(res_num_a,res_num_b)):
    sa = atom_cache_a('resseq ' + x)
    sb = atom_cache_b('resseq ' + y)
    set_a = set(sa.iselection())
    set_b = set(sb.iselection())
    d_res_num = min(set_a) - min(set_b)
    # temporarily adjust atoms numbers to count for shifts
    set_b = {(x + d_res_num) for x in set_b}
    # update selections
    sel_a,sel_b = update_with_common_atoms(sel_a,set_a,sel_b,set_b,d_res_num)
  not_sel_a = ~sel_a
  not_sel_b = ~sel_b
  sel_a = sel_a.iselection()
  not_sel_a = not_sel_a.iselection()
  sel_b = sel_b.iselection()
  not_sel_b = not_sel_b.iselection()
  return sel_a, sel_b, not_sel_a, not_sel_b, chain_id_a, chain_id_b

def update_with_common_atoms(sel_a,set_a,sel_b,set_b,d_res_num):
  """
  update selection (flex.bool) according to the size of the chain (set_x)
  and whether or not atoms should be included in selection

  Args:
    sel_a,sel_b (flex.bool): selected atoms
    set_a,set_b (set of int): atoms in current chain
    d_res_num (int):  min(set_a) - min(set_b)

  Returns:
    sel_a,sel_b (flex.bool): Updated selection
  """
  exist_in_both = set_a.intersection(set_b)
  for i in exist_in_both:
    if (sel_a[i] == False) and (sel_b[i - d_res_num] == False):
      sel_a[i] = True
      sel_b[i -  d_res_num] = True
  return sel_a, sel_b
