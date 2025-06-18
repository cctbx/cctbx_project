"""algorithms for aligning two protein
sequences, where each sequence is represented as a string of one-letter
amino-acid codes"""

from __future__ import absolute_import, division, print_function
import scitbx.array_family.flex
from scitbx.array_family import flex

import boost_adaptbx.boost.python as bp
from six import string_types
from six.moves import zip
from six.moves import range
ext = bp.import_ext("mmtbx_alignment_ext")

"""
Written by Tom Ioerger (http://faculty.cs.tamu.edu/ioerger).
Originally in Python, Oleg Sobolev moved to C++ for better performance.
Send comments or suggestions to: ioerger@cs.tamu.edu

This implementation provides algorithms for aligning two protein
sequences, where each sequence is represented as a string of one-letter
amino-acid codes. The implementation is based on the ideas of Gotoh
(1982) and runs in quadratic time O(M*N), where M and N are the
sequence lengths. It does both global (Needleman & Wunsch, 1970) and
local (Smith & Waterman, 1981) alignments, assuming affine (linear) gap
penalties (for which default gap-cost parameters may be changed by the
user). Alignments are based on maximizing similarity. Similarity scores
between amino acids are specified via symmetric matrices. Similarity
matrices of Dayhoff (1978) and BLOSUM50 (Henikoff & Henikoff (1992),
http://en.wikipedia.org/wiki/BLOSUM) are provided. User-supplied
matrices are also supported (this feature also enables alignment of
non-amino-acid sequences).

Dayhoff, M.O. (1978).
Atlas of Protein Sequence and Structure, Vol. 5 suppl. 3, 345-352.

Gotoh, O. (1982). J. Mol. Biol. 162, 705-708.

Henikoff & Henikoff (1992). PNAS 89, 10915-10919

Needleman, S. & Wunsch, C. (1970). J. Mol. Biol. 48(3), 443-53.

Smith, T.F. & Waterman M.S. (1981). J. Mol. Biol. 147, 195-197.
"""

from libtbx import adopt_init_args
import sys

# based on maximizing similarity (rather than minimizing distance)
# give gap weights as positive; they will be subtracted

# similarity matrices:
#   identity by default
#   Dayhoff (1978) matrix (like pam250)
#   blosum50
#   or arbitrary similarity function specified by user
# See the definitions of dayhoff(), blosum50() for example
# similarity functions.

# matrices go from (0,0) to (M,N) inclusive; M,N are lengths of A,B
# (1,1) represents (A[0],B[0])
# edge entries (i,0) and (0,j) represent initial gaps

# default gap weights (appropriate for identity similarity):
#   gap_opening_penalty = 1
#   gap_extension_penalty = 1
# suggestion for dayhoff:
#   gap_opening_penalty = 150
#   gap_extension_penalty = 20
#   log-likelihoods, scaled up by a factor of 10, with mean=0
#   so gaps cost approx 1.5+0.22*len per "match"

class align(ext.align):

  def __init__(self,
        seq_a,
        seq_b,
        style="global",
        gap_opening_penalty=1,
        gap_extension_penalty=1,
        similarity_function="identity",
        masking_a=None):
    assert style in ["global", "local", "no_end_gaps"]
    sim_fun_str = similarity_function
    self.similarity_function_call = None
    if (   similarity_function is None
        or similarity_function == "identity"):
      self.similarity_function_call = identity
    elif (similarity_function == "dayhoff"):
      self.similarity_function_call = dayhoff
    elif (similarity_function == "blosum50"):
      self.similarity_function_call = blosum50
    elif (isinstance(similarity_function, str)):
      raise RuntimeError(
        'Unknown similarity_function: "%s"' % similarity_function)
    from six import string_types
    if isinstance(seq_a, string_types):
      seq_a = seq_a.upper()
      seq_b = seq_b.upper()
    else:
      seq_a = flex.std_string(seq_a)
      seq_b = flex.std_string(seq_b)

    if masking_a:  #  Gap penalty for insertion scaled by masking_a after
                   #   corresponding position in sequence a and scaled by
                   #   masking_a for deletion before corresponding position in
                   #   sequence a.  Allows increasing mask everywhere
                   #   except where you want a gap to occur
      assert len(list(masking_a))==len(seq_a)
    else:
      masking_a=[1] * len(seq_a)  #  no masking; standard gap penalty everywhere
    m = flex.float(masking_a)
    super(align, self).__init__(
        seq_a=seq_a,
        seq_b=seq_b,
        masking=m,
        style=style,
        gap_opening_penalty=gap_opening_penalty,
        gap_extension_penalty=gap_extension_penalty,
        similarity_function=sim_fun_str,
        )
    self.seq_a = seq_a
    self.seq_b = seq_b
    self.style = style
    self.similarity_function = sim_fun_str
    m,n = self.m,self.n = len(seq_a),len(seq_b)

  def show_matrix(self, data, label=None, out=None):
    if (out is None): out = sys.stdout
    if (label is not None):
      print(label, file=out)
    for a in '  '+self.seq_b: print("%5c" % a, end=' ', file=out)
    print(file=out)
    seq = ' '+self.seq_a
    for i in range(self.m):
      print("%5c" % seq[i], end=' ', file=out)
      for j in range(self.n):
        print("%5.1f" % data[i,j], end=' ', file=out)
      print(file=out)

  def show_matrices(self, out=None):
    for label, data in [("D", self.D),
                        ("I", self.I),
                        ("M", self.M),
                        ("E", self.E)]:
      self.show_matrix(data=data, label=label)
    return self

  def score(self):
    (i,j) = self.endpt()
    return self.M[i,j]

  def endpt(self):
    if self.style=="global":
      return (self.m,self.n)

    elif self.style=="local":
      (best,ii,jj) = (self.M[0,0],0,0)
      for i in range(self.m+1):
        for j in range(self.n+1):
          if self.M[i,j]>best: (best,ii,jj) = (self.M[i,j],i,j)
      return (ii,jj)

    else: # NO_END_GAPS, search edges of matrix
      (best,ii,jj) = (self.M[0,0],0,0)
      for i in range(self.m+1):
        j = self.n
        if self.M[i,j]>best: (best,ii,jj) = (self.M[i,j],i,j)
      for j in range(self.n+1):
        i = self.m
        if self.M[i,j]>best: (best,ii,jj) = (self.M[i,j],i,j)
      return (ii,jj)

  def extract_alignment(self):
    # if GLOBAL, match_codes must go from corner to corner
    #   (though there may be end gaps - still have to search edges)
    # if LOCAL, search for best-scoring end pair anywhere in the matrix,
    #   trace back to entry where score=0 (also use 0 in max)
    # if NO_END_GAPS: like global but initialize rows to 0,
    #   truncate terminal gaps, do NOT take max of each matrix entry with 0
    (M,D,I,E) = (self.M,self.D,self.I,self.E)
    (m,n) = (self.m,self.n)

    match_codes = []
    mcap = match_codes.append
    (i,j) = self.endpt()
    if self.style=="global":
      while i>0 and j>0:
        if E[i,j]==-1: mcap('i'); j -= 1
        elif E[i,j]==1: mcap('d'); i -= 1
        elif E[i,j]==0: mcap('m'); i -= 1; j -= 1
      while i>0: mcap('d'); i -= 1
      while j>0: mcap('i'); j -= 1
      F,G = list(range(len(self.seq_a))), list(range(len(self.seq_b)))
    else:
      (p,q) = (i,j)
      while M[i,j]>0:
        if E[i,j]==-1: mcap('i'); j -= 1
        elif E[i,j]==1: mcap('d'); i -= 1
        elif E[i,j]==0: mcap('m'); i -= 1; j -= 1
      F,G = list(range(i,p+1)),list(range(j,q+1)) # sub-sequences
    match_codes.reverse()
    match_codes = "".join(match_codes)

    sa,sb = [], []
    ia,ib = [], []
    u, v = 0, 0
    for a in match_codes:
      if a=='d':
        i = F[u]
        u += 1
        ia.append(i)
        sa.append(self.seq_a[i])
        ib.append(None)
        sb.append('-')
      elif a=='i':
        ia.append(None)
        sa.append('-')
        i = G[v]
        v += 1
        ib.append(i)
        sb.append(self.seq_b[i])
      elif a=='m':
        i = F[u]
        u += 1
        ia.append(i)
        sa.append(self.seq_a[i])
        i = G[v]
        v += 1
        ib.append(i)
        sb.append(self.seq_b[i])
    from six import string_types
    if isinstance(self.seq_a, string_types):
      sa = "".join(sa)
      sb = "".join(sb)
    return alignment(
      similarity_function=self.similarity_function_call,
      a=sa, b=sb,
      i_seqs_a=ia, i_seqs_b=ib,
      match_codes=match_codes)

class alignment(object):

  def __init__(self,
        similarity_function,
        a, b,
        i_seqs_a, i_seqs_b,
        match_codes):
    adopt_init_args(self, locals())

  def matches(self, similarity_function=None, is_similar_threshold=0):
    if (similarity_function is None):
      similarity_function = self.similarity_function
    result = []
    ap = result.append
    for a,b in zip(self.a, self.b):
      if (a == b):
        ap("|")
      elif (similarity_function(a, b) > is_similar_threshold):
        ap("*")
      else:
        ap(" ")
    return "".join(result)

  def identity_matches(self):
    return self.matches(similarity_function=identity)

  def dayhoff_matches(self, is_similar_threshold=0):
    return self.matches(
      similarity_function=dayhoff,
      is_similar_threshold=is_similar_threshold)

  def blosum50_matches(self, is_similar_threshold=0):
    return self.matches(
      similarity_function=blosum50,
      is_similar_threshold=is_similar_threshold)

  def calculate_sequence_identity(self, skip_chars=()):
    """
    Returns fractional sequence identity, defined here as the number of matches
    between the aligned sequences divided by the number of valid residues in
    the first sequence.  The optional argument skip_chars may be used to
    pass over null residues (e.g. 'X' for a gap in a PDB chain).
    """
    skip_chars = list(skip_chars)
    skip_chars.append("-")
    n_matches = n_total = 0
    for a, b in zip(self.a, self.b):
      # XXX should gaps in 'b' be discounted?
      if (not a in skip_chars) : # and (not b in skip_cars)
        n_total += 1
        if (a == b):
          n_matches += 1
    if (n_total == 0) or (n_matches == 0):
      return 0.
    return n_matches / n_total

  def pretty_print(self,
        matches=None,
        out=None,
        block_size=20,
        n_block=1,
        top_name="reference",
        bottom_name="query",
        comment = None,
        show_ruler=True):
    if (matches is None): matches = self.matches()
    if (out is None): out = sys.stdout

    top_str = (top_name+" "*8)[0:8]
    bot_str = (bottom_name+" "*8)[0:8]
    ruler = ""
    count=0
    for ii in range(n_block):
      for jj in range(block_size):
        count += 1
        ruler += "%s"%( count%10 )
      ruler+="     "
    print(file=out)
    print(file=out)
    if comment is not None:
      print(comment, file=out)
    if (show_ruler):
      print("              "+ruler, file=out)
      print(file=out)

    done=False
    n=len(self.a)
    count=0
    while not done:
      # top
      offset=count*block_size*n_block

      # top
      print(top_str+"     ", end=' ', file=out)
      for ii in range(n_block):
        start=offset+ii*block_size
        stop=offset+(ii+1)*block_size
        if stop > n:
          stop = n
        if start < n:
          tmp=self.a[start:stop]
          print(tmp, "   ", end=' ', file=out)
      print(file=out)

      #middle
      print("             ", end=' ', file=out)
      for ii in range(n_block):
        start=offset+ii*block_size
        stop=offset+(ii+1)*block_size
        if stop > n:
          stop = n
        if start < n:
          tmp=matches[start:stop]
          print(tmp, "   ", end=' ', file=out)
      count += 1
      print(file=out)

      # bottom
      print(bot_str+"     ", end=' ', file=out)
      for ii in range(n_block):
        start=offset+ii*block_size
        stop=offset+(ii+1)*block_size
        if stop > n:
          stop = n
        if start < n:
          tmp=self.b[start:stop]
          print(tmp, "   ", end=' ', file=out)
      print(file=out)
      print(file=out)
      if count*block_size*n_block>n:
        done=True

    return out

  def exact_match_selections(self):
    i_seqs = flex.size_t()
    j_seqs = flex.size_t()
    for a, b, i, j in zip(self.a, self.b, self.i_seqs_a, self.i_seqs_b):
      if(a == b and a not in ["-", None]):
        if(i is not None and j is not None):
          i_seqs.append(i)
          j_seqs.append(j)
    return i_seqs, j_seqs

  def exact_mismatch_selection(self):
    """Returns i_seqs and j_seqs of residues that do not match, excluding
    insertions and deletions.

    Returns:
        (flex.size_t, flex.size_t): i_seqs that mismatch
    """
    i_seqs = flex.size_t()
    j_seqs = flex.size_t()
    for a, b, i, j in zip(self.a, self.b, self.i_seqs_a, self.i_seqs_b):
      if(a != b and a not in ["-", None] and b not in ["-", None]):
        if(i is not None and j is not None):
          i_seqs.append(i)
          j_seqs.append(j)
    return i_seqs, j_seqs

# amino acid similarity scores from Dayhoff's 1978 paper; like PAM250?
dayhoff_mdm78_similarity_scores = [
  [ 18,-20,  3,  3,-35, 13,-14, -5,-12,-19,-11,  2, 11, -4,-15, 11, 12,  2,-58,-35],
  [-20,119,-51,-53,-43,-34,-34,-23,-54,-60,-52,-36,-28,-54,-36,  0,-22,-19,-78,  3],
  [  3,-51, 39, 34,-56,  6,  7,-24,  1,-40,-26, 21,-10, 16,-13,  3, -1,-21,-68,-43],
  [  3,-53, 34, 38,-54,  2,  7,-20, -1,-34,-21, 14, -6, 25,-11,  0, -4,-18,-70,-43],
  [-35,-43,-56,-54, 91,-48,-18, 10,-53, 18,  2,-35,-46,-47,-45,-32,-31,-12,  4, 70],
  [ 13,-34,  6,  2,-48, 48,-21,-26,-17,-41,-28,  3, -5,-12,-26, 11,  0,-14,-70,-52],
  [-14,-34,  7,  7,-18,-21, 65,-24,  0,-21,-21, 16, -2, 29, 16, -8,-13,-22,-28, -1],
  [ -5,-23,-24,-20, 10,-26,-24, 45,-19, 24, 22,-18,-20,-20,-20,-14,  1, 37,-51, -9],
  [-12,-54,  1, -1,-53,-17,  0,-19, 47,-29,  4, 10,-11,  7, 34, -2,  0,-24,-35,-44],
  [-19,-60,-40,-34, 18,-41,-21, 24,-29, 59, 37,-29,-25,-18,-30,-28,-17, 19,-18, -9],
  [-11,-52,-26,-21,  2,-28,-21, 22,  4, 37, 64,-17,-21,-10, -4,-16, -6, 18,-42,-24],
  [  2,-36, 21, 14,-35,  3, 16,-18, 10,-29,-17, 20, -5,  8,  0,  7,  4,-17,-42,-21],
  [ 11,-28,-10, -6,-46, -5, -2,-20,-11,-25,-21, -5, 59,  2, -2,  9,  3,-12,-56,-49],
  [ -4,-54, 16, 25,-47,-12, 29,-20,  7,-18,-10,  8,  2, 40, 13, -5, -8,-19,-48,-40],
  [-15,-36,-13,-11,-45,-26, 16,-20, 34,-30, -4,  0, -2, 13, 61, -3, -9,-25, 22,-42],
  [ 11,  0,  3,  0,-32, 11, -8,-14, -2,-28,-16,  7,  9, -5, -3, 16, 13,-10,-25,-28],
  [ 12,-22, -1, -4,-31,  0,-13,  1,  0,-17, -6,  4,  3, -8, -9, 13, 26,  3,-52,-27],
  [  2,-19,-21,-18,-12,-14,-22, 37,-24, 19, 18,-17,-12,-19,-25,-10,  3, 43,-62,-25],
  [-58,-78,-68,-70,  4,-70,-28,-51,-35,-18,-42,-42,-56,-48, 22,-25,-52,-62,173, -2],
  [-35,  3,-43,-43, 70,-52, -1, -9,-44, -9,-24,-21,-49,-40,-42,-28,-27,-25, -2,101]]

blosum50_similarity_scores = [
  [ 5,-2,-1,-2,-1,-3, 0,-2,-1,-1,-2,-1,-1,-1,-1,-2, 1, 0, 0,-3,-1,-2,-1],
  [-2, 5,-3, 5, 1,-4,-1, 0,-4, 0,-4,-3, 4,-2, 0,-1, 0, 0,-4,-5,-1,-3, 2],
  [-1,-3,13,-4,-3,-2,-3,-3,-2,-3,-2,-2,-2,-4,-3,-4,-1,-1,-1,-5,-2,-3,-3],
  [-2, 5,-4, 8, 2,-5,-1,-1,-4,-1,-4,-4, 2,-1, 0,-2, 0,-1,-4,-5,-1,-3, 1],
  [-1, 1,-3, 2, 6,-3,-3, 0,-4, 1,-3,-2, 0,-1, 2, 0,-1,-1,-3,-3,-1,-2, 5],
  [-3,-4,-2,-5,-3, 8,-4,-1, 0,-4, 1, 0,-4,-4,-4,-3,-3,-2,-1, 1,-2, 4,-4],
  [ 0,-1,-3,-1,-3,-4, 8,-2,-4,-2,-4,-3, 0,-2,-2,-3, 0,-2,-4,-3,-2,-3,-2],
  [-2, 0,-3,-1, 0,-1,-2,10,-4, 0,-3,-1, 1,-2, 1, 0,-1,-2,-4,-3,-1, 2, 0],
  [-1,-4,-2,-4,-4, 0,-4,-4, 5,-3, 2, 2,-3,-3,-3,-4,-3,-1, 4,-3,-1,-1,-3],
  [-1, 0,-3,-1, 1,-4,-2, 0,-3, 6,-3,-2, 0,-1, 2, 3, 0,-1,-3,-3,-1,-2, 1],
  [-2,-4,-2,-4,-3, 1,-4,-3, 2,-3, 5, 3,-4,-4,-2,-3,-3,-1, 1,-2,-1,-1,-3],
  [-1,-3,-2,-4,-2, 0,-3,-1, 2,-2, 3, 7,-2,-3, 0,-2,-2,-1, 1,-1,-1, 0,-1],
  [-1, 4,-2, 2, 0,-4, 0, 1,-3, 0,-4,-2, 7,-2, 0,-1, 1, 0,-3,-4,-1,-2, 0],
  [-1,-2,-4,-1,-1,-4,-2,-2,-3,-1,-4,-3,-2,10,-1,-3,-1,-1,-3,-4,-2,-3,-1],
  [-1, 0,-3, 0, 2,-4,-2, 1,-3, 2,-2, 0, 0,-1, 7, 1, 0,-1,-3,-1,-1,-1, 4],
  [-2,-1,-4,-2, 0,-3,-3, 0,-4, 3,-3,-2,-1,-3, 1, 7,-1,-1,-3,-3,-1,-1, 0],
  [ 1, 0,-1, 0,-1,-3, 0,-1,-3, 0,-3,-2, 1,-1, 0,-1, 5, 2,-2,-4,-1,-2, 0],
  [ 0, 0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1, 0,-1,-1,-1, 2, 5, 0,-3, 0,-2,-1],
  [ 0,-4,-1,-4,-3,-1,-4,-4, 4,-3, 1, 1,-3,-3,-3,-3,-2, 0, 5,-3,-1,-1,-3],
  [-3,-5,-5,-5,-3, 1,-3,-3,-3,-3,-2,-1,-4,-4,-1,-3,-4,-3,-3,15,-3, 2,-2],
  [-1,-1,-2,-1,-1,-2,-2,-1,-1,-1,-1,-1,-1,-2,-1,-1,-1, 0,-1,-3,-1,-1,-1],
  [-2,-3,-3,-3,-2, 4,-3, 2,-1,-2,-1, 0,-2,-3,-1,-1,-2,-2,-1, 2,-1, 8,-2],
  [-1, 2,-3, 1, 5,-4,-2, 0,-3, 1,-3,-1, 0,-1, 4, 0, 0,-1,-3,-2,-1,-2, 5]]

amino_acid_codes = [
  "A",
  "C",
  "D",
  "E",
  "F",
  "G",
  "H",
  "I",
  "K",
  "L",
  "M",
  "N",
  "P",
  "Q",
  "R",
  "S",
  "T",
  "V",
  "W",
  "Y",
  ]

blosum62_similarity_scores = [
  [  4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2 ],
  [  0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2 ],
  [ -2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3 ],
  [ -1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2 ],
  [ -2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3 ],
  [  0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3 ],
  [ -2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2 ],
  [ -1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1 ],
  [ -1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2 ],
  [ -1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1 ],
  [ -1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1 ],
  [ -2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2 ],
  [ -1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3 ],
  [ -1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1 ],
  [ -1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2 ],
  [  1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2 ],
  [  0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2 ],
  [  0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1 ],
  [ -3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2 ],
  [ -2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7 ],
  ]

def identity(a, b):
  return int(a == b)

def dayhoff(a,b):
  AAs = "ACDEFGHIKLMNPQRSTVWY"
  (i,j) = (AAs.find(a),AAs.find(b))
  if i==-1 or j==-1: return 0 # should be mean value?
  return dayhoff_mdm78_similarity_scores[i][j]

def blosum50(a,b):
  AAs = "ABCDEFGHIKLMNPQRSTVWXYZ"
  (i,j) = (AAs.find(a),AAs.find(b))
  if i==-1 or j==-1: return 0 # should be mean value?
  return blosum50_similarity_scores[i][j]

def blosum62(left, right):
  try:
    index_left = amino_acid_codes.index( left )
    index_right = amino_acid_codes.index( right )
  except ValueError:
    return 0
  return blosum62_similarity_scores[index_left][index_right]

def exercise_similarity_scores():
  from scitbx.array_family import flex
  for m in [dayhoff_mdm78_similarity_scores, blosum50_similarity_scores]:
    assert flex.double(m).matrix_is_symmetric(relative_epsilon=1e-15)

class pairwise_global(ext.pairwise_global):
  def __init__(self, seq1, seq2):
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    ext.pairwise_global.__init__(self, seq1, seq2)

class pairwise_global_wrapper(pairwise_global):

  def range_matches_from_aligned_sequences(self):
    a1 = self.result1
    a2 = self.result2

    #in one pass, algorithmically determine the range matches from the given alignment
    assert len(a1) == len(a2)
    #use three pointers, ptr0:overall alignment; ptr1:sequence 1; ptr2:sequence 2
    ptr0 = -1; ptr1 = -1; ptr2 = -1

    in_an_aligned_range = False
    overall_ranges1 = []
    overall_ranges2 = []
    while 1:
      ptr0+=1
      if ptr0 == len(a1): break #done with whole alignment
      if a1[ptr0]!='-':  ptr1+=1
      if a2[ptr0]!='-':  ptr2+=1
      if a1[ptr0]!='-' and a2[ptr0]!='-':
        # we are inside a matching range
        if not in_an_aligned_range:
          overall_ranges1.append([ptr1,ptr1+1])
          overall_ranges2.append([ptr2,ptr2+1])
          in_an_aligned_range = True
        else:
          overall_ranges1[-1][1]+=1
          overall_ranges2[-1][1]+=1

      else:
        in_an_aligned_range = False

    return overall_ranges1,overall_ranges2

  def calculate_sequence_identity(self, skip_chars=()):
    a1 = self.result1
    a2 = self.result2
    assert len(a1) == len(a2)
    n_aligned_residues = 0
    n_matching = 0
    skip_chars = list(skip_chars)
    skip_chars.append("-")
    for i in range(len(a1)):
      if (not a1[i] in skip_chars) and (not a2[i] in skip_chars):
        n_aligned_residues += 1
        if a1[i] == a2[i] :
          n_matching += 1
    if n_matching == 0 or n_aligned_residues == 0 :
      return 0
    return float(n_matching) / float(n_aligned_residues)
