# written by Tom Ioerger (http://faculty.cs.tamu.edu/ioerger)
# send comments or suggestions to ioerger@cs.tamu.edu

from libtbx import adopt_init_args
import sys, types

# This implementation is based on the ideas in Gotoh (1982; JMB, 162:705-708),
#   which runs in quadratic time O(MN), i.e. product of sequence lengths.
# It does both global (Needleman/Wunsch) and local (Smith/Waterman) alignments.
# See usage examples at the end of this file.

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

class align(object):

  def __init__(self,
        seq_a,
        seq_b,
        style="global",
        gap_opening_penalty=1,
        gap_extension_penalty=1,
        similarity_function="identity"):
    assert style in ["global", "local", "no_end_gaps"]
    if (   similarity_function is None
        or similarity_function == "identity"):
      similarity_function = identity
    elif (similarity_function == "dayhoff"):
      similarity_function = dayhoff
    elif (similarity_function == "blosum50"):
      similarity_function = blosum50
    elif (isinstance(similarity_function, str)):
      raise RuntimeError(
        'Unknown similarity_function: "%s"' % similarity_function)
    adopt_init_args(self, locals())
    A,B = seq_a, seq_b
    m,n = self.m,self.n = len(A),len(B)

    # Mij is score of align of A[1..i] with B[1..j] ending in a match
    # Dij is score of align of A[1..i] with B[1..j] ending in a deletion (Ai,gap)
    # Iij is score of align of A[1..i] with B[1..j] ending in an insertion (gap,Bj)
    # E is direction matrix: -1=del, 0=match, +1=ins
    M = self.make_matrix(m+1,n+1)
    D = self.make_matrix(m+1,n+1)
    I = self.make_matrix(m+1,n+1)
    E = self.make_matrix(m+1,n+1)

    # initialize the matrices
    if style=="global":
      for i in range(1,m+1): M[i][0] = -self.gap_cost(i)
      for i in range(1,n+1): M[0][i] = -self.gap_cost(i)
      M[0][0] = 0
    # else (LOCAL, NO_END_GAPS) whole matrix initialized to 0 by default

    # fill in the matrices using dynamic programming
    for i in range(1,m+1):
      for j in range(1,n+1):

        if i==1: D[i][j] = M[i-1][j]-self.gap_cost(1)
        else: D[i][j] = max(M[i-1][j]-self.gap_cost(1),D[i-1][j]-gap_extension_penalty)

        if j==1: I[i][j] = M[i][j-1]-self.gap_cost(1)
        else: I[i][j] = max(M[i][j-1]-self.gap_cost(1),I[i][j-1]-gap_extension_penalty)

        M[i][j] = max(M[i-1][j-1]+similarity_function(A[i-1],B[j-1]),D[i][j],I[i][j])
        if style=="local": M[i][j] = max(M[i][j],0)

        if M[i][j]==D[i][j]: E[i][j] = 1    # deletion, i.e. of A[i]
        elif M[i][j]==I[i][j]: E[i][j] = -1 # insertion, i.e. of B[j]
        else: E[i][j] = 0

    (self.M,self.D,self.I,self.E) = (M,D,I,E)

  def gap_cost(self, width):
    return self.gap_opening_penalty + width * self.gap_extension_penalty

  def make_matrix(self,m,n):
    R = []
    for i in xrange(m):
      R.append([0]*n)
    return R

  def show_matrix(self, data, label=None, out=None):
    if (out is None): out = sys.stdout
    if (label is not None):
      print >> out, label
    for a in '  '+self.seq_b: print >> out, "%5c" % a,
    print >> out
    seq = ' '+self.seq_a
    for i,row in enumerate(data):
      print >> out, "%5c" % seq[i],
      for x in row: print >> out, "%5.1f" % x,
      print >> out

  def show_matrices(self, out=None):
    for label, data in [("D", self.D),
                        ("I", self.I),
                        ("M", self.M),
                        ("E", self.E)]:
      self.show_matrix(data=data, label=label)
    return self

  def score(self):
    (i,j) = self.endpt()
    return self.M[i][j]

  def endpt(self):
    if self.style=="global":
      return (self.m,self.n)

    elif self.style=="local":
      (best,ii,jj) = (self.M[0][0],0,0)
      for i in xrange(self.m+1):
        for j in xrange(self.n+1):
          if self.M[i][j]>best: (best,ii,jj) = (self.M[i][j],i,j)
      return (ii,jj)

    else: # NO_END_GAPS, search edges of matrix
      (best,ii,jj) = (self.M[0][0],0,0)
      for i in xrange(self.m+1):
        j = self.n
        if self.M[i][j]>best: (best,ii,jj) = (self.M[i][j],i,j)
      for j in xrange(self.n+1):
        i = self.m
        if self.M[i][j]>best: (best,ii,jj) = (self.M[i][j],i,j)
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
        if E[i][j]==-1: mcap('i'); j -= 1
        elif E[i][j]==1: mcap('d'); i -= 1
        elif E[i][j]==0: mcap('m'); i -= 1; j -= 1
      while i>0: mcap('d'); i -= 1
      while j>0: mcap('i'); j -= 1
      F,G = range(len(self.seq_a)), range(len(self.seq_b))
    else:
      (p,q) = (i,j)
      while M[i][j]>0:
        if E[i][j]==-1: mcap('i'); j -= 1
        elif E[i][j]==1: mcap('d'); i -= 1
        elif E[i][j]==0: mcap('m'); i -= 1; j -= 1
      F,G = range(i,p+1),range(j,q+1) # sub-sequences
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

    return alignment(
      similarity_function=self.similarity_function,
      a="".join(sa), b="".join(sb),
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

  def pretty_print(self,
        matches=None,
        out=None,
        block_size=20,
        n_block=1,
        top_name="reference",
        bottom_name="query",
        comment = None):
    if (matches is None): matches = self.matches()
    if (out is None): out = sys.stdout

    top_str = (top_name+" "*8)[0:8]
    bot_str = (bottom_name+" "*8)[0:8]
    ruler = ""
    count=0
    for ii in xrange(n_block):
      for jj in xrange(block_size):
        count += 1
        ruler += "%s"%( count%10 )
      ruler+="     "
    print >> out
    print >> out
    if comment is not None:
      print >> out, comment
    print >> out, "              "+ruler
    print >> out

    done=False
    n=len(self.a)
    count=0
    while not done:
      # top
      offset=count*block_size*n_block

      # top
      print >> out, top_str+"     ",
      for ii in xrange(n_block):
        start=offset+ii*block_size
        stop=offset+(ii+1)*block_size
        if stop > n:
          stop = n
        if start < n:
          tmp=self.a[start:stop]
          print >> out, tmp, "   ",
      print >> out

      #middle
      print >> out, "             ",
      for ii in xrange(n_block):
        start=offset+ii*block_size
        stop=offset+(ii+1)*block_size
        if stop > n:
          stop = n
        if start < n:
          tmp=matches[start:stop]
          print >> out, tmp, "   ",
      count += 1
      print >> out

      # bottom
      print >> out, bot_str+"     ",
      for ii in xrange(n_block):
        start=offset+ii*block_size
        stop=offset+(ii+1)*block_size
        if stop > n:
          stop = n
        if start < n:
          tmp=self.b[start:stop]
          print >> out, tmp, "   ",
      print >> out
      print >> out
      if count*block_size*n_block>n:
        done=True

    return out

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

def exercise_similarity_scores():
  from scitbx.array_family import flex
  for m in [dayhoff_mdm78_similarity_scores, blosum50_similarity_scores]:
    # exception if not symmetric
    flex.double(m).matrix_symmetric_as_packed_l(relative_epsilon=0)

def exercise():
  A = "AAAGGTT"
  B = "AAATT"
  obj = align(A,B)
  obj.show_matrices()

  print "score=%.1f" % obj.score()
  alignment = obj.extract_alignment()
  print alignment.match_codes
  print alignment.a
  print alignment.identity_matches()
  print alignment.b

  # 1rra vs. 1bli
  A = "AESSADKFKRQHMDTEGPSKSSPTYCNQMMKRQGMTKGSCKPVNTFVHEPLEDVQAICSQGQVTCKNGRNNCHKSSSTLRITDCRLKGSSKYPNCDYTTTDSQKHIIIACDGNPYVPVHFDASV"
  B = "DNSRYTHFLTQHYDAKPQGRDDRYCESIMRRRGLTSPCKDINTFIHGNKRSIKAICENKNGNPHRENLRISKSSFQVTTCKLHGGSPWPPCQYRATAGFRNVVVACENGLPVHLDQSIFRRP"
  obj = align(A,B,gap_opening_penalty=150,gap_extension_penalty=20,similarity_function=dayhoff,style="global")

  print "\n1rra vs. 1bli; GLOBAL allignment; mdm78"
  print "score=%.1f" % obj.score()
  alignment = obj.extract_alignment()

  print alignment.match_codes
  print alignment.a
  print alignment.dayhoff_matches()
  print alignment.b


  # 1rra vs. 1bli
  A = "AESSADKFKRQHMDTEGPSKSSPTYCNQMMKRQGMTKGSCKPVNTFVHEPLEDVQAICSQGQVTCKNGRNNCHKSSSTLRITDCRLKGSSKYPNCDYTTTDSQKHIIIACDGNPYVPVHFDASV"
  B = "DNSRYTHFLTQHYDAKPQGRDDRYCESIMRRRGLTSPCKDINTFIHGNKRSIKAICENKNGNPHRENLRISKSSFQVTTCKLHGGSPWPPCQYRATAGFRNVVVACENGLPVHLDQSIFRRP"
  obj = align(A,B,gap_opening_penalty=150,gap_extension_penalty=20,similarity_function="dayhoff",style="local")

  print "\n1rra vs. 1bli; LOCAL allignment; mdm78"
  print "score=%.1f" % obj.score()
  alignment = obj.extract_alignment()

  print alignment.match_codes
  print alignment.a
  print alignment.dayhoff_matches()
  print alignment.b



  # 1rra vs. 1bli
  A = "AESSADKFKRQHMDTEGPSKSSPTYCNQMMKRQGMTKGSCKPVNTFVHEPLEDVQAICSQGQVTCKNGRNNCHKSSSTLRITDCRLKGSSKYPNCDYTTTDSQKHIIIACDGNPYVPVHFDASV"
  B = "DNSRYTHFLTQHYDAKPQGRDDRYCESIMRRRGLTSPCKDINTFIHGNKRSIKAICENKNGNPHRENLRISKSSFQVTTCKLHGGSPWPPCQYRATAGFRNVVVACENGLPVHLDQSIFRRP"
  obj = align(A,B,gap_opening_penalty=10,gap_extension_penalty=2,similarity_function=blosum50,style="global")

  print "\n1rra vs. 1bli; GLOBAL allignment; blosum50"
  print "score=%.1f" % obj.score()
  alignment = obj.extract_alignment()

  print alignment.match_codes
  print alignment.a
  print alignment.matches()
  print alignment.b

  # 1rra vs. 1bli
  A = "AESSADKFKRQHMDTEGPSKSSPTYCNQMMKRQGMTKGSCKPVNTFVHEPLEDVQAICSQGQVTCKNGRNNCHKSSSTLRITDCRLKGSSKYPNCDYTTTDSQKHIIIACDGNPYVPVHFDASV"
  B = "DNSRYTHFLTQHYDAKPQGRDDRYCESIMRRRGLTSPCKDINTFIHGNKRSIKAICENKNGNPHRENLRISKSSFQVTTCKLHGGSPWPPCQYRATAGFRNVVVACENGLPVHLDQSIFRRP"
  obj = align(A,B,gap_opening_penalty=10,gap_extension_penalty=2,similarity_function="blosum50",style="local")

  print "\n1rra vs. 1bli; LOCAL allignment; blosum50"
  print "score=%.1f" % obj.score()
  alignment = obj.extract_alignment()

  print alignment.match_codes
  print alignment.a
  print alignment.matches(similarity_function=blosum50, is_similar_threshold=0)
  print alignment.b
  print
  alignment.pretty_print(
    matches = None,
    out = None,
    block_size = 50,
    n_block = 1,
    top_name = "1rra",
    bottom_name = "1bli",
    comment = """pretty_print is pretty pretty""")

  print "OK" # necessary for auto_build checking

if __name__=="__main__":
  exercise_similarity_scores()
  exercise()
