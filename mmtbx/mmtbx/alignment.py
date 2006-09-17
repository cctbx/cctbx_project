# written by Tom Ioerger (http://faculty.cs.tamu.edu/ioerger)
# send comments or suggestions to ioerger@cs.tamu.edu

import sys,types

# This implementation is based on the ideas in Gotoh (1982; JMB, 162:705-708),
#   which runs in quadratic time O(MN), i.e. product of sequence lengths.
# It does both global (Needleman/Wunsch) and local (Smith/Waterman) alignments.
# See usage examples at the end of this file.

# based on maximizing similarity (rather than minimizing distance)
# give gap weights as positive; they will be subtracted

# options/parameters:
#   gap weights:
#     affine gap costs: cost = gop + gep*length
#     gop = gap open penalty
#     gep = gap extension penalty
#   similarity matrices:
#     identity by default, Dayhoff (1978) matrix (like pam250) provided
#     or arbitrary similarity function specified by user

# matrices go from (0,0) to (M,N) inclusive; M,N are lengths of A,B
# (1,1) represents (A[0],B[0])
# edge entries (i,0) and (0,j) represent initial gaps

# sim parameter can be:
#   a pair of (match score,mismatch score), e.g. sim=(1,-0.5)
#   a function that takes a pair and returns score
#   default is identity: sim=(1,0)

# default gap weights: gop=1, gep=1 (appropriate for identity sim)
# suggestion for dayhoff: gop=150, gep=20
#   log-likelihoods, scaled up by a factor of 10, with mean=0
#   so gaps cost approx 1.5+0.22*len per "match"

class align:

  # for specifying alignment styles, e.g. align(...,style=align.LOCAL)
  LOCAL       = 1 # Smith/Waterman
  GLOBAL      = 2 # Needleman/Wunsch (default)
  NO_END_GAPS = 3 # like NW, except end gaps don't count

  def __init__(self,A,B,gop=1,gep=1,sim=(1,0),style=GLOBAL,verbose=False):
    (self.sim,self.style) = (sim,style)
    (self.A,self.B) = (A,B)
    (m,n) = (self.m,self.n) = (len(A),len(B))
    (self.gop,self.gep) = (gop,gep)

    # Mij is score of align of A[1..i] with B[1..j] ending in a match
    # Dij is score of align of A[1..i] with B[1..j] ending in a deletion (Ai,gap)
    # Iij is score of align of A[1..i] with B[1..j] ending in an insertion (gap,Bj)
    # E is direction matrix: -1=del, 0=match, +1=ins
    M = self.make_array(m+1,n+1)
    D = self.make_array(m+1,n+1)
    I = self.make_array(m+1,n+1)
    E = self.make_array(m+1,n+1)

    # initialize the matrices
    if style==self.GLOBAL:
      for i in range(1,m+1): M[i][0] = -self.gap(i)
      for i in range(1,n+1): M[0][i] = -self.gap(i)
      M[0][0] = 0
    # else (LOCAL, NO_END_GAPS) whole matrix initialized to 0 by default

    # fill in the matrices using dynamic programming
    for i in range(1,m+1):
      for j in range(1,n+1):

        if i==1: D[i][j] = M[i-1][j]-self.gap(1)
        else: D[i][j] = max(M[i-1][j]-self.gap(1),D[i-1][j]-gep)

        if j==1: I[i][j] = M[i][j-1]-self.gap(1)
        else: I[i][j] = max(M[i][j-1]-self.gap(1),I[i][j-1]-gep)

        M[i][j] = max(M[i-1][j-1]+self.eval_sim(A[i-1],B[j-1],sim),D[i][j],I[i][j])
        if style==self.LOCAL: M[i][j] = max(M[i][j],0)

        if M[i][j]==D[i][j]: E[i][j] = 1    # deletion, i.e. of A[i]
        elif M[i][j]==I[i][j]: E[i][j] = -1 # insertion, i.e. of B[j]
        else: E[i][j] = 0

    (self.M,self.D,self.I,self.E) = (M,D,I,E)

    if verbose==True:
      print "D"; self.print_array(D)
      print "I"; self.print_array(I)
      print "M"; self.print_array(M)
      print "E"; self.print_array(E)

  def eval_sim(self,a,b,sim):
    if isinstance(sim,tuple):
      if a==b: return sim[0]
      return sim[1]
    elif isinstance(sim,types.FunctionType):
      return sim(a,b)

  def gap(self,i): return self.gop+i*self.gep

  # num digits?

  def print_array(self,R):
    for a in '  '+self.B: print "%5c" % a,
    print
    seq = ' '+self.A
    for i in xrange(len(R)):
      row = R[i]
      print "%5c" % seq[i],
      for x in row: print "%5.1f" % x,
      print

  def make_array(self,m,n):
    R = []
    for i in xrange(m):
      R.append([0]*n)
    return R

  def score(self):
    (i,j) = self.endpt()
    return self.M[i][j]

  def endpt(self):
    if self.style==self.GLOBAL:
      return (self.m,self.n)

    elif self.style==self.LOCAL:
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

  # if GLOBAL, traceback must go from corner to corner
  #   (though there may be end gaps - still have to search edges)
  # if LOCAL, search for best-scoring end pair anywhere in the matrix,
  #   trace back to entry where score=0 (also use 0 in max)
  # if NO_END_GAPS: like global but initialize rows to 0, truncate terminal gaps,
  #   do NOT take max of each matrix entry with 0

  # returns a pair of strings
  # if local, how does caller know coordinates? return (i,j) too?

  def extract_alignment(self):
    (M,D,I,E) = (self.M,self.D,self.I,self.E)
    (m,n) = (self.m,self.n)

    traceback = ""
    (i,j) = self.endpt()
    if self.style==self.GLOBAL:
      while i>0 and j>0:
        if E[i][j]==-1: traceback += 'i'; j -= 1
        elif E[i][j]==1: traceback += 'd'; i -= 1
        elif E[i][j]==0: traceback += 'm'; i -= 1; j -= 1
      while i>0: traceback += 'd'; i -= 1
      while j>0: traceback += 'i'; j -= 1
      (F,G) = (self.A,self.B)
    else:
      (p,q) = (i,j)
      while M[i][j]>0:
        if E[i][j]==-1: traceback += 'i'; j -= 1
        elif E[i][j]==1: traceback += 'd'; i -= 1
        elif E[i][j]==0: traceback += 'm'; i -= 1; j -= 1
      (F,G) = (self.A[i:p+1],self.B[j:q+1]) # sub-sequences
    traceback = self.reverse(traceback)
    print traceback

    (r,s) = ("","")
    (u,v) = (0,0)
    for a in traceback:
      if a=='d': r += F[u]; u += 1; s += '-'
      if a=='i': r += '-';          s += G[v]; v += 1
      if a=='m': r += F[u]; u += 1; s += G[v]; v += 1

    return (r,s)

  def reverse(self,s):
    foo = ""
    for i in xrange(len(s)): foo += s[-i-1]
    return foo

  # amino acid similarity scores from Dayhoff's 1978 paper; like PAM250?

def dayhoff(a,b):
  AAs = "ACDEFGHIKLMNPQRSTVWY"
  mdm78 = [ \
     [ 18,-20,  3,  3,-35, 13,-14, -5,-12,-19,-11,  2, 11, -4,-15, 11, 12,  2,-58,-35], \
     [-20,119,-51,-53,-43,-34,-34,-23,-54,-60,-52,-36,-28,-54,-36,  0,-22,-19,-78,  3], \
     [  3,-51, 39, 34,-56,  6,  7,-24,  1,-40,-26, 21,-10, 16,-13,  3, -1,-21,-68,-43], \
     [  3,-53, 34, 38,-54,  2,  7,-20, -1,-34,-21, 14, -6, 25,-11,  0, -4,-18,-70,-43], \
     [-35,-43,-56,-54, 91,-48,-18, 10,-53, 18,  2,-35,-46,-47,-45,-32,-31,-12,  4, 70], \
     [ 13,-34,  6,  2,-48, 48,-21,-26,-17,-41,-28,  3, -5,-12,-26, 11,  0,-14,-70,-52], \
     [-14,-34,  7,  7,-18,-21, 65,-24,  0,-21,-21, 16, -2, 29, 16, -8,-13,-22,-28, -1], \
     [ -5,-23,-24,-20, 10,-26,-24, 45,-19, 24, 22,-18,-20,-20,-20,-14,  1, 37,-51, -9], \
     [-12,-54,  1, -1,-53,-17,  0,-19, 47,-29,  4, 10,-11,  7, 34, -2,  0,-24,-35,-44], \
     [-19,-60,-40,-34, 18,-41,-21, 24,-29, 59, 37,-29,-25,-18,-30,-28,-17, 19,-18, -9], \
     [-11,-52,-26,-21,  2,-28,-21, 22,  4, 37, 64,-17,-21,-10, -4,-16, -6, 18,-42,-24], \
     [  2,-36, 21, 14,-35,  3, 16,-18, 10,-29,-17, 20, -5,  8,  0,  7,  4,-17,-42,-21], \
     [ 11,-28,-10, -6,-46, -5, -2,-20,-11,-25,-21, -5, 59,  2, -2,  9,  3,-12,-56,-49], \
     [ -4,-54, 16, 25,-47,-12, 29,-20,  7,-18,-10,  8,  2, 40, 13, -5, -8,-19,-48,-40], \
     [-15,-36,-13,-11,-45,-26, 16,-20, 34,-30, -4,  0, -2, 13, 61, -3, -9,-25, 22,-42], \
     [ 11,  0,  3,  0,-32, 11, -8,-14, -2,-28,-16,  7,  9, -5, -3, 16, 13,-10,-25,-28], \
     [ 12,-22, -1, -4,-31,  0,-13,  1,  0,-17, -6,  4,  3, -8, -9, 13, 26,  3,-52,-27], \
     [  2,-19,-21,-18,-12,-14,-22, 37,-24, 19, 18,-17,-12,-19,-25,-10,  3, 43,-62,-25], \
     [-58,-78,-68,-70,  4,-70,-28,-51,-35,-18,-42,-42,-56,-48, 22,-25,-52,-62,173, -2], \
     [-35,  3,-43,-43, 70,-52, -1, -9,-44, -9,-24,-21,-49,-40,-42,-28,-27,-25, -2,101]]
  (i,j) = (AAs.find(a),AAs.find(b))
  if i==-1 or j==-1: return 0 # should be mean value?
  return mdm78[i][j]

def ident_matches(alignment):
  matches = ""
  (a,b) = alignment
  for i in xrange(len(a)):
    if a[i]==b[i]: matches += "|"
    else: matches += ' '
  return matches

def dayhoff_matches(alignment):
  matches = ""
  (a,b) = alignment
  for i in xrange(len(a)):
    if a[i]==b[i]: matches += "|"
    elif dayhoff(a[i],b[i])>0: matches += '*'
    else: matches += ' '
  return matches

def exercise():
  A = "AAAGGTT"
  B = "AAATT"
  obj = align(A,B,verbose=True)

  print "score=%0.1f" % obj.score()
  alignment = obj.extract_alignment()
  print alignment[0]
  print ident_matches(alignment)
  print alignment[1]

  # 1rra vs. 1bli
  A = "AESSADKFKRQHMDTEGPSKSSPTYCNQMMKRQGMTKGSCKPVNTFVHEPLEDVQAICSQGQVTCKNGRNNCHKSSSTLRITDCRLKGSSKYPNCDYTTTDSQKHIIIACDGNPYVPVHFDASV"
  B = "DNSRYTHFLTQHYDAKPQGRDDRYCESIMRRRGLTSPCKDINTFIHGNKRSIKAICENKNGNPHRENLRISKSSFQVTTCKLHGGSPWPPCQYRATAGFRNVVVACENGLPVHLDQSIFRRP"
  obj = align(A,B,gop=150,gep=20,sim=dayhoff,style=align.LOCAL)

  print "\n1rra vs. 1bli"
  print "score=%0.1f" % obj.score()
  alignment = obj.extract_alignment()
  print alignment[0]
  print dayhoff_matches(alignment)
  print alignment[1]

  print "OK" # necessary for auto_build checking

if __name__=="__main__":

  exercise()
