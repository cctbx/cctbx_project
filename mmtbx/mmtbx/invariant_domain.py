import iotbx.phil
import libtbx.phil.command_line
from libtbx.utils import Sorry
from mmtbx import amino_acid_codes
from mmtbx import alignment
import sys, os, math
from iotbx import pdb
from cctbx.array_family import flex
from scitbx.math import euler_angles_as_matrix
from scitbx.math import superpose, matrix
import random
from libtbx.test_utils import approx_equal



class find_domain(object):
  """This class tries to find pseudo invariant domain given a set of sites.
It requiers an empirical parameter named max_level to be set to a sensible
value. For large movements, havingf this value set to 2.5 is probably okai.
If the movements are more subtle, more problem dedicated software might be
the best course of action"""

  def __init__(self,
               set_a,
               set_b,
               initial_rms=0.5,
               max_level=2.0,
               overlap_thres=0.75,
               minimum_size=25):
    assert set_a.size() == set_b.size()
    assert set_a.size() > 3
    assert minimum_size <= set_a.size()
    assert minimum_size <= set_b.size()
    assert initial_rms < max_level

    self.matches = []
    self.set_a = set_a
    self.set_b = set_b
    self.n = set_a.size()
    self.max_iter = 10
    tmp_match = self.zipper(initial_rms,max_level)
    self.matches = self.process( tmp_match,
                                 overlap_thres=overlap_thres,
                                 minimum_size=minimum_size)

  def show(self, out=None):
    if out is None:
      out = sys.stdout
    count=1
    print >> out, "%i invariant-like domains have been found. "% len(self.matches)
    if len(self.matches)>0:
      print >> out, "Listing transformations below."
    else:
      print "Change settings in improve results."

    for item in self.matches:
      r = item[1]
      t = item[2]
      rmsd = item[3]
      n = item[4]
      print >> out
      print >> out, "Operator set %i: "%(count)
      print >> out, r.mathematica_form(label="r", one_row_per_line=True, format="%8.5f")
      print >> out
      print >> out, t.mathematica_form(label="t", format="%8.5f")
      print >> out
      print >> out, "rmsd: %5.3f; number of sites: %i"%(rmsd,n)
      print >> out
      count += 1


  def process(self,matches, overlap_thres=0.75,minimum_size=10):
    #convert the booleans iubnto doubles
    used = flex.bool( len(matches), False )
    tmp_matches = []
    final_matches = []
    done = False
    while not done:
      # find the largest domain please
      size, index = self.find_largest( matches, used )
      if size>minimum_size:
        final_matches.append( matches[index] )
      # find all other sequences that share more than n% similarity
      sims = self.find_similar_matches( matches[index], matches, used, overlap_thres )
      used = used.set_selected( sims, True )
      if used.count( True ) == used.size() :
        done = True
    return final_matches

  def find_similar_matches( self, target, matches, used, overlap_thres ):
    tmp_a =  flex.double(  self.set_a.size() , 0.0 ).set_selected( target[0], 1.0 )
    result = flex.size_t()
    for ii in xrange(len(matches) ):
      if not used[ii]:
        match = matches[ii][0]
        tmp_b = flex.double( self.set_a.size(), 0.0 ).set_selected( match, 1.0 )
        similar = flex.sum( tmp_a*tmp_b )/flex.sum( tmp_a )
        if similar > overlap_thres:
          result.append( ii )
    return( result )

  def pair_sites(self, r, t, cut_off):
    new_sites = r.elems*self.set_b+t.elems
    deltas = self.set_a - new_sites
    deltas = flex.sqrt( deltas.dot(deltas) )
    select = flex.bool( deltas < cut_off )
    tmp_a = self.set_a.select( select.iselection() )
    tmp_b = self.set_b.select( select.iselection() )
    return tmp_a, tmp_b, select

  def zipper(self,initial_rms, level):
    matches = []
    for jj in xrange( 1,self.n-1 ):
      ii = jj - 1
      kk = jj + 1
      #make triplets of sequence related sites
      xi = self.set_a[ii] ; xpi = self.set_b[ii]
      xj = self.set_a[jj] ; xpj = self.set_b[jj]
      xk = self.set_a[kk] ; xpk = self.set_b[kk]
      #get the lsq matrix
      ref = flex.vec3_double( [xi,xj,xk] )
      mov = flex.vec3_double( [xpi,xpj, xpk] )

      lsq = superpose.least_squares_fit(ref,mov)
      #here we have the rotation and translation operators
      r = lsq.r
      t = lsq.t
      rmsd = 10.0
      #we would like to know the rmsd on the coords used for superposition
      new_sites = lsq.other_sites_best_fit()
      deltas = ref - new_sites
      rmsd = deltas.rms_length()
      if rmsd < initial_rms:
        # please apply this rotation to the full set
        converged = False
        count=0
        match_size = 0
        previous_match_size = 0
        tmp_a = None
        tmp_b = None
        select = flex.bool()
        while not converged:
          previous_match_size = match_size
          tmp_a, tmp_b, select = self.pair_sites(r,t,level)
          #print count, tmp_a.size()
          match_size = tmp_a.size()
          if match_size <= previous_match_size:
            converged=True
            break
          if count>self.max_iter:
            converged=True
            break
          if tmp_b.size()>0:
            lsq = superpose.least_squares_fit(tmp_a,tmp_b)
            tmp_sites = lsq.other_sites_best_fit()
            rmsd = tmp_a.rms_difference(tmp_sites)
            r = lsq.r
            t = lsq.t
            count += 1
        if converged:

          matches.append( [select.deep_copy().iselection(),
                           r,
                           t,
                           rmsd,
                           select.deep_copy().iselection().size() ]  )
    return matches

  def find_largest(self,matches,used_flags=None):
    sizes = flex.double()
    if used_flags is None:
      used_flags = flex.bool( len(matches), False )
    for match in matches:
      sizes.append( match[0].size() )
    multi = flex.double( len(matches), 1 )
    multi = multi.set_selected( used_flags.iselection(), 0)
    sizes = sizes*multi
    max_size = flex.max( sizes )
    max_loc = flex.max_index( sizes )
    return max_size, max_loc

  def write_ensemble(self, pdb_inp, out):
    for index in xrange( len(self.matches) ):
      name = out+"_%s.pdb"%(index)
      file = open(name, "w")
      for ser, lbl, atom in zip( pdb_inp.atom_serial_number_strings(),
                                 pdb_inp.input_atom_labels_list(),
                                 pdb_inp.atoms() ):
        r = self.matches[index][1]
        t = self.matches[index][2]
        print >> file, iotbx.pdb.format_atom_record(
          record_name={False: "ATOM", True: "HETATM"}[atom.hetero],
          serial=int(ser),
          name=lbl.name(),
          altLoc=lbl.altloc(),
          resName=lbl.resname(),
          resSeq=lbl.resseq,
          chainID=lbl.chain(),
          iCode=lbl.icode(),
          site=r * matrix.col(atom.xyz) + t,
          occupancy=atom.occ,
          tempFactor=atom.b,
          segID=atom.segid,
          element=atom.element,
          charge=atom.charge)


def run(n=10):
  # make two random sets of sites please
  c = euler_angles_as_matrix([random.uniform(0,360) for i in xrange(3)])
  set_1 = flex.vec3_double(flex.random_double(n*3)*10-2)
  set_2 = flex.vec3_double(flex.random_double(n*3)*10-2)
  set_3 = tuple(c)*set_2

  set_a = set_1.concatenate( set_2 )
  set_b = set_1.concatenate( set_3 )

  tmp = find_domain(set_a,
                    set_b,
                    initial_rms=0.02,
                    max_level=0.03,
                    minimum_size=1)
  #tmp.show()
  assert len(tmp.matches)==2
  assert approx_equal(tmp.matches[0][3],0,eps=1e-5)
  assert approx_equal(tmp.matches[0][4],n,eps=1e-5)

if( __name__ == "__main__"):
  run(4)
  run(10)
  run(20)
  run(100)
  run(200)
  print "OK"
