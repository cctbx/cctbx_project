from scitbx.array_family import flex
from scitbx.math import correlation
from stdlib import math as smath
from scitbx import fftpack, simplex

def get_mean_sigma( nlm_array ):
  coef = nlm_array.coefs()
  mean = abs( coef[0] )
  var = flex.sum( flex.norm(coef) )
  sigma = smath.sqrt( var-mean*mean )
  return mean, sigma


class align(object):
  def __init__( self, fixed, moving, nmax=10, n_beta=21, ngrid=21, topn=10, refine=False, check_inversion=False,  show_result=False ):
    self.nmax = nmax
    self.fixed = fixed
    self.moving = moving
    self.beta = smath.pi*2.0/float(n_beta-1)*flex.double( range(n_beta) )
    self.nb = n_beta
    self.pad = max(0, (ngrid-1)//2 - nmax )
    self.ngrid = (self.pad+nmax) * 2 + 1
    self.dx = smath.pi*2.0/(self.ngrid*10)
    self.topn = topn
    self.refine = refine
    self.check_inversion = check_inversion
    self.inversion = False
    self.show_result=show_result
    self.top_align=[]
    self.top_scores = flex.double()
    self.scores = flex.double()
    self.cc_obj = correlation( fixed, moving, nmax, 0 ) # make default beta=0
    self.scan()
    ea = self.best_ea
    self.moving_nlm = self.cc_obj.rotate_moving_obj( ea[0],ea[1], ea[2], self.inversion )

  def get_cc( self ):
    fix_mean, fix_s = get_mean_sigma( self.fixed )
    mov_mean, mov_s = get_mean_sigma( self.moving)
    self.cc = ( self.best_score - fix_mean*mov_mean ) / ( fix_s*mov_s )
    return self.cc

  def scan( self ):
    fft = fftpack.complex_to_complex_2d( self.ngrid, self.ngrid )
    inversion = False
    for beta in self.beta:
      self.cc_obj.set_beta( beta )
      mm = self.cc_obj.mm_coef(0,inversion)
      if( self.pad > 0):
        mm = self.cc_obj.mm_coef(self.pad, inversion)
      fft_input= mm
      scores = fft.backward( fft_input ).as_1d()
      self.scores = self.scores.concatenate( -flex.norm( scores )  )
    self.best_indx = flex.min_index( self.scores )
    self.best_score = smath.sqrt( -self.scores[ self.best_indx ])


    if self.check_inversion:
    ### Inversion of the Spherical Harmonics ###
      inversion = True
      inversion_scores = flex.double()
      for beta in self.beta:
        self.cc_obj.set_beta( beta )
        mm = self.cc_obj.mm_coef(0,inversion)
        if( self.pad > 0):
          mm = self.cc_obj.mm_coef(self.pad, inversion)
        fft_input= mm.deep_copy()
        scores = fft.backward( fft_input ).as_1d()
        inversion_scores = inversion_scores.concatenate( -flex.norm( scores )  )
      inv_best_indx = flex.min_index( inversion_scores )
      inv_best_score = smath.sqrt(-inversion_scores[ inv_best_indx ] )

      if( inv_best_score < self.best_score ):
        self.score = inversion_scores
        self.best_indx = inv_best_indx
        self.best_score = inv_best_score
        self.inversion =  True
      else:
        self.inversion = False



    b=self.best_indx//(self.ngrid*self.ngrid)
    a=(self.best_indx - self.ngrid*self.ngrid*b ) // self.ngrid
    g=self.best_indx - self.ngrid*self.ngrid*b - self.ngrid*a

    b = self.beta[b]
    g = smath.pi*2.0 *( float(g)/(self.ngrid-1) )
    a = smath.pi*2.0 *( float(a)/(self.ngrid-1) )

    self.best_ea = (a, b, g )

    self.find_top( self.topn )
    if( self.refine ):
      self.refined = []
      self.refined_score = flex.double()
      for t in self.top_align:
        r = self.run_simplex( t )
        self.refined.append ( r )
        self.refined_score.append( self.target( r ) )

      orders=flex.sort_permutation( self.refined_score )
      self.best_score = -self.refined_score[orders[0]]


# show the refined results
      if( self.show_result ):
        print "refined results:"
        for ii in range( self.topn ):
          o = orders[ii]
          o = ii
          print ii, ":", list( self.refined[o] ), ":", self.refined_score[o]
      ea = self.refined[ orders[0] ]
      self.best_ea = (ea[0], ea[1], ea[2] )
      self.moving_nlm = self.cc_obj.rotate_moving_obj( ea[0],ea[1], ea[2], self.inversion )


  def find_top( self, topn ):
    orders = flex.sort_permutation( self.scores )
    for ii in range( topn ):
      o = orders[ii]
      b=o//(self.ngrid*self.ngrid)
      a=(o - self.ngrid*self.ngrid*b ) // self.ngrid
      g=o - self.ngrid*self.ngrid*b - self.ngrid*a

      b = self.beta[b]
      g = smath.pi*2.0 *( float(g)/(self.ngrid-1) )
      a = smath.pi*2.0 *( float(a)/(self.ngrid-1) )
      self.top_align.append( flex.double( (a, b, g) ) )
      self.top_scores.append( self.scores[o] )
      #print ii, ":", a, b, g, ":", self.scores[o]

  def run_simplex( self, start, max_iter=500 ):
    dim = 3
    starting_matrix = [ start ]
    for ii in range( dim ):
      starting_matrix.append( start + (flex.random_double(dim)*2-1)*self.dx )
    optimizer = simplex.simplex_opt( dimension = dim,
                                     matrix = starting_matrix,
                                     evaluator = self,
                                     max_iter = max_iter,
                                     tolerance=1e-5)
    result = optimizer.get_solution()
    return result

  def target( self, v ):
    cc = self.cc_obj.calc_correlation( v[0], v[1], v[2], self.inversion )
    return -abs( cc )


def tst():
  xyz = flex.vec3_double(  [ (-1,-1,-1) , (1,1,1) ] )
