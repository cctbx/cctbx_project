import sys, os, time
from scitbx.array_family import flex
from scitbx.math import correlation
from stdlib import math as smath
from libtbx import easy_pickle
from scitbx import math, fftpack, simplex
import iotbx.phil

from sastbx.interface import get_input


master_params = iotbx.phil.parse("""\
align{
  fix = None
  .type=path
  .help="pickle Cnlm coef of the fixed object"

  mov = None
  .type=path
  .help="pickle Cnlm coef of the fixed object"

  num_grid = 41
  .type=int
  .help="number of point in each euler angle dimension"

  nmax = 20
  .type=int
  .help="maximum order of zernike polynomial:fixed for the existing database"

  topn = 10
  .type=int
  .help="top N alignments will be further refined if required"

  refine = False
  .type=bool
  .help="Refine the initial alignments or not"

  show_result = False
  .type=bool
  .help="show results to stdout or not"

}
""")

banner = "-----------------------Align two structures (zernike coefficients pickle file) ----------------"

class align(object):
  def __init__( self, fixed, moving, nmax=20, n_beta=41, ngrid=41, topn=10, refine=False, show_result=False ):
    self.nmax = nmax
    self.fixed = fixed
    self.moving = moving
    self.beta = smath.pi*2.0/float(n_beta)*flex.double( range(n_beta) )
    self.nb = n_beta
    self.na = nmax*2+1
    self.pad = (ngrid-1)/2 - nmax
    self.ngrid = (self.pad+nmax) * 2 + 1
    self.dx = smath.pi*2.0/self.ngrid
    self.topn = topn
    self.refine = refine
    self.show_result=show_result
    self.top_align=[]
    self.top_scores = flex.double()
    self.scores = flex.double()
    self.cc_obj = correlation( fixed, moving, nmax, 0 ) # make default beta=0
    self.scan()
    ea = self.best_ea
    self.moving_nlm = self.cc_obj.rotate_moving_obj( ea[0],ea[1], ea[2] )

  def scan( self ):
    fft_input = flex.complex_double( flex.grid(self.ngrid,self.ngrid) )
    fft = fftpack.complex_to_complex_2d( self.ngrid, self.ngrid )
    for beta in self.beta:
      self.cc_obj.set_beta( beta )
      mm = self.cc_obj.mm_coef(0)
      if( self.pad > 0):
        mm = self.cc_obj.mm_coef(self.pad)
      fft_input= mm
      scores = fft.backward( fft_input ).as_1d()
      self.scores = self.scores.concatenate( -flex.abs( scores )  )
    self.best_indx = flex.min_index( self.scores )
    self.best_score = self.scores[ self.best_indx ]

    b=self.best_indx/(self.ngrid*self.ngrid)
    a=(self.best_indx - self.ngrid*self.ngrid*b ) / self.ngrid
    g=self.best_indx - self.ngrid*self.ngrid*b - self.ngrid*a

    b = self.beta[b]
    g = smath.pi*2.0 *( float(g)/(self.na-1) )
    a = smath.pi*2.0 *( float(a)/(self.na-1) )

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
# show the refined results
      if( self.show_result ):
        print "refined results:"
        for ii in range( self.topn ):
          o = orders[ii]
          o = ii
          print ii, ":", list( self.refined[o] ), ":", self.refined_score[o]
      ea = self.refined[ orders[0] ]
      self.moving_refined = self.cc_obj.rotate_moving_obj( ea[0],ea[1], ea[2] )


  def find_top( self, topn ):
    orders = flex.sort_permutation( self.scores )
    for ii in range( topn ):
      o = orders[ii]
      b=o/(self.ngrid*self.ngrid)
      a=(o - self.ngrid*self.ngrid*b ) / self.ngrid
      g=o - self.ngrid*self.ngrid*b - self.ngrid*a

      g = smath.pi*2.0 *( float(g)/(self.na-1) )
      b = smath.pi*2.0 *( float(b)/(self.na-1) )
      a = smath.pi*2.0 *( float(a)/(self.na-1) )
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
    cc = self.cc_obj.calc_correlation( v[0], v[1], v[2] )
    return -abs( cc )


################## main function #################
def run(args):
  params = get_input(args, master_params, "align", banner, help)
  fix = params.align.fix
  mov = params.align.mov
  num_grid = params.align.num_grid
  nmax = params.align.nmax
  fixed = math.nlm_array(nmax)
  moving = math.nlm_array(nmax)
  coefs = easy_pickle.load( fix )
  fixed.load_coefs( fixed.nlm(), coefs )
  coefs = easy_pickle.load( mov )
  moving.load_coefs( moving.nlm(), coefs )
  base = mov.split('\n')[0].split('.')[0]

  align_obj = align( fixed, moving, nmax=nmax, ngrid=num_grid, refine=params.align.refine, show_result=params.align.show_result  )

  print align_obj.best_indx
  print align_obj.best_score
  print align_obj.best_ea

  easy_pickle.dump(base+'_aligned.nlm.pickle', align_obj.moving_nlm.coefs() )
  if(params.align.refine):
    easy_pickle.dump(base+'_refined.nlm.pickle', align_obj.moving_refined.coefs() )

def help():
  print "Usage: fft_align.py fix=fix_pickle_file mov=moving_pickle_file"

if __name__ == "__main__":
  t1 = time.time()
  args=sys.argv[1:]
  run(args)
  t2 = time.time()
  print "#OK","total time: ", t2-t1
