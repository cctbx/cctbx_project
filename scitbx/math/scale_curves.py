from scitbx import differential_evolution as de
from scitbx.array_family import flex
import math, sys
from libtbx.test_utils import approx_equal

class linear_scaler(object):
  """This class scales together varios curves. means and varainces should be lists of flex arrays"""
  def __init__(self, means, variances, reference_id=0,factor=50,f=0.7,eps=1e-12,out=None,show_progress=False):
    self.out = None
    if self.out is None:
      self.out = sys.stdout

    self.means = means
    self.vars  = variances
    self.ref = reference_id
    self.n_sets = len( self.means)
    self.map = self.setup_coeff_map()

    self.n = (len(self.means)-1)*2
    self.x = None
    self.domain = [ ( -10,10 ) ]*self.n
    self.optimizer =  de.differential_evolution_optimizer(self,population_size=self.n*factor,show_progress=show_progress,eps=eps, f=f,n_cross=2,cr=0.8)

  def setup_coeff_map(self):
    result = []
    for ii in xrange(self.n_sets):
      if ii != self.ref:
        result.append( ii )
    return result

  def get_scales_offsets(self,vector):
    scales = []
    offsets = []
    for scale in xrange(self.n_sets-1):
      tmp_scale = abs(vector[scale])
      scales.append( tmp_scale )

    for offset in xrange(self.n_sets-1,self.n):
      offsets.append( vector[offset] )
    return scales , offsets

  def target(self,vector):
    scales, offsets = self.get_scales_offsets(vector)
    dr = self.means[self.ref]
    vr = self.vars[self.ref]
    result = 0
    for jj in xrange(self.n_sets-1):
      ii = self.map[jj]
      dj =  scales[jj]*(self.means[ii]+offsets[jj])
      vj =  self.vars[ii]*scales[jj]*scales[jj]
      t  =  flex.pow((dj-dr),2)/( 1e-13 + vj + vr)
      result += flex.sum( t )
    return result

  def print_status(self, best,mean,vector,count):
    scales,offsets = self.get_scales_offsets(vector)
    print >> self.out, count, best, mean
    print >> self.out, scales
    print >> self.out, offsets
    print >> self.out


  def retreive_results(self):
    scales,offsets = self.get_scales_offsets( self.x )
    return scales,offsets


def test_curve_scaler():
   x = flex.double( range(100) )/10.0
   y = flex.exp( -x )
   y0 = y
   y1 = y*100+50
   y2 = y*1000+500
   v0 = y0*0+1
   v1 = y0*0+1.0
   v2 = v1

   scaler_object = linear_scaler( [y0,y1,y2], [v0,v1,v2], 0 , factor=10,show_progress=False)
   s,o = scaler_object.retreive_results()
   assert approx_equal(s[0],0.01,eps=1e-3)
   assert approx_equal(s[1],0.001,eps=1e-3)
   assert approx_equal(o[0],-50,eps=1e-2)
   assert approx_equal(o[1],-500,eps=1e-2)

if __name__ == "__main__":
  test_curve_scaler()
  print "OK"
