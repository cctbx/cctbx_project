from scitbx import differential_evolution as de
from scitbx.array_family import flex
import math, sys
from libtbx.test_utils import approx_equal


class linear_scaler(object):
  """This class scales together varios curves. means and varainces should be lists of flex arrays"""
  def __init__(self, means, variances, reference_id=0,spread=3.0, factor=50,f=0.7,eps=1e-12,out=None,show_progress=False,insert_solution_vector=None):
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
    self.domain = [ ( -spread,spread ) ]*self.n
    self.optimizer =  de.differential_evolution_optimizer(self,population_size=self.n*factor,show_progress=show_progress,eps=eps, f=f,n_cross=2,cr=0.8,insert_solution_vector=insert_solution_vector)

  def setup_coeff_map(self):
    """ ii: datset number ; result[ii]=associated scale_factor index """
    result = []
    count = 0
    for ii in xrange(self.n_sets):
      if ii != self.ref:
        result.append( count )
        count += 1
      else:
        result.append( None )
    return result

  def get_mean(self,scales,offsets):
    if self.n_sets>2:
      result = self.means[self.ref]*0
      weights = self.means[self.ref]*0
      for m,v,s,o in zip(self.means,self.vars,scales,offsets):
        mm = s*(m+o)
        w = 1/(flex.sqrt(v)+1e-8)
        weights += w
        result = result + mm*w
      result = result/weights
      return result
    else:
      return self.means[self.ref]

  def get_scales_offsets(self,vector):
    scales = []
    offsets = []
    for scale in xrange(self.n_sets):
      tmp_scale = 1.0
      tmp_offset = 0.0
      if scale != self.ref:
        tmp_scale = abs(vector[self.map[scale]])
        tmp_offset= vector[self.map[scale]+self.n_sets-1]
      scales.append( tmp_scale )
      offsets.append( tmp_offset )
    return scales , offsets

  def target(self,vector):
    scales, offsets = self.get_scales_offsets(vector)
    dr = self.get_mean(scales,offsets)
    result = 0
    for jj in xrange(self.n_sets):
      dj =  scales[jj]*(self.means[jj]+offsets[jj])
      vj =  self.vars[jj]*scales[jj]*scales[jj]
      t  =  flex.pow((dj-dr),2)/( 1e-13 + vj )
      if self.n_sets != 2:
        result += flex.sum( t )
      else:
        if jj != self.ref:
          vr = self.vars[self.ref]
          result += flex.sum(flex.pow((dj-dr),2)/( 1e-13 + vj + vr ))
    return result

  def print_status(self, best,mean,vector,count):
    scales,offsets = self.get_scales_offsets(vector)
    print >> self.out, "PROGRESS"
    print >> self.out, count, best, mean
    print >> self.out, scales
    print >> self.out, offsets
    print >> self.out


  def retrieve_results(self):
    scales,offsets = self.get_scales_offsets( self.x )
    return scales,offsets,self.map

def scale_it_global(m,v,ref_id=0,factor=10,show_progress=False):
  scaler_object = linear_scaler( m,v,ref_id,factor=factor,show_progress=show_progress)
  s,o,m = scaler_object.retrieve_results()
  return s,o

def scale_it_pairwise(m,v,ref_id=0,factor=10,show_progress=False):
  n = len(m)
  scales = []
  ofsets = []
  for ii in xrange(n):
    if ii != ref_id:
      ms = [ m[ref_id],m[ii] ]
      vs = [ m[ref_id],m[ii] ]
      scaler_object = linear_scaler( ms,vs,0,factor=factor,show_progress=show_progress)
      s,o,map = scaler_object.retrieve_results()
      scales.append( s[1] )
      ofsets.append( o[1] )
    else:
      scales.append( 1.0 )
      ofsets.append( 0.0 )
  return scales, ofsets

def test_curve_scaler():
   flex.set_random_seed( 12345 )
   x = flex.double( range(10) )/10.0
   y = flex.exp( -x )
   y0 = y
   y1 = y*100+50
   y2 = y*1000+500
   v0 = y0*0+1
   v1 = y0*0+1.0
   v2 = v1

   scaler_object = linear_scaler( [y0,y1,y2], [v0,v1,v2], 0 , factor=5,show_progress=False)
   s,o,m = scaler_object.retrieve_results()
   assert approx_equal(s[1],0.01,eps=1e-3)
   assert approx_equal(s[2],0.001,eps=1e-3)
   assert approx_equal(o[1],-50,eps=1e-2)
   assert approx_equal(o[2],-500,eps=1e-2)

   s,o = scale_it_pairwise([y0,y1,y2],[v0,v1,v2],0, show_progress=False)
   assert approx_equal(s[1],0.01,eps=1e-3)
   assert approx_equal(s[2],0.001,eps=1e-3)
   assert approx_equal(o[1],-50,eps=1e-2)
   assert approx_equal(o[2],-500,eps=1e-2)


if __name__ == "__main__":
  test_curve_scaler()
  print "OK"
