from __future__ import absolute_import, division, print_function
from scitbx import differential_evolution as de
from scitbx.array_family import flex
import sys
from libtbx.test_utils import approx_equal
from six.moves import range
from six.moves import zip

class curve_interpolator(object):
  def __init__(self, start, stop, n_points=100):
    self.start    = start
    self.stop     = stop
    self.n_points = n_points
    self.target_x = flex.double( range(n_points) )/float(n_points-1)
    self.target_x = self.target_x*(self.stop-self.start)+self.start
    self.delta    = self.target_x[1]-self.target_x[0]

  def interpolate(self, x_array, y_array):
    index_array = []
    result_array = []

    start_index_user = None
    end_index_user = None

    start_index_target = None
    end_index_target = None

    user_min = flex.min( x_array )
    user_max = flex.max( x_array )

    for jj,x in enumerate(self.target_x):
      this_index = None
      break_again = False
      for index,this_x in enumerate(x_array):
        if this_x - x >= 0:
          if x >= user_min:
            if x <= user_max:
              this_index = index
              if start_index_user is None:
                 start_index_user = this_index
              if start_index_target is None:
                start_index_target = jj
              end_index_user = this_index
              end_index_target = jj
              break
      index_array.append( this_index )
      y = None
      if this_index is not None:
        if this_index == 0:
          y = self.two_point_interpolate( x,
                                          x_array[this_index  ], y_array[this_index  ],
                                          x_array[this_index+1], y_array[this_index+1] )
        elif this_index == len(x_array)-1:
          y = self.two_point_interpolate( x,
                                          x_array[this_index-1], y_array[this_index-1],
                                          x_array[this_index], y_array[this_index] )

        else:
          y = self.parabolic_interpolate( x,
                                          x_array[this_index-1], y_array[this_index-1],
                                          x_array[this_index  ], y_array[this_index  ],
                                          x_array[this_index+1], y_array[this_index+1] )

        result_array.append( y )


    n = len(result_array)
    x = flex.double(self.target_x[start_index_target:end_index_target+1])
    y = flex.double(result_array)
    return x,y,(start_index_user,end_index_user),(start_index_target,end_index_target)

  def two_point_interpolate(self,x, xo, fxo, xp, fxp):
    ph = x-xo
    h = xp-xo
    p = ph
    if ph >0:
      p = ph/h
    result = (1-p)*fxo + p*fxp
    return result

  def parabolic_interpolate(self, x, xm, fxm, xo, fxo, xp, fxp):
    result = fxm*( (x-xo)*(x-xp) )/ ( (xm-xo)*(xm-xp) ) + \
             fxo*( (x-xm)*(x-xp) )/ ( (xo-xm)*(xo-xp) ) + \
             fxp*( (x-xm)*(x-xo) )/ ( (xp-xm)*(xp-xo) )
    return result









class linear_scaler(object):
  """This class scales together varios curves. means and varainces should be lists of flex arrays"""
  def __init__(self, means, variances, reference_id=0,init_mean=1.0,spread=0.1, factor=50,f=0.7,eps=1e-12,out=None,show_progress=False,insert_solution_vector=None,add=True):
    self.out = None
    if self.out is None:
      self.out = sys.stdout

    self.add=add

    self.means = means
    self.vars  = variances
    self.ref = reference_id
    self.n_sets = len( self.means)
    self.map = self.setup_coeff_map()

    self.n = (len(self.means)-1)*2
    self.x = None
    self.domain = [ ( -spread+init_mean,spread+init_mean ) ]*self.n
    self.optimizer =  de.differential_evolution_optimizer(self,population_size=self.n*factor,show_progress=show_progress,eps=eps, f=f,n_cross=2,cr=0.8,insert_solution_vector=insert_solution_vector)

  def setup_coeff_map(self):
    """ ii: datset number ; result[ii]=associated scale_factor index """
    result = []
    count = 0
    for ii in range(self.n_sets):
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
    for scale in range(self.n_sets):
      tmp_scale = 1.0
      tmp_offset = 0.0
      if scale != self.ref:
        tmp_scale = abs(vector[self.map[scale]])
        tmp_offset= vector[self.map[scale]+self.n_sets-1]
        if not self.add:
          tmp_offset=0.0
      scales.append( tmp_scale )
      offsets.append( tmp_offset )
    return scales , offsets

  def target(self,vector):
    scales, offsets = self.get_scales_offsets(vector)
    dr = self.get_mean(scales,offsets)
    result = 0
    for jj in range(self.n_sets):
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
    print("PROGRESS", file=self.out)
    print(count, best, mean, file=self.out)
    print(scales, file=self.out)
    print(offsets, file=self.out)
    print(file=self.out)


  def retrieve_results(self):
    scales,offsets = self.get_scales_offsets( self.x )
    return scales,offsets,self.map

def scale_it_global(m,v,ref_id=0,factor=10,show_progress=False,add=True):
  scaler_object = linear_scaler( m,v,ref_id,factor=factor,show_progress=show_progress,add=add)
  s,o,m = scaler_object.retrieve_results()
  return s,o

def scale_it_pairwise(m,v,ref_id=0,factor=10,show_progress=False,add=True):
  n = len(m)
  scales = []
  ofsets = []
  for ii in range(n):
    if ii != ref_id:
      ms = [ m[ref_id],m[ii] ]
      vs = [ m[ref_id],m[ii] ]
      scaler_object = linear_scaler( ms,vs,0,factor=factor,show_progress=show_progress,add=add)
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



def tst_curve_interpolator():
  x = flex.double( range(25) )/24.0
  y = x*x
  ip = curve_interpolator(0,2.0,200)
  x_target = ip.target_x
  y_ref = x_target*x_target
  nx,ny,a,b = ip.interpolate(x,y)
  count = 0
  for xx in x_target:
    if flex.max(x) >= xx:
      count += 1
  assert count==len(nx)


  for yy,yyy in zip(ny,y_ref):
    assert approx_equal(yy,yyy,eps=1e-3)
  assert a[0]==0
  assert a[1]==24
  assert b[0]==0
  assert b[1] in (99,100)


  x = flex.double( range(5,23) )/24.0
  y = x*x
  ip = curve_interpolator(0,2.0,200)
  nx,ny,a,b = ip.interpolate(x,y)
  assert nx[0] >= flex.min(x)
  assert nx[-1] <= flex.max(x)
  y_ref= nx*nx
  for yy,yyy in zip(ny,y_ref):
    assert approx_equal(yy,yyy,eps=1e-3)


if __name__ == "__main__":
  test_curve_scaler()
  tst_curve_interpolator()
  print("OK")
