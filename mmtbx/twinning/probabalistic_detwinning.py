from cctbx import maptbx
from cctbx import miller
from cctbx import crystal
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import xray
from cctbx import eltbx
from cctbx import adptbx
from scitbx import lbfgs
from mmtbx import masks
from libtbx import adopt_init_args
from cctbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
from iotbx.pdb import xray_structure
import mmtbx.scaling
import scitbx.math as sm
from mmtbx.scaling import absolute_scaling, relative_scaling
from mmtbx.scaling import matthews, twin_analyses
from mmtbx.scaling import basic_analyses, pair_analyses
from mmtbx import masks
from libtbx import table_utils
import scitbx.lbfgs
import libtbx.phil.command_line
from cStringIO import StringIO
from scitbx import differential_evolution
import sys, os, math, time
from libtbx.test_utils import approx_equal


class detwin(object):
  def __init__(self,
               i1,
               s1,
               i2,
               s2,
               alpha,
               eps=1e-13):
    assert s1 > 0
    assert s2 > 0
    self.i1=i1
    self.s1=s1
    self.i2=i2
    self.s2=s2
    self.a=alpha
    self.eps=eps
    self.x=flex.double( [math.sqrt(math.fabs(i1)),
                         math.sqrt(math.fabs(i2))] )
    self.minimizer = lbfgs.run(
      target_evaluator=self)
    self.xm =self.x[0]
    self.ym =self.x[1]
    self.vcv=self.snd(self.xm,self.ym)

  def compute_functional_and_gradients(self,h=0.00000001):
    f = -self.log_p(self.x[0], self.x[1])
    g = -self.d_log_p(self.x[0], self.x[1])
    return f,g

  def log_p(self, xm, ym):
    if xm<=0:
      xm=self.eps
    if ym<=0:
      ym=self.eps

    tmp1 = math.log(xm) + math.log(ym)
    a =  self.a
    ic1 = (1-a)*xm*xm + a*ym*ym
    ic2 = (1-a)*ym*ym + a*xm*xm
    tmp2 = self.i1-ic1
    tmp2 = tmp2*tmp2/(2.0*self.s1*self.s1)
    tmp3 = self.i2-ic2
    tmp3 = tmp3*tmp3/(2.0*self.s2*self.s2)
    return tmp1 - tmp2 - tmp3

  def d_log_p(self, xm, ym):
    if xm<=0:
      xm=self.eps
    if ym<=0:
      ym=self.eps

    a    = self.a
    df1 = -2*(-1+a)*xm*((-1+a)*xm*xm-a*ym*ym+self.i1)/(self.s1**2) +\
          -2*a*xm*(ym*ym+a*(xm-ym)*(xm+ym)-self.i2)/(self.s2**2) +\
          1.0/xm

    df2 = -2*a*ym*(-(-1+a)*xm*xm+a*ym*ym-self.i1)/(self.s1**2) +\
          -2*(-1+a)*ym*(-ym*ym+a*(-xm*xm+ym*ym)+self.i2)/(self.s2**2) +\
          1.0/ym

    return flex.double([df1,df2])

  def snd(self, xm, ym ):
    if xm<=0:
      xm=self.eps
    if ym<=0:
      ym=self.eps
    a = self.a
    sndx = 2*a*xm*xm*(3*a*xm*xm+ym*ym-a*ym*ym-self.i2)*self.s1*self.s1
    sndx+= (2*(-1+a)*xm*xm*(3*(-1+a)*xm*xm-a*ym*ym+self.i1)
            +self.s1*self.s1)*self.s2**2
    sndx = -sndx/(xm*xm*self.s1*self.s1*self.s2*self.s2)

    sndy = 2*(-1+a)*ym*ym*(-a*xm*xm-3*ym*ym+3*a*ym*ym+self.i2)*self.s1**2
    sndy+= (2*a*ym*ym*(-(-1+a)*xm*xm+3*a*ym*ym-self.i1)+self.s1**2)*self.s2**2
    sndy = -sndy/(ym*ym*self.s1*self.s1*self.s2*self.s2)

    dxy = 4*(-1+a)*a*xm*ym*(self.s1**2 +self.s2**2)/(
      (self.s1**2) * (self.s2**2))

    det = sndx*sndy-dxy*dxy
    s11 = sndy/det
    s22 = sndx/det
    s12 = -dxy/det
    return( (-s11,-s22,-s12) )


def detwin_miller_array(miller_obs,
                        twin_law,
                        twin_fraction):
  # for the moment, ignore incompleten twin pairs please
  cb_op = sgtbx.change_of_basis_op( twin_law )
  twin_related_miller  = miller_obs.change_basis( cb_op ).set_observation_type(
    miller_obs )
  set1, set2 = miller_obs.common_sets( twin_related_miller )\
               .set_observation_type(miller_obs )

  assert miller_obs.observation_type() is not None
  assert miller_obs.sigmas() is not None
  if set1.is_xray_amplitude_array():
    set1 = set1.f_as_f_sq()
    set2 = set1.f_as_f_sq()

  detwinned_f     = flex.double()
  detwinned_sigma = flex.double()
  if set1.is_xray_intensity_array():
    if set2.is_xray_intensity_array():
      for i1,s1,i2,s2 in zip( set1.data(), set1.sigmas(),
                              set2.data(), set2.sigmas() ):
        tmp_detwinner = detwin(i1,s1,i2,s2)
        # we do some double work here actually
        ni1 = tmp_detwinner.xm
        ns1 = math.sqrt( math.abs(self.vcv[0]) )
        detwinned_f.append( ni1 )
        detwinned_s.append( ns1 )

  set1 = set1.f_sq_as_f()
  new_f = set1.customized_copy


def test_detwin():
  # test the detwinning
  i1=3.5
  i2=2.5
  s1=0.01
  s2=0.01
  a=0.25
  tmp = detwin( i1,s1,i2,s2,a)
  assert approx_equal( tmp.xm, 2.0, eps=1e-4)
  assert approx_equal( tmp.ym, math.sqrt(2.0), eps=1e-4)

  # test the variance covariance matrix.
  # when the twin fraction is 0, we have two independent observations!
  i1=4.0
  i2=9.0
  s1=0.1
  s2=0.1
  a=0.0
  tmp = detwin( i1,s1,i2,s2,a)
  ttt = (i2*i2+2*s2*s2)**0.25
  ttt = s2/(2*ttt)
  assert approx_equal( tmp.vcv[1], ttt*ttt, eps=1e-5)
  ttt = (i1*i1+2*s1*s1)**0.25
  ttt = s1/(2*ttt)
  assert approx_equal( tmp.vcv[0], ttt*ttt, eps=1e-5)
  assert approx_equal( tmp.vcv[2], 0 , eps=1e-8)


  i1=4000.0
  i2=3500.0
  s1=10
  s2=12
  a=0.5
  tmp = detwin( i1,s1,i2,s2,a)
  assert approx_equal([tmp.xm, tmp.ym], [61.6048502913, 61.6038732295])

def run():
  test_detwin()

if __name__ == '__main__':
  run()
  print "OK"
