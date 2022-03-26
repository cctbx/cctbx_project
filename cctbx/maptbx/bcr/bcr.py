from __future__ import absolute_import, division, print_function
import math
from scipy.optimize import minimize
from libtbx import adopt_init_args

from cctbx.array_family import flex
from mmtbx.ncs import tncs

# Adopted by Pavel Afonine. This is a verbatim copy of the original code
# supplied by A. Urzhumtsev on 11/5/21, 08:19.

######################################################
#
#    decomposition of oscillating curves
#    (atomic images at a given resolution)
#                                by A.Urzhumtsev, L.Urzhumtseva, 2021
#
#    python remake of the respective fortran program
#
#    input curves are obtained by 'image-atom.py'
#
######################################################

def gauss(B, C, r, b_iso):
  assert B>0
  B = B + b_iso
  fpob = 4 * math.pi / B
  return C * fpob**1.5 * math.exp(-fpob*math.pi*r**2)

def chi(B, R, r, b_iso):
  assert B>0
  assert R>0
  assert r>=0
  B = B + b_iso
  mfpsob = (-4*math.pi**2)/B
  if(r==0):
    return ((4*math.pi/B)**1.5) * math.exp(mfpsob*R**2)
  else:
    e1 = math.exp(mfpsob*(r-R)**2)
    e2 = math.exp(mfpsob*(r+R)**2)
    return (4*math.pi*B)**(-0.5) * (1/(r*R)) * (e1-e2)

class calculator(object):

  def __init__(self, npeak,dens,yc,dist,nfmes, x):
    adopt_init_args(self, locals())
    self.x = flex.double(x)

  def target_and_gradients(self, x):
    self.x = x
    target = FuncFit(self.x, self.npeak, self.dens, self.yc, self.dist,
      self.nfmes)
    gradients = GradFit(self.x, self.npeak, self.dens, self.yc, self.dist,
      self.nfmes)
    return target, flex.double(gradients)

class minimizer_bound(object):
  def __init__(self,
               calculator,
               use_bounds,
               lower_bound,
               upper_bound,
               max_iterations):
    adopt_init_args(self, locals())
    self.x = self.calculator.x
    self.n = self.x.size()
    self.n_func_evaluations = 0
    self.max_iterations = max_iterations

  def run(self):
    self.minimizer = tncs.lbfgs_run(
      target_evaluator = self,
      use_bounds       = self.use_bounds,
      lower_bound      = self.lower_bound,
      upper_bound      = self.upper_bound,
      max_iterations   = self.max_iterations)
    self()
    return self

  def __call__(self):
    self.n_func_evaluations += 1
    f, g = self.calculator.target_and_gradients(x = self.x)
    self.f = f
    self.g = g
    return self.x, self.f, self.g

#-------------------------------------------------------------------------------

def PreciseData(kpres,epsc,peak,dstep,nfmes):

#   ceps - defines accuracy in absolute values
#   bmin - defines sharpest drop < 10**(-ndrop) of the 3D-exponential
#            at a distance = grid step
#   cmin - defines minimal peak height (64 majorates (4*pi)**1.5 )
#   epsp - drop when approximating internal peaks

    epsp = 0.001

    if kpres > 0:
      ceps=epsc
    else:
      ceps=peak*epsc

    ndrop = 2
    bmin  = 16.*dstep*dstep/ndrop
    cmin  = ceps*(bmin**1.5)/64.
    rmin  = 0.0

    if(nfmes is not None):
      print ('',file=nfmes)
      print ('estimated accuracy parameters :',file=nfmes)
      print ('   absolute max allowed error',ceps,file=nfmes)
      print ('   bmin, cmin :              ',bmin,cmin,file=nfmes)
      print ('---------------------------',file=nfmes)
      print ('',file=nfmes)

    return ceps,bmin,cmin,rmin,epsp

#============================
def PeakOrigin(dens,dist,epsp,bmin,cmin,nfmes):

    ndist = len(dist) - 1

#   invert the curve if negative peak

    d0 = dens[0]
    if d0 > 0.0:
      kneg=0
    else:
      kneg=1
      dens,d0 = CurveInvert(dens,d0)

#   parameters of the Gaussian peak in the origin

    b0,c0 = GaussBC(dist,epsp,bmin,cmin,dens,nfmes)

#   invert curve back if necessary

    if kneg >0 :
      dens,d0 = CurveInvert(dens,d0)

    return b0,c0

#============================
def CurveInvert(dens,d0):

    ndist = len(dens) - 1

    for i in range(ndist+1):
      dens[i] = -dens[i]

    d0 = -d0

    return dens,d0

#============================
def GaussBC(dist,epsp,bmin,cmin,curres,nfmes):

    ndist = len(dist) -1

    cpi   = math.pi
    pisq4 = 4.*cpi*cpi

#   inflection point and approximate curve till this point

    kpeak = 0
    minf = NumInflRight(curres,kpeak,epsp,nfmes)
    if minf == 0:
       b0 = bmin
       c0 = curres[kpeak] * (b0/4.0*math.pi)**1.5
    else:
       ug,vg = uvFitGauss(dist,curres,minf)
       b0 = pisq4/vg
       c0 = ug*(cpi/vg)**1.5
       if b0 < bmin:
          b0 = bmin

    if c0 < cmin:
       c0 = cmin

    return b0,c0

#============================
def RippleBCR(dist,kpeak,epsp,bmin,cmin,rmin,curres,nfmes):

#   inflection point and approximate curve till this point

    ndist = len(dist) -1

    minfr = NumInflRight(curres,kpeak,epsp,nfmes)
    minfl = NumInflLeft(curres,kpeak,epsp,nfmes)
    if minfr-kpeak > kpeak-minfl:
      kref = minfr
    else:
      kref = minfl

    if kref == kpeak:
       b0 = bmin
       c0 = curres[kpeak] * math.sqrt(4.0*math.pi*b0) * dist[kpeak]**2
       if c0 < cmin:
          c0 = cmin

    b0,c0,r0 = FitBCR(dist,kpeak,kref,epsp,bmin,cmin,rmin,curres,nfmes)

    return b0,c0,r0

#============================
def TailBCR(dist,epsp,bmin,cmin,rmin,curres,nfmes):

#   inflection point for the peak at the upper interval bound

    ndist = len(dist) -1

    kref = NumInflLeft(curres,ndist,epsp,nfmes)

#   parameters for this peak

    b0,c0,r0 = FitBCR(dist,ndist,kref,epsp,bmin,cmin,rmin,curres,nfmes)

    return b0,c0,r0

#============================
def FitBCR(dist,kpeak,kref,epsp,bmin,cmin,rmin,curres,nfmes):

    ndist = len(dist) - 1

#   BCR parameters for an internal ripple

    pi4   = 4.  * math.pi
    pisq4 = pi4 * math.pi

#   scorr is an artificial correcting factor against B overestimation
#   when r0 is taken without correction B/(8*peak*pi**2)

    scorr=1.5

    if kref == kpeak :
       q0, r0  = curres[kpeak] , dist[kpeak]
       b0 = bmin

    else :

       q2, r2 = curres[kpeak] , dist[kpeak]
       q1, r1 = curres[kref]  , dist[kref]

       arg = (q2*r2) / (q1*r1)
       b2  = pisq4 * (r1-r2)**2 /math.log(arg)
       if b2 < bmin:
          b2 = bmin

#   correction of the peak center

       r2c = r2 + b2/(scorr*2.*pisq4*r2)
       q0, r0  = curres[ndist] , dist[ndist]

       for k in range(kpeak,ndist+1):
           if dist[k] >= r2c or curres[k] <= q1:
              q0, r0 = curres[k-1] , dist[k-1]
              break

       arg = (q0*r0) / (q1*r1)
       b0  = pisq4 * (r1-r0)**2 / math.log(arg)
       if b0 < bmin:
          b0 = bmin

    c0 = q0 * math.sqrt(pi4*b0) * r0**2
    if c0 < cmin:
       c0 = cmin

    return b0,c0,r0

#============================
def NumInflRight(curres,kpeak,epsp,nfmes):

#   right inflection point ; always 0 <= kpeak < ndist

    ndist = len(curres) - 1

    clim = epsp * curres[kpeak]

    k = ndist
    if kpeak < ndist-1 :
       for i in range (kpeak,ndist-1) :
           if curres[i+1] <= clim or curres[i+2]+curres[i]-2.*curres[i+1] >= 0.0 :
              k = i
              break

    if k < ndist or curres[ndist] > clim :
       return k
    else :
       return ndist-1

#============================
def NumInflLeft(curres,kpeak,epsp,nfmes):

#   left inflection point; always 0 < kpeak <= ndist ;
#   return 0 (r=0) is forbidden giving 0 in fitBCR

    ndist = len(curres) -1

    clim = epsp * curres[kpeak]

    k = 1
    if kpeak > 1:
       for i in range (kpeak,1,-1):
           if curres[i-1] <= clim or curres[i-2]+curres[i]-2.*curres[i-1] >= 0.0 :
              k = i
              break

    return k

#============================
def uvFitGauss(dist,curres,minf):

#   used below : w = ln(u)

    sumr2, sumr4, sumdl, sumrdl = 0.0 , 0.0 , 0.0 , 0.0

    for i in range(minf+1):
      curvei = math.log(curres[i])
      ri2    = dist[i] * dist[i]
      ri4    = ri2     * ri2
      sumr2  = sumr2  + ri2
      sumr4  = sumr4  + ri4
      sumdl  = sumdl  + curvei
      sumrdl = sumrdl + curvei * ri2

    det  = -(minf+1) * sumr4  + sumr2 * sumr2
    detw =    -sumdl * sumr4  + sumrdl*sumr2
    detv =  (minf+1) * sumrdl - sumr2 * sumdl

    wg = detw / det
    vg = detv / det
    ug = math.exp(wg)

    return ug,vg

#============================
def CurveDiff(dens,dist,nfmes,bpeak,cpeak,rpeak,npeak):

    ndist = len(dist) - 1

    curve  = [ [ 0 for j in range(npeak+1)] for i in range(ndist+1) ]
    curres = [ dens[i] for i in range(ndist+1) ]

#   cycle over components

    for ipeak in range(npeak+1):

#     select the  parameters

      b0, c0, r0 = bpeak[ipeak] , cpeak[ipeak] , rpeak[ipeak]

#     calculate the contribution
#         ripple in point r0 or Gaussian in the origin

      if r0 > 0.0 :
        tcurve = CurveRipple(dist,b0,c0,r0)
      else:
        tcurve = CurveGauss(dist,b0,c0)

#     remove contribution

      for ix in range(ndist):
        curve [ix][ipeak] = tcurve[ix]
        curres[ix]        = curres[ix]-tcurve[ix]

#   end of cycle over components; residual error

    epsres=0.0
    for ix in range (ndist):
      curabs=abs(curres[ix])
      if curabs > epsres:
         epsres = curabs

    if(nfmes is not None):
      print('',file=nfmes)
      #print(f'with {npeak+1:4} peaks modeled max residual peak is {epsres:12.7f}',
      #      file=nfmes)

    return curres,curve,epsres

#============================
def CurveGauss(dist,b0,c0):

    ndist = len(dist) -1

    argmax = 30.

#   curve for a Gaussian in the origin

    curve = [ 0.0 for i in range(ndist+1) ]

    cpi   = math.pi
    pisq4 = 4.*cpi*cpi

    vg = pisq4 / b0
    ug = c0 * (vg/cpi)**1.5

    for ix in range(ndist+1):
      arg = vg * dist[ix]**2
      if arg < argmax :
        curve[ix] = ug * math.exp(-arg)

    return curve

#============================
def CurveRipple(dist,b0,c0,r0):

    ndist = len(dist) -1

    argmax = 30.

#   curve for a rippe in point r0

    curve = [ 0.0 for i in range(ndist+1) ]

    cpi   = math.pi
    pi4   = 4.*cpi
    pisq4 = pi4*cpi

    scaler = c0 / (r0*math.sqrt(pi4*b0))
    coef   = pisq4 / b0

    arg0 = coef * r0**2
    if arg0 < argmax :
      curve[0] = c0 * math.exp(-arg0) * (pi4/b0)**1.5

    for ix in range(1,ndist+1):
      r = dist[ix]
      argm, argp = coef * (r-r0)**2 , coef * (r+r0)**2

      if argm < argmax :
        eargm = math.exp(-argm)
      else:
        eargm = 0.0

      if argp < argmax :
        eargp = math.exp(-argp)
      else:
        eargp = 0.0

      curve[ix] = scaler * (eargm-eargp) / r

    return curve

#============================
def OutCurves(dens,dist,curres,curve,filres):

    nfcurv = open(filres, 'w')

    ndist = len(dist) - 1

    #if(nfmes is not None):
    #  print(f'    dist       curve         modeled         resid',file=nfcurv)
    for i in range(ndist):
        densum = dens[i]-curres[i]
        #if(nfmes is not None):
        #  print(f'{dist[i]:8.3f}{dens[i]:15.8f}{densum:15.8f}{curres[i]:15.8f}',file=nfcurv)

    return

#============================
def OutCoefficients(bpeak,cpeak,rpeak,npeak,nfmes):
    if(nfmes is None): return
    space=' '

    print(' final parameter values:',file=nfmes)
    print(' Npeak     Bpeak',space*13,'Cpeak',space*13,'Rpeak',file=nfmes)
    #for ipeak in range(npeak+1):
    #   print(f'{ipeak:6}{bpeak[ipeak]:20.14f}{cpeak[ipeak]:20.14f}{rpeak[ipeak]:20.14f}',
    #         file=nfmes)

    return

#============================
def MorePeaks(curres,dist,nfmes,ceps,bpeak,cpeak,rpeak,npeak,epsp,bmin,cmin,rmin,mxp):

#   first iteration: the peak in the origin has been already treated
#   get starting point in the curve for further search

    if npeak == 0:
       kpeak=1
    else:
       kpeak=0

    ndist = len(dist) - 1

#   cycle over peaks

    for ipeak in range(npeak+1,mxp+1):

#       next maximal residual peak if yet significant

        kpeak,kneg = NextPeak(curres,kpeak,nfmes,ceps)

#       no more significant peak found

        if kpeak < 0:
           break

#       peak in the origin - model it by a Gaussian

        elif kpeak == 0:
           b0 , c0 = GaussBC(dist,epsp,bmin,cmin,curres,nfmes)
           r0 = 0.0

#       internal peak

        elif kpeak < ndist:
          b0, c0, r0 = RippleBCR(dist,kpeak,epsp,bmin,cmin,rmin,curres,nfmes)

#       peak at the right border of the interval

        else:
          b0, c0, r0 = TailBCR(dist,epsp,bmin,cmin,rmin,curres,nfmes)

        if kneg > 0:
          curres, c0 = CurveInvert(curres,c0)

        bpeak[ipeak] , cpeak[ipeak] , rpeak[ipeak] = b0 , c0, r0
        kpeak = kpeak+1

    if kpeak < 0:
       npeak = ipeak - 1
    else:
       npeak = ipeak

    return bpeak,cpeak,rpeak,npeak

#============================
def NextPeak(curres,kpeak,nfmes,ceps):

    ndist = len(curres) -1

#   check the peak in the origin

    if kpeak == 0:

#     check if the peak in the origin if flat

      cx = curres[0]
      i  = 1
      while curres[i]-cx == 0.0:
         i = i+1

      if cx*(curres[i]-cx) < 0.0 and abs(cx) > ceps:
         ix = 0

#     there is no peak in the origin; check next points

      else:
         kpeak = 1

#   (we cannot use 'else' below since kpeak may change after its previous check)

    if kpeak > 0:

#     already checked the whole interval and no significant peak found

      if kpeak >= ndist:
         ix = -1
         cx = 0.0

#     search for an internal peak

      else:
         kfound = 0
         for ix in range(kpeak,ndist):
             cx   = curres[ix]
             delm = cx-curres[ix-1]
             delp = curres[ix+1]-cx
             if delm*delp <= 0. and delp*cx < 0. and abs(cx) > ceps:
                kfound = 1
                break

#       no internal peak found ; check the upper interval bound

         if kfound == 0 and abs(curres[ndist]) > ceps:
            ix = ndist
            cx = curres[ndist]

#   flip the curve temporarily if the peak is negative

    kpeak = ix
    if cx >= 0.0:
       kneg = 0
    else:
       kneg = 1
       curres,cx = CurveInvert(curres,cx)

    return kpeak,kneg

#============================
def RefineBCR(dens,dist,bpeak,cpeak,rpeak,npeak,bmin,cmin,rmin,nfmes):

    ndist = len(dist) - 1

    xc        = [ 0 for i in range(3*(npeak+1)) ]
    yc        = [ 0 for i in range(ndist+1) ]
    bcrbounds = [ (None,None) for i in range(3*(npeak+1)) ]

    sp = ' '

#   define bounds
#   cmin is the bound for |c| and not for c itself

    for i in range(npeak+1):
        i3 = i*3
        xc[i3]   = bpeak[i]
        xc[i3+1] = cpeak[i]
        xc[i3+2] = rpeak[i]
        bcrbounds[i3]   = (bmin,1000.)
        if cpeak[i] > 0.0 :
           bcrbounds[i3+1] = (cmin,1000.)
        else :
           bcrbounds[i3+1] = (-1000.,-cmin)
        bcrbounds[i3+2] = (rmin,dist[ndist]*1.2)

#   minimization

    for it in range(1,3):
      if it > 1:
        xc = res.x
      CALC = calculator(npeak,dens,yc,dist,nfmes, xc)
      lbound = []
      ubound = []
      for b in bcrbounds:
        lbound.append(b[0])
        ubound.append(b[1])
      res = minimizer_bound(
        calculator  = CALC,
        use_bounds  = 2,
        lower_bound = lbound,
        upper_bound = ubound,
        max_iterations = 500).run().calculator
    res.x = list(res.x)

    if 0: # SciPy analogue. Works with Python 3 only.
      res = minimize(FuncFit,xc,args=(npeak,dens,yc,dist,ndist,nfmes),
                 method='L-BFGS-B',jac=GradFit, bounds = bcrbounds)

#   recover refined values
    if(nfmes is not None):
      print(' parameters before and after their refinement:',file=nfmes)
      print(' Npeak   Bpeak',sp*9,'Cpeak',sp*9,'Rpeak',sp*6,
                     'Bpeak',sp*9,'Cpeak',sp*9,'Rpeak',file=nfmes)

    for i in range(npeak+1):
        i3 = i*3
        b0 , c0, r0 = res.x[i3] , res.x[i3+1] , res.x[i3+2]

        #if(nfmes is not None):
        #  print(f'{i:6}{bpeak[i]:16.10f}{cpeak[i]:16.10f}{rpeak[i]:10.5f}   ',
        #        f'{b0:16.10f}{c0:16.10f}{r0:10.5f}',file=nfmes)

        bpeak[i] , cpeak[i] , rpeak[i] = b0 , c0 , r0

    return bpeak,cpeak,rpeak

#============================
def FuncFit(xc,npeak,y0,yc,dist,nfmes):

#   calculate model curve

    yc = CurveCalc(xc,npeak,yc,dist,nfmes)

#   comparer model and control curves

    fvalue = LSFunc(yc,y0)

    return fvalue

#============================
def GradFit(xc,npeak,y0,yc,dist,nfmes):

#   calculate the gradient with respect to the curve yc

    gy = GradCurve(yc,y0)

#   calculate the gradient with respect to the parameters

    gx = GradX(xc,npeak,yc,gy,dist,nfmes)

    return gx

#============================
def CurveCalc(xc,npeak,yc,dist,nfmes):

    ndist = len(dist) - 1

    for ix in range(ndist+1):
        yc[ix] = 0.0

    for i in range(npeak+1):
        i3 = i*3
        b0 , c0 , r0 = xc[i3] , xc[i3+1] , xc[i3+2]

#       contribution from a ripple or from a Gaussian in the origin

        if r0 > 0.0 :
           tcurve = CurveRipple(dist,b0,c0,r0)
        else:
           tcurve = CurveGauss(dist,b0,c0)

        for ix in range(ndist):
            yc[ix] = yc[ix] + tcurve[ix]

    return yc

#============================
def LSFunc(yc,y0):

    fvalue = 0.0
    ndist = len(yc) -1

    for i in range(ndist+1):
        fvalue = fvalue + (yc[i]-y0[i])**2

    fvalue = fvalue * 0.5

    return fvalue

#============================
def GradCurve(yc,y0):

    ndist = len(yc) -1

    gy = [ yc[i]-y0[i] for i in range(ndist+1) ]

    return gy

#============================
def GradX(xc,npeak,yc,gy,dist,nfmes):

    ndist = len(dist) -1

    gx = [ 0 for i in range(3*(npeak+1)) ]

    for i in range(npeak+1):
       i3 = i * 3
       b0 , c0, r0 = xc[i3] , xc[i3+1] , xc[i3+2]

       if r0 > 0.0:
          gb, gc, gr = GradRipple(yc,gy,dist,b0,c0,r0)
       else:
          gb, gc = GradGauss(yc,gy,dist,b0,c0)
          gr    = 0.0

       gx[i3] , gx[i3+1] , gx[i3+2]  = gb , gc , gr

    return gx

#============================
def GradRipple(yc,gy,dist,b0,c0,r0):

    ndist = len(dist) -1

    argmax = 30.

    cpi    = math.pi
    pisq4b = 4.0 * cpi * cpi/b0
    pisq8r = pisq4b * 2.0 * r0

    scc   = 1. / (2.0*r0*math.sqrt(b0*cpi))
    c0r0 , c0b0 = c0/r0 , c0/b0
    c0b02       = c0b0 / 2.0

#   gradient in the origin, x = 0

    arg0 = pisq4b * r0*r0
    if arg0 < argmax:
       gadd = 8.0 * math.exp(-arg0) * gy[0] * (cpi/b0)**1.5
       gb , gc , gr  =  gadd * c0b0 * (arg0 - 1.5) , gadd , -gadd * c0 * pisq8r
    else:
       gb , gc , gr = 0.0 , 0.0 , 0.0

#   gradient in internal points

    for ix in range (1,ndist):
        rx   = dist[ix]
        sccx = scc/rx
        gi   = gy[ix]

        xm   = rx - r0
        argm = pisq4b * xm*xm
        if argm < argmax:
           gadd = math.exp(-argm) * gi * sccx
           gb   = gb + gadd * c0b02 * (2.0*argm-1.0)
           gc   = gc + gadd
           gr   = gr + gadd * c0r0 * (pisq8r*xm-1.0)

        xp   = rx + r0
        argp = pisq4b * xp*xp
        if argp < argmax:
           gadd = math.exp(-argp) * gi * sccx
           gb   = gb - gadd * c0b02 * (2.0*argp -1.0)
           gc   = gc - gadd
           gr   = gr + gadd * c0r0  * (pisq8r*xp+1.0)

    return gb,gc,gr

#============================
def GradGauss(yc,gy,dist,b0,c0):

    ndist = len(dist) -1

    argmax = 30.

    cpi = math.pi

    vg     = 4.0 * cpi / b0
    ugc    = vg**1.5
    ugb    = ugc * c0/b0
    pisq4b = vg  * cpi

    gb , gc = 0.0 , 0.0

    for ix in range(ndist+1):
        xsq = dist[ix]**2
        arg = pisq4b * xsq
        if arg < argmax:
           gadd = math.exp(-arg) * gy[ix]
           gb   = gb + ugb*gadd*(arg-1.5)
           gc   = gc + ugc*gadd

    return gb,gc

#============================
def get_BCR(dens,dist,mxp,epsc,kpres,nfmes=None):

    bpeak  = [ 0 for i in range(mxp+1) ]
    cpeak  = [ 0 for i in range(mxp+1) ]
    rpeak  = [ 0 for i in range(mxp+1) ]

#
#   definition / estimations of the internal parameters
#
    ceps,bmin,cmin,rmin,epsp = PreciseData(kpres,epsc,dens[0],dist[1],nfmes)

#--------------------------------------
#
#   treat the principal peak in the origin - model it by a Gaussian
#   and remove its contribution
#

    bpeak[0],cpeak[0] = PeakOrigin(dens,dist,epsp,bmin,cmin,nfmes)

#   calculate residual curve

    npeak = 0
    curres,curve,epsres = CurveDiff(dens,dist,nfmes,bpeak,cpeak,rpeak,npeak)

#--------------------------------------
#
#   main cycle over peaks including the residual one in the origin

    while epsres >= ceps and npeak < mxp:

#       find next group of peaks

        bpeak,cpeak,rpeak,npeak = MorePeaks(curres,dist,nfmes,ceps,
                                            bpeak,cpeak,rpeak,npeak,epsp,bmin,cmin,rmin,mxp)

#       refine estimated parameters

        bpeak,cpeak,rpeak = RefineBCR(dens,dist,bpeak,cpeak,rpeak,npeak,bmin,cmin,rmin,nfmes)

#       remove the contribution of the modeled peaks

        curres,curve,epsres = CurveDiff(dens,dist,nfmes,bpeak,cpeak,rpeak,npeak)

#   output final coefficients

    OutCoefficients(bpeak,cpeak,rpeak,npeak,nfmes)

    return bpeak,cpeak,rpeak,npeak,curve,curres
