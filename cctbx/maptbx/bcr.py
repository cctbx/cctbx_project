from __future__ import absolute_import, division, print_function
import math
from scipy.optimize import minimize
from libtbx import adopt_init_args

from cctbx.array_family import flex
from mmtbx.ncs import tncs

# Adopted by Pavel Afonine. This is a verbatim copy of the original code
# supplied by A. Urzhumtsev on 11/5/21, 08:19.

#    to call the function :
#
#    bpeak,cpeak,rpeak,npeak,curve,curres = get_BCR(dens,dist,ndist,mxp,epsc,kpres,nfmes)
#
#    here input parameters
#
#    dens, dist - two arrays with the distance and the function values
#                 (I think the grid step may be NOT constant but I did not test this)
#    ndist      - maximal index of these arrays
#                 (minimal index equal to 0 and dens[0] should be also equal to 0)
#                 if you need to model the curve up to a shorter limit, assign smaller ndist value
#    mxp        - maximal allowed total number of terms (Gaussians + chi) in the sums
#    epsc       - accuracy (e.g. 0.001)
#    kpres      - flag (= 1 means accuracy in the absolute scale;
#                       = 0 means relative accuracy, with respect to the dens[0] value)
#    nfmes      - message file ; if filemes is the filename, define the file by the operator
#                 nfmes = open(filemes, 'w')
#
#    output parameters
#
#    bpeak,cpeak,rpeak - arrays of coefficiens; Gaussians are indicated by rpeak[i]=0 ;
#                        for example the first set is always the Gaussian, rpeak[0] = 0,
#                        this is the principal contribution for the peak in the origin.
#    npeak             - the total number of triplets of coefficients (the terms in the sums)
#    curve             - contributioj of each of terms in the sum ;
#                        array [ [ 0 for j in range(mxp+1)] for i in range(ndist+1) ]
#    curres            - residual curve ; array [ [0] for i in range(ndist+1) ]
#
######################################################
#
#    decomposition of oscillating curves (e.g. atomic images at a given resolution)
#    by A.Urzhumtsev, L.Urzhumtseva, 2021
#
######################################################

class calculator(object):

  def __init__(self, npeak,dens,yc,dist,ndist,nfmes, x):
    adopt_init_args(self, locals())
    self.x = flex.double(x)

  def target_and_gradients(self, x):
    self.x = x
    target = FuncFit(self.x, self.npeak, self.dens, self.yc, self.dist,
      self.ndist, self.nfmes)
    gradients = GradFit(self.x, self.npeak, self.dens, self.yc, self.dist,
      self.ndist, self.nfmes)
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


def PreciseData(kpres,epsc,peak,dstep,nfmes):

#   ceps - defines accuracy in absolute values
#   bmin - defines sharpest drop < 10**(-ndrop) of the 3D-exponential
#            at a distance = grid step
#   cmin - defines minimal peak height (64 majorates (4*pi)**1.5 )
#   epsp - drop when approximating internal peaks

    epsp = 0.1

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
def PeakOrigin(dens,dist,ndist,epsp,bmin,cmin,nfmes):

#   invert the curve if negative peak

    d0 = dens[0]
    if d0 > 0.0:
      kneg=0
    else:
      kneg=1
      dens,d0 = CurveInvert(dens,ndist,d0)

#   parameters of the Gaussian peak in the origin

    b0,c0 = GaussBC(dist,ndist,epsp,bmin,cmin,dens,nfmes)

#   invert curve back if necessary

    if kneg >0 :
      dens,d0 = CurveInvert(dens,ndist,d0)

    return b0,c0

#============================
def CurveInvert(dens,ndist,d0):

    for i in range(ndist+1):
      dens[i] = -dens[i]

    d0 = -d0

    return dens,d0

#============================
def GaussBC(dist,ndist,epsp,bmin,cmin,curres,nfmes):

    cpi   = math.pi
    pisq4 = 4.*cpi*cpi
#
#   inflection point and approximate curve till this point
#
    kpeak = 0
    minf = NumInflRight(curres,ndist,kpeak,epsp,nfmes)
    if minf == 0:
       if(nfmes is not None):
         print(' too sharp peak in the origin',curres[0],curres[1],file=nfmes)
       exit()

    ug,vg = uvFitGauss(dist,curres,minf)
    b0 = pisq4/vg
    c0 = ug*(cpi/vg)**1.5

    if b0 < bmin:
       b0 = bmin
    if c0 < cmin:
       c0 = cmin

    return b0,c0

#============================
def RippleBCR(dist,ndist,kpeak,epsp,bmin,cmin,rmin,curres,nfmes):

#   inflection point and approximate curve till this point

    minfr = NumInflRight(curres,ndist,kpeak,epsp,nfmes)
    minfl = NumInflLeft(curres,ndist,kpeak,epsp,nfmes)
    if minfr-kpeak > kpeak-minfl:
      kref = minfr
    else:
      kref = minfl

    if kref == kpeak:
      if(nfmes is not None):
        print(' too sharp peak ',kpeak,curres[kpeak-1],curres[kpeak],
          curres[kpeak+1],file=nfmes)
      exit()

    b0,c0,r0 = FitBCR(dist,ndist,kpeak,kref,bmin,cmin,rmin,curres,nfmes)

    return b0,c0,r0

#============================
def TailBCR(dist,ndist,epsp,bmin,cmin,rmin,curres,nfmes):

#   inflection point for the peak at the upper interval bound

    kref = NumInflLeft(curres,ndist,ndist,epsp,nfmes)

    if kref == 0:
      if(nfmes is not None):
        print(' no inflection point found for the tail peak',file=nfmes)
      exit()

#   parameters for this peak

    b0,c0,r0 = FitBCR(dist,ndist,ndist,kref,bmin,cmin,rmin,curres,nfmes)

    return b0,c0,r0

#============================
def FitBCR(dist,ndist,kpeak,kref,bmin,cmin,rmin,curres,nfmes):

#   BCR parameters for an internal ripple

    pi4   = 4.  * math.pi
    pisq4 = pi4 * math.pi
    epsp  = 0.1

#   scorr is artificial correcting factor against B overestimation
#   when r0 is taken without correction B/(8*peak*pi**2)

    scorr=1.5

#   check if the next curve values is negative

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
def NumInflRight(curres,ndist,kpeak,epsp,nfmes):

#   right inflrction point

    clim = epsp * curres[kpeak]
    if curres[kpeak+1] <= clim :
       i = kpeak
    else:
       for i in range (kpeak,ndist):
           if curres[i+1] <= clim or curres[i+1]+curres[i-1]-2.*curres[i] >= 0.0 :
              break

    return i

#============================
def NumInflLeft(curres,ndist,kpeak,epsp,nfmes):

#   left inflection point

    clim = epsp * curres[kpeak]

    if curres[kpeak-1] <= clim :
       i = kpeak
    else:
       for i in range (kpeak-1,1,-1):
           if curres[i-1] <= clim or curres[i+1]+curres[i-1]-2.*curres[i] >= 0.0 :
              break

    return i

#============================
def uvFitGauss(dist,curres,minf):

#   used below : w = ln(u)

    sumr2, sumr4, sumdl, sumrdl = 0.0 , 0.0 , 0.0 , 0.0

    for i in range(minf):
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
def CurveDiff(dens,dist,ndist,nfmes,bpeak,cpeak,rpeak,npeak):

    curve  = [ [ 0 for j in range(npeak+1)] for i in range(ndist+1) ]
    curres = [ dens[i] for i in range(ndist+1) ]

#   cycle over components

    for ipeak in range(npeak+1):

#     select the  parameters

      b0, c0, r0 = bpeak[ipeak] , cpeak[ipeak] , rpeak[ipeak]

#     calculate the contribution
#         ripple in point r0 or Gaussian in the origin

      if r0 > 0.0 :
        tcurve = CurveRipple(dist,ndist,b0,c0,r0)
      else:
        tcurve = CurveGauss(dist,ndist,b0,c0)

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
      print('with {npeak+1:4} peaks modeled max residual peak is {epsres:12.7f}',
            file=nfmes)

    return curres,curve,epsres

#============================
def CurveGauss(dist,ndist,b0,c0):

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
def CurveRipple(dist,ndist,b0,c0,r0):

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
def OutCurves(dens,dist,curres,curve,ndist,filres):

    nfcurv = open(filres, 'w')

    print('   dist   curve   modeled   resid',file=nfcurv)
    for i in range(ndist):
        densum = dens[i]-curres[i]
        print('{dist[i]:7.3f}{dens[i]:9.5f}{densum:9.5f}{curres[i]:9.5f}',
            file=nfcurv)

    return

#============================
def OutCoefficients(bpeak,cpeak,rpeak,npeak,nfmes):
    if(nfmes is None): return
    space=' '

    print(' final parameter values:',file=nfmes)
    print(' Npeak     Bpeak',space*13,'Cpeak',space*13,'Rpeak',file=nfmes)
    for ipeak in range(npeak+1):
       print('{ipeak:6}{bpeak[ipeak]:20.14f}{cpeak[ipeak]:20.14f}{rpeak[ipeak]:20.14f}',
             file=nfmes)

    return

#============================
def MorePeaks(curres,dist,ndist,nfmes,ceps,bpeak,cpeak,rpeak,npeak,epsp,bmin,cmin,rmin,mxp):

#   first iteration: the peak in the origin has been already treated
#   get starting point in the curve for further search

    if npeak == 0:
       kpeak=1
    else:
       kpeak=0

#   cycle over peaks

    for ipeak in range(npeak+1,mxp+1):

#       next maximal residual peak if yet significant

        kpeak,kneg = NextPeak(curres,ndist,kpeak,nfmes,ceps)

#       no more significant peak found

        if kpeak < 0:
           break

#       peak in the origin - model it by a Gaussian

        elif kpeak == 0:
           b0 , c0 = GaussBC(dist,ndist,epsp,bmin,cmin,curres,nfmes)
           r0 = 0.0

#       internal peak

        elif kpeak < ndist:
          b0, c0, r0 = RippleBCR(dist,ndist,kpeak,epsp,bmin,cmin,rmin,curres,nfmes)

#       peak at the right border of the interval

        else:
          b0, c0, r0 = TailBCR(dist,ndist,epsp,bmin,cmin,rmin,curres,nfmes)

        if kneg > 0:
          curres, c0 = CurveInvert(curres,ndist,c0)

        bpeak[ipeak] , cpeak[ipeak] , rpeak[ipeak] = b0 , c0, r0
        kpeak = kpeak+1

    if kpeak < 0:
       npeak = ipeak - 1
    else:
       npeak = ipeak

    return bpeak,cpeak,rpeak,npeak

#============================
def NextPeak(curres,ndist,kpeak,nfmes,ceps):

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
       curres,cx = CurveInvert(curres,ndist,cx)

    return kpeak,kneg

#============================
def RefineBCR(dens,dist,ndist,bpeak,cpeak,rpeak,npeak,bmin,cmin,rmin,nfmes):

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
           bcrbounds[i3+1] = (-1000,-cmin)
        bcrbounds[i3+2] = (rmin,1000)

#   minimization

    for it in range(1,3):
      if it > 1:
        xc = res.x
      CALC = calculator(npeak,dens,yc,dist,ndist,nfmes, xc)
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

        if(nfmes is not None):
          print('{i:6}{bpeak[i]:16.10f}{cpeak[i]:16.10f}{rpeak[i]:10.5f}   ',
                '{b0:16.10f}{c0:16.10f}{r0:10.5f}',file=nfmes)

        bpeak[i] , cpeak[i] , rpeak[i] = b0 , c0 , r0

    return bpeak,cpeak,rpeak

#============================
def FuncFit(xc,npeak,y0,yc,dist,ndist,nfmes):

#   calculate model curve

    yc = CurveCalc(xc,npeak,yc,dist,ndist,nfmes)

#   comparer model and control curves

    fvalue = LSFunc(yc,y0,ndist)

    return fvalue

#============================
def GradFit(xc,npeak,y0,yc,dist,ndist,nfmes):

#   calculate the gradient with respect to the curve yc

    gy = GradCurve(yc,y0,ndist)

#   calculate the gradient with respect to the parameters

    gx = GradX(xc,npeak,yc,gy,dist,ndist,nfmes)

    return gx

#============================
def CurveCalc(xc,npeak,yc,dist,ndist,nfmes):

    for ix in range(ndist+1):
        yc[ix] = 0.0

    for i in range(npeak+1):
        i3 = i*3
        b0 , c0 , r0 = xc[i3] , xc[i3+1] , xc[i3+2]

#       contribution from a ripple or from a Gaussian in the origin

        if r0 > 0.0 :
           tcurve = CurveRipple(dist,ndist,b0,c0,r0)
        else:
           tcurve = CurveGauss(dist,ndist,b0,c0)

        for ix in range(ndist):
            yc[ix] = yc[ix] + tcurve[ix]

    return yc

#============================
def LSFunc(yc,y0,ndist):

    fvalue = 0.0

    for i in range(ndist+1):
        fvalue = fvalue + (yc[i]-y0[i])**2

    fvalue = fvalue * 0.5

    return fvalue

#============================
def GradCurve(yc,y0,ndist):

    gy = [ yc[i]-y0[i] for i in range(ndist+1) ]

    return gy

#============================
def GradX(xc,npeak,yc,gy,dist,ndist,nfmes):

    gx = [ 0 for i in range(3*(npeak+1)) ]

    for i in range(npeak+1):
       i3 = i * 3
       b0 , c0, r0 = xc[i3] , xc[i3+1] , xc[i3+2]

       if r0 > 0.0:
          gb, gc, gr = GradRipple(yc,gy,dist,ndist,b0,c0,r0)
       else:
          gb, gc = GradGauss(yc,gy,dist,ndist,b0,c0)
          gr    = 0.0

       gx[i3] , gx[i3+1] , gx[i3+2]  = gb , gc , gr

    return gx

#============================
def GradRipple(yc,gy,dist,ndist,b0,c0,r0):

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
def GradGauss(yc,gy,dist,ndist,b0,c0):

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

    ndist = dens.size()-1

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

    bpeak[0],cpeak[0] = PeakOrigin(dens,dist,ndist,epsp,bmin,cmin,nfmes)

#   calculate residual curve

    npeak = 0
    curres,curve,epsres = CurveDiff(dens,dist,ndist,nfmes,bpeak,cpeak,rpeak,npeak)

#--------------------------------------
#
#   main cycle over peaks including the residual one in the origin

    while epsres >= ceps and npeak < mxp:

#       find next group of peaks

        bpeak,cpeak,rpeak,npeak = MorePeaks(curres,dist,ndist,nfmes,ceps,
          bpeak,cpeak,rpeak,npeak,epsp,bmin,cmin,rmin,mxp)

#       refine estimated parameters

        bpeak,cpeak,rpeak = RefineBCR(dens,dist,ndist,bpeak,cpeak,rpeak,npeak,
          bmin,cmin,rmin,nfmes)

#       remove the contribution of the modeled peaks

        curres,curve,epsres = CurveDiff(dens,dist,ndist,nfmes,bpeak,cpeak,rpeak,
          npeak)

    return bpeak,cpeak,rpeak,npeak,curve,curres
