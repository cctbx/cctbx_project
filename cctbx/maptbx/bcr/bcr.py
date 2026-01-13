from __future__ import absolute_import, division, print_function
import math
from scipy.optimize import minimize
from libtbx import adopt_init_args
from cctbx import maptbx

from cctbx.array_family import flex
import scitbx.minimizers

import json

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

def curve(B, C, R, radii, b_iso=0):
  result = flex.double()
  # FILTER
  B_, C_, R_ = [],[],[]
  for bi, ci, ri in zip(B, C, R):
    if(abs(bi)<1.e-6 or abs(ci)<1.e-6): continue
    else:
      B_.append(bi)
      C_.append(ci)
      R_.append(ri)
  #
  for r in radii:
    first = 0
    second = 0
    for bi, ci, ri in zip(B_, C_, R_):
      if(abs(ri)<1.e-6):
        first += gauss(B=bi, C=ci, r=r, b_iso=b_iso)
      else:
        second += ci*chi(B=bi, R=ri, r=r, b_iso=b_iso)
    result.append(first + second)
  return result

class calculator(object):

  def __init__(self, npeak,dens,yc,dist,nfmes, x, mdist,edist,
               bound_flags, lower_bound, upper_bound):
    adopt_init_args(self, locals())
    self.x = flex.double(x)

  def update(self, x):
    self.x = x

  def gradients(self):
    return flex.double(GradFit(self.x, self.npeak, self.dens, self.yc, self.dist,
      self.mdist,self.edist, self.nfmes))

  def target(self):
    return FuncFit(self.x, self.npeak, self.dens, self.yc, self.dist,
      self.mdist,self.edist, self.nfmes)

#=================================================================
def get_BCR(dens,dist,dmax,mxp,epsc,epsp,edist,kpres,kprot,nfmes):


#    dens   - array of the curve values,
#             in in creasing order of the argument starting from 0
#    dist   - array of the argument values ; normally this is a regular grid
#    dmax   - maximal distance value till which the curve will be approximated
#    mxp    - maximal number of terms ; except special cases, it is recommended
#             to choose a large value (e.g., 1000) and limit the number of terms
#             rather referring to the parameter epsc below
#    epsc   - precision with which the decomposition will be constructed;
#             A KEY PARAMETER ;
#             the value recommended is 0.001; at least for atomic images, this gives
#             an accurate approximation with a single term per a local peak for most of
#             atomic images at conventional resolution
#    epsp   - part of the value of the local peak defining the peak limits;
#             use values in the ragne 0.0 - 0.005
#    edist  - value at which the atomic contribution considered as negligibly small;
#             1.0E-13 ie recommended;
#    kpres  - defines if the precision is given in absolute values (=1) or as a part of
#             peak in the origin, dens[0]
#    kprot  - type of the peak analysis processing
#             111 - single iteration; refinement at the end
#             112 - single iteration; refinement at the end; pre-refinement of the first term
#             122 - single iteration; instant refinement of each term
#             211 - multiple iterations; refinement at the end
#             212 - multiple iterations; refinement at the end; pre-refinement of the first term
#             222 - multiple iterations; instant refinement of each term at the first iteration

    if nfmes is not None:
      print('',file=nfmes)
      print(30*'=',' NEW CURVE IN PROCESS ',30*'=',file=nfmes)

      print('')
      print(30*'=',' NEW CURVE IN PROCESS ',30*'=')

    bpeak  = [ 0 for i in range(2*mxp+2) ]
    cpeak  = [ 0 for i in range(2*mxp+2) ]
    rpeak  = [ 0 for i in range(2*mxp+2) ]

#   definition / estimations of the internal parameters

    ceps,bmin,cmin,rmin,edist,mdist = PreciseData(kpres,epsc,epsp,edist,
                                                  dmax,dens[0],dist,nfmes)
#--------------------------------------
#
#   starting point to search for the peaks;
#   number of terms found previously

    kpeak = 0
    npeak = -1
    curres,curve,epsres =  CurveDiff(dens,dist,mdist,edist,nfmes,bpeak,cpeak,rpeak,npeak,mxp)

#--------------------------------------

#   main cycle over peaks including the residual one in the origin

    while (epsres >= ceps) and (npeak < mxp) :

#       find next group of peaks

        kpeak,kneg,curres = NextPeak(curres,mdist,kpeak,nfmes,ceps)

        if kpeak >= 0 and kpeak <= mdist :

           bpeak,cpeak,rpeak,npeak,kdist,kpeak = PeakBCR(curres,dist,mdist,ceps,kneg,kpeak,
                                            bpeak,cpeak,rpeak,npeak,epsp,bmin,cmin,rmin,mxp,nfmes)

           if npeak == mxp or (kprot == 122 or kprot == 222) or \
              (kprot == 112 and npeak == 0) or (kprot == 212 and npeak == 0) :

#             if the last (mxp) peak found                       or
#             if each term refined instantly                     or
#             if this is the first peak to be refined instantly
#             - refine parameters in the extended interval
#             - and remove contribution

              bpeak,cpeak,rpeak =        \
                    RefineBCR(dens,dist,mdist,edist,bpeak,cpeak,rpeak,npeak,bmin,cmin,rmin,nfmes)

              curres,curve,epsres = \
                               CurveDiff(dens,dist,mdist,edist,nfmes,bpeak,cpeak,rpeak,npeak,mxp)
        else :

#          no more peaks ; the whole interval has been analysed
#          - refine estimated parameters, all together, in the required interval
#          - remove the contribution of the modeled peaks and update start search point

           if kprot == 111 or kprot == 112 or kprot == 211 or kprot == 212 :

              bpeak,cpeak,rpeak =        \
                    RefineBCR(dens,dist,mdist,edist,bpeak,cpeak,rpeak,npeak,bmin,cmin,rmin,nfmes)

              curres,curve,epsres = \
                               CurveDiff(dens,dist,mdist,edist,nfmes,bpeak,cpeak,rpeak,npeak,mxp)

           if kprot == 111 or kprot == 112 or kprot == 122 :
              break
           elif kprot == 222 :
              kprot = 211

           kpeak = 0
#
#   end of the main cycle ;
#
#   filter out weak terms if possible

    bpeak,cpeak,rpeak,mpeak,curve,curres,epsres = FilterWeak(dens,curve,curres,dist,mdist, \
                       edist,epsres,bpeak,cpeak,rpeak,npeak,mxp,bmin,cmin,rmin,nfmes,ceps)

    return bpeak,cpeak,rpeak,mpeak,curve,curres,epsres

#============================
def RefineBCR(dens,dist,mdist,edist,bpeak,cpeak,rpeak,npeak,bmin,cmin,rmin,nfmes):

    ndist = len(dist) - 1

    xc        = [ 0 for i in range(3*(npeak+1)) ]
    yc        = [ 0 for i in range(ndist+1) ]
    bcrbounds = [ (None,None) for i in range(3*(npeak+1)) ]

    sp = ' '

#   define bounds
#   cmin is the bound for |c| and not for c itself

#!!!!!!!!!!!!!! inserted by P.Afonine instead of the original call of "minimize"

    for i in range(npeak+1):
        i3 = i*3
        xc[i3]   = bpeak[i]
        xc[i3+1] = cpeak[i]
        xc[i3+2] = rpeak[i]
        bcrbounds[i3]   = (bmin,1000.)
        bcrbounds[i3+2] = (rmin,dist[ndist]*1.2)
        if cpeak[i] > 0.0 :
           bcrbounds[i3+1] = (cmin,1000.)
        else :
           bcrbounds[i3+1] = (-1000.,-cmin)

#   minimization
    if nfmes is not None:
      print('')
      print(' L-BFGS-B refinement in progress ; parameters before refinement : ')
    for i in range(npeak+1):
       if nfmes is not None:
         print(f'{i+1:6}{rpeak[i]:16.10f}{bpeak[i]:16.10f}{cpeak[i]:11.5f}')

    if 0: # lbfgsb from cctbx
      #for it in range(1,10):
      #for it in range(1,2):
      for it in [0,]:
        lbound = []
        ubound = []
        for b in bcrbounds:
          lbound.append(b[0])
          ubound.append(b[1])
        CALC = calculator(npeak,dens,yc,dist,nfmes, xc, mdist,edist,
          bound_flags = flex.int(len(xc), 2),
          lower_bound = lbound,
          upper_bound = ubound)
        res = scitbx.minimizers.lbfgs(
             mode='lbfgsb', max_iterations=500, calculator=CALC)
        res.x = list(res.x)
        xc = res.x
        if 1: # mix
          res = minimize(FuncFit,xc,args=(npeak,dens,yc,dist,ndist,edist,nfmes),
               method='L-BFGS-B',jac=GradFit, bounds = bcrbounds)
          xc = res.x

    else: # SciPy analogue. Works with Python 3 only.
      res = minimize(FuncFit,xc,args=(npeak,dens,yc,dist,ndist,edist,nfmes),
               method='L-BFGS-B',jac=GradFit, bounds = bcrbounds)

#   recover refined values
    if nfmes is not None:
      print('',file=nfmes)
      print(' parameters before and after their refinement:'     ,file=nfmes)
      print(' Nterm    R (M)   ',6*' ','B (N)   ',6*' ','C (K)   ',5*' ',
                     'R (M)   ',6*' ','B (N)   ',6*' ','C (K)   ',file=nfmes)

      print('')
      print(' parameters before and after their refinement:')
      print(' Nterm    R (M)   ',6*' ','B (N)   ',6*' ','C (K)   ',5*' ',
                     'R (M)   ',6*' ','B (N)   ',6*' ','C (K)   ')

    for i in range(npeak+1):
        i3 = i*3
        b0 , c0, r0 = res.x[i3] , res.x[i3+1] , res.x[i3+2]
        if nfmes is not None:
          print(f'{i+1:6}{rpeak[i]:16.10f}{bpeak[i]:16.10f}{cpeak[i]:11.5f}   ',
                f'{r0:16.10f}{b0:16.10f}{c0:11.5f}',file=nfmes)

          print(f'{i+1:6}{rpeak[i]:16.10f}{bpeak[i]:16.10f}{cpeak[i]:11.5f}   ',
                f'{r0:16.10f}{b0:16.10f}{c0:11.5f}')

        bpeak[i] , cpeak[i] , rpeak[i] = b0 , c0 , r0

    return bpeak,cpeak,rpeak

#============================
def PreciseData(kpres,epsc,epsp,edist,dmax,peak,dist,nfmes):

#   ceps  - defines accuracy in absolute values
#   bmin  - defines sharpest drop < 10**(-ndrop) of the 3D-exponential
#            at a distance = grid step
#   cmin  - defines minimal peak height (64 majorates (4*pi)**1.5 )
#   epsp  - drop when approximating internal peaks
#   edist - limit of the term contribution

    ndist = len(dist) -1

#   precision constants

    if kpres > 0:
      ceps  = epsc
    else:
      ceps  = abs(peak) * epsc
      edist = abs(peak) * edist

#   find the point in 'dist' for the required distance limit

    mdist = ndist
    if dist[ndist] < dmax :
       if nfmes is not None:
         print('*** warning : curve interval is defined for x <= ',dist[ndist])
         print('              shorter that the required distance  ',dmax)
    else :
       for ix in range(ndist+1) :
           if dist[ix] > dmax :
              mdist = ix - 1
              break

#   estimates for a regular grid with step = dist[1] - dist[0]

    dstep = dist[1]
    ndrop = 2.
    bmin  = 16.*dstep*dstep/ndrop
    cmin  = ceps*(bmin**1.5)/64.
    rmin  = 0.0

    if nfmes is not None:
      print ('',file=nfmes)
      print (' INTERNAL DECOMPOSITION PARAMETERS ',file=nfmes)
      print ('',file=nfmes)
      print (' absolute max allowed error     ',f'{ceps:13.5e}' ,file=nfmes)
      print (' drop to estimate the peak width',f'{epsp:13.5e}' ,file=nfmes)
      if edist > 0.0 :
         print (' term extension limit           ',f'{edist:13.5e}',file=nfmes)
      else :
         print (' term extension limit :            no limit applied',file=nfmes)
      print (' estimated accuracy parameters :',file=nfmes)
      print ('    bmin                        ',f'{bmin:13.5e}' ,file=nfmes)
      print ('    cmin                        ',f'{cmin:13.5e}' ,file=nfmes)
      print ('---------------------------',file=nfmes)

    return ceps,bmin,cmin,rmin,edist,mdist

#============================
def CurveInvert(dens,d0):

    ndist = len(dens) - 1

    for i in range(ndist+1):
      dens[i] = -dens[i]

    d0 = -d0

    return dens,d0

#============================
def GaussBC(curres,dist,epsp,bmin,cmin,nfmes):

    ndist = len(dist) -1

    cpi   = math.pi

#   inflection point and approximate curve till this point

    kpeak = 0
    r0    = dist[kpeak]

    minf = NumInflRight(curres,kpeak,epsp)
    if minf == 0:
       b0 = 2.0*bmin
       c0 = curres[kpeak] * (b0/(4.0*cpi))**1.5
    else:
       ug,vg = uvFitGaussR(curres,dist,r0,kpeak,minf)
       b0 = 4.0 * cpi**2 / vg
       c0 = ug  * (cpi/vg)**1.5
       if b0 < 2.0*bmin:
          b0 = 2.0*bmin

    if c0 < 2.0*cmin:
       c0 = 2.0*cmin

    return b0,c0,r0,minf

#============================
def RippleBCR(curres,dist,mdist,kpeak,epsp,bmin,cmin,rmin,nfmes):

#   inflection point and approximate curve till this point

    ndist = len(dist) -1
    cpi   = math.pi

    minfr = NumInflRight(curres,kpeak,epsp)
    minfl = NumInflLeft(curres,kpeak,epsp)

    if minfr-minfl <= 1:
       b0 = 2.0*bmin
       c0 = curres[kpeak] * math.sqrt(4.0*cpi*b0) * dist[kpeak]**2
    else:
       ug,vg = uvFitGaussR(curres,dist,dist[kpeak],minfl,minfr)
       b0    = 4.0 * cpi * cpi / vg
       if b0 < 2.0*bmin:
          b0 = 2.0*bmin
       c0    = ug * math.sqrt(b0)

    if c0 < 2.0*cmin:
       c0 = 2.0*cmin
    r0 = dist[kpeak]

    kdist = minfr

    return b0,c0,r0,kdist

#============================
def TailBCR(curres,dist,mdist,epsp,bmin,cmin,rmin,nfmes):

#   inflection point for the peak at the upper interval bound

    ndist = len(dist) -1

#   the peak at the upper bound of the required interval

    kpeak = mdist
    if mdist < ndist :
       kpeak = NumPeakRight(curres,mdist)

    if kpeak < ndist :

#      an internal peak in the extended interval

       b0, c0, r0, kdist = RippleBCR(curres,dist,mdist,kpeak,epsp,bmin,cmin,rmin,nfmes)

    else :

#      otherwise find the inflection point to the left
#      and approximate the curve till this point

       b0, c0, r0, kdist = border(dist,mdist,epsp,bmin,cmin,rmin,curres,nfmes)

    return b0, c0, r0, kdist

#============================
def border(dist,kpeak,epsp,bmin,cmin,rmin,curres,nfmes):

    ndist = len(dist) - 1

    kdist = ndist
    r0    = dist[ndist]
    kref = NumInflLeft(curres,ndist,epsp)

    if kref == ndist :
          b0 = 2.0*bmin
          c0 = curres[ndist] * math.sqrt(4*math.pi*b0) * r0**2
    elif kref == ndist-1 :
          b0 = 4.0*bmin
          c0 = curres[ndist] * math.sqrt(4*math.pi*b0) * r0**2
    else :
       ug,vg = uvFitGaussR(curres,dist,r0,kref,ndist)
       b0    = 4.0 * math.pi**2 / vg
       if b0 < 2.0*bmin:
          b0 = 2.0*bmin
       c0    = ug * math.sqrt(b0)

    if c0 < 2.0*cmin:
       c0 = 2.0*cmin

    return b0,c0,r0,kdist

#============================
def NumPeakRight(curres,npeak):

#   peak in the interval npeak <= ndist
#   (if nreak is on the upper bound of the search interval)

    ndist = len(curres) - 1

    npeakr = npeak
    if npeak < ndist :
       for i in range (npeak+1,ndist+1) :
           if curres[i] < curres[i-1] :
              npeakr = i - 1
              break

    return npeakr

#============================
def NumInflRight(curres,kpeak,epsp):

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
def NumInflLeft(curres,kpeak,epsp):

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
def uvFitGaussR(curres,dist,d0,minx,maxx):

#   used below : w = ln(u)

    sumr2, sumr4, sumdl, sumrdl = 0.0 , 0.0 , 0.0 , 0.0
    sumx                        = maxx - minx + 1

    if d0 > 0.0 :
       scalf = d0 * d0 * math.sqrt(4.0 * math.pi)
    else :
       scalf = 1.0

    for i in range(minx,maxx+1):
      curvei = math.log(curres[i] * scalf)
      dr = dist[i] - d0
      ri2    = dr * dr
      ri4    = ri2    * ri2
      sumr2  = sumr2  + ri2
      sumr4  = sumr4  + ri4
      sumdl  = sumdl  + curvei
      sumrdl = sumrdl + curvei * ri2

    det  =  -sumx * sumr4  + sumr2 * sumr2
    detw = -sumdl * sumr4  + sumrdl*sumr2
    detv =   sumx * sumrdl - sumr2 * sumdl

    wg = detw / det
    vg = detv / det
    ug = math.exp(wg)

    return ug,vg

#============================
def CurveDiff(dens,dist,mdist,edist,nfmes,bpeak,cpeak,rpeak,npeak,mxp):

    ndist = len(dist) - 1

    curve  = [ [ 0 for j in range(mxp+1)] for i in range(ndist+1) ]
    curres = [ dens[i] for i in range(ndist+1) ]

#   cycle over components

    for ipeak in range(npeak+1):

#     select the  parameters

      b0, c0, r0 = bpeak[ipeak] , cpeak[ipeak] , rpeak[ipeak]

#     calculate the contribution
#         ripple in point r0 or Gaussian in the origin

      if r0 > 0.0 :
        tcurve = CurveRipple(dist,edist,b0,c0,r0)
      else:
        tcurve = CurveGauss(dist,edist,b0,c0)

#     remove contribution

      for ix in range(ndist+1):
        curve [ix][ipeak] = tcurve[ix]
        curres[ix]        = curres[ix] - tcurve[ix]

#   end of cycle over components; residual errors

    epsres = 0.0
    for ix in range (mdist+1):
      curabs = abs(curres[ix])
      if curabs > epsres:
         epsres = curabs

    totres = epsres
    if ndist > mdist :
       for ix in range (mdist+1,ndist+1):
         curabs = abs(curres[ix])
         if curabs > totres:
            totres = curabs

    if nfmes is not None:
      space = 23*' '
      print('',file=nfmes)
      print(f' with {npeak+1:4} terms max residual peaks are',file=nfmes)
      print(space,' inside the interval of modeling',f'{epsres:12.7f}',file=nfmes)
      print(space,' in the whole range of distances',f'{totres:12.7f}',file=nfmes)

      print('')
      print(f' with {npeak+1:4} terms max residual peaks are')
      print(space,' inside the interval of modeling',f'{epsres:12.7f}')
      print(space,' in the whole range of distances',f'{totres:12.7f}')

    return curres,curve,epsres

#============================
def CurveGauss(dist,edist,b0,c0):

    ndist = len(dist) -1

#   curve for a Gaussian in the origin

    curve = [ 0.0 for i in range(ndist+1) ]

    cpi  = math.pi
    pi4b = 4. * cpi / b0
    vg   = pi4b * cpi
    ug   = c0 * pi4b**1.5

    if edist > 0.0 :
       distmx = math.sqrt(-math.log(edist/abs(ug)) / vg)
    else :
       distmx = dist[ndist]

    for ix in range(ndist+1):
      distx = dist[ix]
      if distx <= distmx :
         arg = vg * distx**2
         curve[ix] = ug * math.exp(-arg)
      else :
         break

    return curve

#============================
def CurveRipple(dist,edist,b0,c0,r0):

    ndist = len(dist) -1

#   curve for a rippe in point r0

    curve = [ 0.0 for i in range(ndist+1) ]

    cpi   = math.pi
    pi4   = 4.0 * cpi
    pisq4 = pi4 * cpi
    ur    = c0  * (pi4/b0)**1.5
    vr    = pisq4 / b0
    vr2   = 2.0 * vr

    scaler = c0 / (r0*math.sqrt(pi4*b0))
    if edist > 0.0 :
       arg0mx = -math.log(edist/abs(ur))
    else :
       arg0mx = vr * dist[ndist]**2

    arg0 = vr * r0**2
    if arg0 < arg0mx :
       curve[0] = c0 * math.exp(-arg0) * (pi4/b0)**1.5

#   find the distance point for r0

    ir = 0
    for ix in range(1,ndist+1) :
        if dist[ix] > r0 :
           ir = ix
           break
    if ir == 0 :
       ir = ndist
    else :
       ir = ir - 1

#   left part of the curve

    if ir > 0 :
       for ix in range(ir,0,-1) :
           rx         = dist[ix]
           argm, argp = vr * (rx-r0)**2 , vr * (rx+r0)**2
           curvex     = scaler * (math.exp(-argm) - math.exp(-argp)) / rx
           if abs(curvex) >= edist :
              curve[ix]   = curvex
           else :
              break

#   right part of the curve

    if ir < ndist :
       for ix in range(ir+1,ndist+1) :
           rx         = dist[ix]
           argm, argp = vr * (rx-r0)**2 , vr * (rx+r0)**2
           curvex     = scaler * (math.exp(-argm) - math.exp(-argp)) / rx
           if abs(curvex) >= edist :
              curve[ix]   = curvex
           else :
              break

    return curve

#============================
def PeakBCR(curres,dist,mdist,ceps,kneg,kpeak,bpeak,cpeak,rpeak,npeak,       \
              epsp,bmin,cmin,rmin,mxp,nfmes):

    ndist = len(dist) - 1
    kdist = mdist

    if kpeak == 0:

#      peak in the origin - model it by a Gaussian

       b0, c0, r0, minf = GaussBC(curres,dist,epsp,bmin,cmin,nfmes)
       kpeak = minf

    elif kpeak < mdist:

#      internal peak

       b0, c0, r0, kdist = \
                   RippleBCR(curres,dist,mdist,kpeak,epsp,bmin,cmin,rmin,nfmes)
       kpeak = kdist

#      peak at the right border of the interval

    else:
       b0, c0, r0, kdist = TailBCR(curres,dist,mdist,epsp,bmin,cmin,rmin,nfmes)
       kpeak = mdist

    if kneg > 0:
       curres, c0 = CurveInvert(curres,c0)

    kpeak = kpeak+1
    npeak = npeak+1
    bpeak[npeak] , cpeak[npeak] , rpeak[npeak] = b0 , c0, r0

    return bpeak,cpeak,rpeak,npeak,kdist,kpeak

#============================
def NextPeak(curres,mdist,kpeak,nfmes,ceps):

    kfound = 0
    kneg   = 0

#   check the peak in the origin

    if kpeak == 0:

#     check if the peak in the origin if flat

      cx = curres[0]
      i  = 1
      while curres[i]-cx == 0.0:
         i = i+1

      if cx*(curres[i]-cx) < 0.0 and abs(cx) > ceps:
         ix     = 0
         kfound = 1

#     there is no peak in the origin; check next points

      else:
         kpeak = 1

#   (we cannot use 'else' below since kpeak may change after its previous check)

    if kfound == 0 and kpeak < mdist:

#        search for an internal peak

         for ix in range(kpeak,mdist):
             cx   = curres[ix]
             delm = cx-curres[ix-1]
             delp = curres[ix+1]-cx
             if delm*delp <= 0. and delp*cx < 0. and abs(cx) > ceps:
                kfound = 2
                break

#        no internal peak found ; check the upper interval bound

         if kfound == 0 :
            cx = curres[mdist]
            if (cx-curres[mdist-1])*cx > 0 and abs(cx) > ceps :
               ix = mdist
               kfound = 3

#   flip the curve temporarily if the peak is negative

    if kfound == 0 :
       kpeak = -1
    else :
       kpeak = ix
       if cx < 0.0 :
          kneg = 1
          curres,cx = CurveInvert(curres,cx)

    return kpeak,kneg,curres

#============================
def FuncFit(xc,npeak,y0,yc,dist,mdist,edist,nfmes):

#   calculate model curve

    yc = CurveCalc(xc,npeak,yc,dist,mdist,edist,nfmes)

#   comparer model and control curves

    fvalue = LSFunc(yc,y0,mdist)

    return fvalue

#============================
def GradFit(xc,npeak,y0,yc,dist,mdist,edist,nfmes):

#   calculate the gradient with respect to the curve yc

    gy = GradCurve(yc,y0,mdist)

#   calculate the gradient with respect to the parameters

    gx = GradX(xc,npeak,yc,gy,dist,mdist,edist,nfmes)

    return gx

#============================
def CurveCalc(xc,npeak,yc,dist,mdist,edist,nfmes):

    ndist = mdist

    for ix in range(ndist+1):
        yc[ix] = 0.0

    for i in range(npeak+1):
        i3 = i*3
        b0 , c0 , r0 = xc[i3] , xc[i3+1] , xc[i3+2]

#       contribution from a ripple or from a Gaussian in the origin

        if r0 > 0.0 :
           tcurve = CurveRipple(dist,edist,b0,c0,r0)
        else:
           tcurve = CurveGauss(dist,edist,b0,c0)

        for ix in range(ndist+1):
            yc[ix] = yc[ix] + tcurve[ix]

    return yc

#============================
def LSFunc(yc,y0,mdist):

    fvalue = 0.0
    ndist = mdist

    for i in range(ndist+1):
        fvalue = fvalue + (yc[i]-y0[i])**2

    fvalue = fvalue * 0.5

    return fvalue

#============================
def GradCurve(yc,y0,mdist):

    ndist = mdist

    gy = [ yc[i]-y0[i] for i in range(ndist+1) ]

    return gy

#============================
def GradX(xc,npeak,yc,gy,dist,mdist,edist,nfmes):

    gx = [ 0 for i in range(3*(npeak+1)) ]

    for i in range(npeak+1):
       i3 = i * 3
       b0 , c0, r0 = xc[i3] , xc[i3+1] , xc[i3+2]

       if r0 > 0.0:
          gb, gc, gr = GradRipple(yc,gy,dist,mdist,edist,b0,c0,r0)
       else:
          gb, gc = GradGauss(yc,gy,dist,mdist,edist,b0,c0)
          gr    = 0.0

       gx[i3] , gx[i3+1] , gx[i3+2]  = gb , gc , gr

    return gx

#============================
def GradRipple(yc,gy,dist,mdist,edist,b0,c0,r0):

    ndist = mdist

    cpi    = math.pi
    pi4    = 4.0 * cpi
    pisq4  = pi4 * cpi
    vr     = pisq4 / b0
    vr2    = 2.0 * vr
    pisq8r = vr2 * r0

    scaler = c0 / (r0 * math.sqrt(pi4*b0))
    ur     = c0 * (pi4/b0)**1.5

    if edist > 0.0 :
       arg0mx = -math.log(edist/abs(ur))
    else :
       arg0mx = vr * dist[ndist]**2

    scc   = 1. / (2.0*r0*math.sqrt(b0*cpi))
    c0r0 , c0b0 = c0/r0 , c0/b0
    c0b02       = c0b0 / 2.0

#   gradient in the origin, x = 0

    arg0 = vr * r0**2
    if arg0 < arg0mx :
       gadd          = 8.0 * math.exp(-arg0) * gy[0] * (cpi/b0)**1.5
       gb , gc , gr  =  gadd * c0b0 * (arg0 - 1.5) , gadd , -gadd * c0 * pisq8r
    else:
       gb , gc , gr  = 0.0 , 0.0 , 0.0

#   find the distance point for r0

    ir = 0
    for ix in range(1,ndist+1) :
        if dist[ix] > r0 :
           ir = ix
           break
    if ir == 0 :
       ir = ndist
    else :
       ir = ir - 1

#   left part of the curve

    if ir > 0 :
       for ix in range(ir,0,-1) :
           rx           = dist[ix]
           xm   , xp    = rx - r0         , rx + r0
           argm , argp  = vr * xm*xm      , vr * xp*xp
           gaddm, gaddp = math.exp(-argm) , math.exp(-argp)
           delmp        = gaddm-gaddp
           curvex       = scaler * (math.exp(-argm) - math.exp(-argp)) / rx
           if abs(curvex) >= edist :
              sccx = scc * gy[ix] / rx
              gc   = gc + sccx * delmp
              gb   = gb + sccx * c0b02 * ((gaddm*(2.0*argm - 1.0) - gaddp*(2.0*argp - 1.0)))
##              gb   = gb + sccx * c0b02 * (2.0 * (gaddm*argm - gaddp*argp) - delmp)
              gr   = gr + sccx * c0r0 * (gaddm*(pisq8r*xm-1.0) + gaddp*(pisq8r*xp+1.0))
##              gr   = gr + sccx * c0r0 * (pisq8r * (gaddm*xm + gaddp*xp) - delmp)
           else :
              break

#   right part of the curve

    if ir < ndist :
       for ix in range(ir+1,ndist+1) :
           rx           = dist[ix]
           xm   , xp    = rx - r0         , rx + r0
           argm , argp  = vr * xm*xm      , vr * xp*xp
           gaddm, gaddp = math.exp(-argm) , math.exp(-argp)
           delmp        = gaddm - gaddp
           curvex       = scaler * delmp / rx
           if abs(curvex) >= edist :
              sccx = scc * gy[ix] / rx
              gc   = gc + sccx * delmp
              gb   = gb + sccx * c0b02 * ((gaddm*(2.0*argm - 1.0) - gaddp*(2.0*argp - 1.0)))
##              gb   = gb + sccx * c0b02 * (2.0 * (gaddm*argm - gaddp*argp) - delmp)
              gr   = gr + sccx * c0r0 * (gaddm*(pisq8r*xm-1.0) + gaddp*(pisq8r*xp+1.0))
##              gr   = gr + sccx * c0r0 * (pisq8r * (gaddm*xm + gaddp*xp) - delmp)
           else :
              break

    return gb,gc,gr

#============================
def GradGauss(yc,gy,dist,mdist,edist,b0,c0):

    ndist = mdist

    cpi  = math.pi
    pi4b = 4.0  * cpi / b0
    vg   = pi4b * cpi
    ugc  = pi4b**1.5
    ug   = ugc * c0
    ugb  = ug  / b0

    gb , gc = 0.0 , 0.0

    if edist > 0.0 :
       distmx = math.sqrt(-math.log(edist/abs(ug)) / vg)
    else :
       distmx = dist[ndist]

    for ix in range(ndist+1):
        distx = dist[ix]
        if distx <= distmx :
           arg = vg * distx**2
           gadd = math.exp(-arg) * gy[ix]
           gb   = gb + ugb*gadd*(arg-1.5)
           gc   = gc + ugc*gadd
        else :
           break

    return gb,gc

#============================
def FilterWeak(dens,curve,curres,dist,mdist,edist,epsres,bpeak,cpeak,rpeak,npeak,mxp,  \
                                 bmin,cmin,rmin,nfmes,ceps):

    bpeakc,cpeakc,rpeakc = CheckLow(bpeak,cpeak,rpeak,npeak,ceps,nfmes)

    mpeak = len(bpeakc) - 1

    if mpeak < npeak :

#      small contributions found; check the result after removing these terms

       bpeakc,cpeakc,rpeakc = \
         RefineBCR(dens,dist,mdist,edist,bpeakc,cpeakc,rpeakc,mpeak,bmin,cmin,rmin,nfmes)

       curres,curve,epsres = \
          CurveDiff(dens,dist,mdist,edist,nfmes,bpeakc,cpeakc,rpeakc,mpeak,mxp)

#      removing small contributions makes the residual error above the limit;
#      include these terms back

       if epsres > ceps :

          mpeak  = npeak
          bpeakc = [ bpeak[ip] for ip in range(npeak+1)]
          cpeakc = [ cpeak[ip] for ip in range(npeak+1)]
          rpeakc = [ rpeak[ip] for ip in range(npeak+1)]

          curres,curve,epsres = \
            CurveDiff(dens,dist,mdist,edist,nfmes,bpeakc,cpeakc,rpeakc,mpeak,mxp)

    return bpeakc,cpeakc,rpeakc,mpeak,curve,curres,epsres

#================================
def CheckLow(bpeak,cpeak,rpeak,npeak,ceps,nfmes):

    pi4    = 4.0 * math.pi
    pisq16 = pi4 * pi4

    bpeakc, cpeakc, rpeakc = [], [], []

    if nfmes is not None:
      print('',file=nfmes)
      print(' Check term values in the points r = r0 :',file=nfmes)
      print(' Nterm    R (M)   ',6*' ','B (N)   ',6*' ','C (K)   ',7*' ','value',file=nfmes)

      print('')
      print(' Check term values in the points r = r0 :')
      print(' Nterm    R (M)   ',6*' ','B (N)   ',6*' ','C (K)   ',7*' ','value')

    for ip in range (npeak+1):
        b0 = bpeak[ip]
        c0 = cpeak[ip]
        r0 = rpeak[ip]
        if r0 == 0.0 :
           fval = c0 * (pi4/b0)**1.5
        else :
           fval = c0 * (1.0 - math.exp(-pisq16*r0**2/b0)) / (math.sqrt(pi4*b0) * r0**2)

        if nfmes is not None:
          print(f'{ip+1:6}{r0:16.10f}{b0:16.10f}{c0:11.5f}   {fval:12.7f}',file=nfmes)
          print(f'{ip+1:6}{r0:16.10f}{b0:16.10f}{c0:11.5f}   {fval:12.7f}')

        if abs(fval) >= ceps :
           bpeakc.append(b0)
           cpeakc.append(c0)
           rpeakc.append(r0)

    return bpeakc, cpeakc, rpeakc

#=== INDEPENDNT COMPUTATIONAL PART STOPS HERE ===

#==================================================================================

def InputScat(filedata, Badd):

    NGauss = 6

    nfdat = open(filedata, 'r')
    Coefs = nfdat.readlines()
    NTypes = int(len(Coefs) / 3)

    TypesAtoms = ['' for itype in range(NTypes)]
    ScatAtoms  = [[0.0 for ig in range(2 * NGauss)] for itype in range(NTypes)]

    Iblock  = 0
    NGauss1 = NGauss - 1
    for itype in range(NTypes) :
       Line1 = Coefs[itype*3]
       Line2 = Line1.replace('{', ' ')
       Line3 = Line2.replace('}', ' ')
       Line4 = Line3.replace(',', ' ')
       Line5 = Line4.replace('"', ' ')
       Line6 = Line5.split()
       AtomType = Line6[0]
       for ig in range(NGauss1) : ScatAtoms[itype][ig] = float(Line6[ig+1])

       Line1 = Coefs[itype*3 + 1]
       Line2 = Line1.replace('{', ' ')
       Line3 = Line2.replace('}', ' ')
       Line4 = Line3.replace(',', ' ')
       Line5 = Line4.replace('"', ' ')
       Line6 = Line5.split()
       for ig in range(NGauss1) : ScatAtoms[itype][ig + NGauss] = float(Line6[ig]) + Badd

       Line1 = Coefs[itype*3 + 2]
       Line2 = Line1.replace('{', ' ')
       Line3 = Line2.replace('}', ' ')
       Line4 = Line3.replace(',', ' ')
       Line5 = Line4.replace('"', ' ')
       Line6 = Line5.split()
       ScatAtoms[itype][NGauss1]          = float(Line6[0])
       ScatAtoms[itype][NGauss1 + NGauss] = Badd

       TypesAtoms[itype] = AtomType
    return ScatAtoms, TypesAtoms

#============================

def SFactG(ScatAtom,Resolution,NSgrid) :
    ScatFunc = [0.0 for ig in range(NSgrid+1)]

    Smax    = 1.0 / Resolution
    dsstep  = Smax / NSgrid
    NGauss  = int(len(ScatAtom) / 2)

    for isg in range(NSgrid+1) :

        sg   = dsstep * isg
        ss24 = sg * sg / 4.0
        fact = 0.0

        for ig in range(NGauss) :
            argexp = ScatAtom[ig + NGauss] * ss24
            fact  += ScatAtom[ig] * math.exp(-argexp)

        ScatFunc[isg] = fact

    return ScatFunc

#============================

def AtomImage(ScatFunc,Resolution,DistImage,StepImage) :

    NImage = int(DistImage/StepImage) + 1
    NSGrid = len(ScatFunc) - 1
    Smax   = 1.0 / Resolution
    SStep  = Smax / NSGrid

    Image = [0.0 for j in range(NImage+1)]

    dx = 2. * math.pi * StepImage

#   integrate scattering curve

#   odd points

    for igs in range(1, NSGrid, 2):
        ss     = SStep * igs
        fatoms = ScatFunc[igs] * ss * 4.

        for ir in range(1,NImage):
            rr   = dx * ir
            arg  = rr * ss
            sarg = math.sin(arg)
            Image[ir] = Image[ir] + fatoms * sarg

        Image[0] = Image[0] + fatoms * ss

#   even points

    for igs in range(2, NSGrid-1, 2):
        ss     = SStep * igs
        fatoms = ScatFunc[igs] * ss * 2.

        for ir in range(1,NImage):
            rr   = dx * ir
            arg  = rr * ss
            sarg = math.sin(arg)
            Image[ir] = Image[ir] + fatoms * sarg

        Image[0] = Image[0] + fatoms * ss

#   terminal point (point s = 0 gives zero contribution and is ignored)

    ss     = SStep * NSGrid
    fatoms = ScatFunc[NSGrid] * ss

    for ir in range(1,NImage):
        rr   = dx * ir
        arg  = rr * ss
        sarg = math.sin(arg)
        Image[ir] = Image[ir] + fatoms * sarg

    Image[0] = Image[0] + fatoms * ss

# ---- normalisation ----

    scal = 2.0 * SStep / 3.0
    for ir in range(1,NImage):
        rr = ir * StepImage
        Image[ir] = Image[ir] * scal / rr

    Image[0] = Image[0] * SStep * 4. * math.pi / 3.

    return Image

#============================

def GetPrecision(Image, StepImage, DistMax) :

    NMax = int(DistMax/StepImage)
    NImage  = len(Image)

    for ig in range(NMax,NImage) :
        AbsImageValue = abs(Image[ig])
        if AbsImageValue >= abs(Image[ig-1]) and AbsImageValue > abs(Image[ig+1]) :
           LastMax = AbsImageValue
           break

    return LastMax

#============================

def CompleteTable(AtomType,Resolution,DistMax,Ngrid,NumTerms,
                  ErrorMin,BCR_Atom,ImageScale,Protocol,fileBCR) :

    ScaleD = Resolution
    ScaleR = Resolution
    ScaleB = Resolution * Resolution
    ScaleC = Resolution * ScaleB * ImageScale

    print(f'Atom {AtomType:4} Res {Resolution:8.4f} RhoMax {abs(ImageScale):9.4f} ',
          f'Dist {DistMax*ScaleD:8.4f} Ngrid {Ngrid:5} Nterms {NumTerms:2} ',
          f'ErrRel {ErrorMin:8.5f} ErrAbs {ErrorMin*ImageScale:8.5f} Prot {Protocol:3}',
          file = fileBCR)

    for iterms in range(NumTerms) :
        print(f'{iterms:3}{BCR_Atom[iterms][0]*ScaleR:16.10f}',
              f'{BCR_Atom[iterms][1]*ScaleB:16.10f}',
              f'{BCR_Atom[iterms][2]*ScaleC:16.10f}', file = fileBCR)

    return

#============================

def UpdatedBCR(BCR_Save_itype,Terms_itype) :

    rpeak = [BCR_Save_itype[iterms][0] for iterms in range(Terms_itype)]
    bpeak = [BCR_Save_itype[iterms][1] for iterms in range(Terms_itype)]
    cpeak = [BCR_Save_itype[iterms][2] for iterms in range(Terms_itype)]

    mpeak = Terms_itype - 1

    return rpeak, bpeak, cpeak, mpeak

#============================

def ext_BCR(kpres,epsc,epsp,edist,DistMax,Image,Distance,MaxTerms,
            BCR_Prev,Terms_Prev,Current,Start,fileBCRlog) :

#   there was no previous trial
    if Current <= Start :
       epsres = 1.0
       mpeak  = -1
       curres = Image
       curve  = Image
       bpeak, cpeak, rpeak  = [], [], []

#   use the BCR values updated from the previous resolution or from the previous atom
    else :

       ImageMax = 1.0

       rpeak, bpeak, cpeak, mpeak = UpdatedBCR(BCR_Prev,Terms_Prev)

       ceps,bmin,cmin,rmin,edist,mdist = \
            PreciseData(kpres,epsc,epsp,edist,DistMax,ImageMax,Distance,fileBCRlog)

       bpeak,cpeak,rpeak =               \
            RefineBCR(Image,Distance,mdist,edist,bpeak,cpeak,rpeak,mpeak,bmin,cmin,rmin,fileBCRlog)

       curres,curve,epsres =             \
            CurveDiff(Image,Distance,mdist,edist,fileBCRlog,bpeak,cpeak,rpeak,mpeak,MaxTerms)

    return bpeak,cpeak,rpeak,mpeak,curve,curres,epsres


#============================

def add_entry(results, d_min,
                       scattering_table,
                       RhoMax,
                       dist,
                       Ngrid,
                       Nterms,
                       ErrRel,
                       ErrAbs,
                       Prot,
                       B,C,R):
  results[str(d_min)] = {
      "scattering_table": scattering_table,
      "RhoMax"          : RhoMax          ,
      "dist"            : dist            ,
      "Ngrid "          : Ngrid           ,
      "Nterms"          : Nterms          ,
      "ErrRel"          : ErrRel          ,
      "ErrAbs"          : ErrAbs          ,
      "Prot"            : Prot            ,
      "R": R,
      "B": B,
      "C": C
  }

def RefinedTable(AtomType,Resolutions,DistMax,Ngrid,Terms_Resol,
                  ErrorMins,BCR_Resol,ImageMaxs,Protocol,fileBCR,
                  results=None, scattering_table=None) :

    Nresol = len(Resolutions)

    for ires in range(Nresol) :

        Resolution = Resolutions[ires]
        ImageScale = ImageMaxs[ires]
        NumTerms   = Terms_Resol[ires]
        ErrorMin   = ErrorMins[ires]

        ScaleD = Resolution
        ScaleR = Resolution
        ScaleB = Resolution * Resolution
        ScaleC = Resolution * ScaleB * ImageScale

        print(f'Atom {AtomType:4} Res {Resolution:8.4f} RhoMax {abs(ImageScale):9.4f} ',
            f'Dist {DistMax*ScaleD:8.4f} Ngrid {Ngrid:5} Nterms {NumTerms:2} ',
            f'ErrRel {ErrorMin:8.5f} ErrAbs {ErrorMin*ImageScale:8.5f} Prot {Protocol:3}',
            file = fileBCR)

        d_min  = Resolution
        RhoMax = abs(ImageScale)
        dist   = DistMax*ScaleD
        Ngrid  = Ngrid
        Nterms = NumTerms
        ErrRel = ErrorMin
        ErrAbs = ErrorMin*ImageScale
        Prot   = Protocol
        B,C,R  = [],[],[]

        for iterms in range(NumTerms) :
            R.append(BCR_Resol[ires][iterms][0]*ScaleR)
            B.append(BCR_Resol[ires][iterms][1]*ScaleB)
            C.append(BCR_Resol[ires][iterms][2]*ScaleC)
            print(f'{iterms:3}{BCR_Resol[ires][iterms][0]*ScaleR:16.10f}',
                f'{BCR_Resol[ires][iterms][1]*ScaleB:16.10f}',
                f'{BCR_Resol[ires][iterms][2]*ScaleC:16.10f}', file = fileBCR)
        if results is not None:
          add_entry(
            results          = results,
            d_min            = d_min,
            scattering_table = scattering_table,
            RhoMax           = RhoMax,
            dist             = dist,
            Ngrid            = Ngrid,
            Nterms           = Nterms,
            ErrRel           = ErrRel,
            ErrAbs           = ErrAbs,
            Prot             = Prot,
            B=B,C=C,R=R)

    return

#####################################################

def compute_tables(
      MinResolution = 4.0,
      MaxResolution = 4.0,
      DistMax       = 5.0,
      scattering_table = "wk1995",
      TypesAtoms = ["C",]):

  fileBCRlog = None #open('BCR.log', 'w')
  #                                        file to save residual curves when decomposition fails
  fileCurves = None #open('DiffBCR.curve', 'w')

  # parameters eventually NOT to be modified

  edist     = 1.0E-13              # term contribution limit value
  epsp      = 0.000                # peak limit value with respect to the central value
  kpres     = 1                    # precision to be defined in absolute values
  Badd      = 0.0                  # Badd (normally zero)
  Ngrid     = 2000                 # grid number for atomic images
  NSgrid    = 2000                 # grid number for scattering function
  Precision = 0.99                 # part of the local max next to the last; to estimate accuracy limit
  MaxTerms  = 50                   # max allowed number of terms
  validate  = 0.0010               # validation indictor : ErrorMax / ImageMax < validate is OK

  Protocols      = [0, 1, 111, 112, 122, 212]              # allowed protocol to be tried
  Nprotocols     = len(Protocols)

  CheckProtAll   = 100
  CheckProtFast  = 20

  RefineFlag = True

  # resolution bounds defining resolution increment for the resolution tables

  ResBound  = [ 0.6   , 0.8   , 1.0   , 2.0   , 100.0  ]
  AddArray  = [ 0.0020, 0.0025, 0.0040, 0.0050, 0.0100 ]


  #====================================================================================================

  #  MAIN PROGRAM

  # parameters which may be modified

  DistMax = DistMax                   # approximation up to DistMax * Resolution

  #==========================================================================

  #   define resolutions

  NumBound   = len(ResBound)
  Resolution = MinResolution
  if MaxResolution > ResBound[NumBound-1] : MaxResolution = ResBound[NumBound-1]

  for kbound in range(NumBound) :
      if MinResolution < ResBound[kbound] :
         UpperBound    = ResBound[kbound]
         AddRes        = AddArray[kbound]
         break

  Resolutions = []
  Resolution  = MinResolution
  while Resolution <= MaxResolution :
      Resolutions.append(Resolution)
      Resolution += AddRes
      if Resolution >= UpperBound * 0.999999 :
         kbound += 1
         if kbound == NumBound : break
         UpperBound = ResBound[kbound]
         AddRes     = AddArray[kbound]

  Nresol = len(Resolutions)

  #================== define common parameters

  BCR_Save    = [ [ [0.0, 0.0, 0.0] for iterm in range(MaxTerms)] for ires in range(Nresol) ]
  Terms_Save  = [ 0 for ires in range(Nresol) ]

  BCR_Resol    = [ [0.0, 0.0, 0.0] for ires in range(-1, Nresol) ]
  Terms_Resol  = [ 0   for ires in range(-1, Nresol) ]
  ErrorMins    = [ 0   for ires in range(Nresol) ]
  ImageMaxs    = [ 0.0 for ires in range(Nresol) ]

  #   atomic image will be generated up to the distance limit higher than that for the decomposition
  #   to eventually use an approximation also to the ripple(s) next beyond the distance limit
  #   atomic images are normalized by its value (max = 1.0) and by distance (resolution -> 1A)

  DistImage  = DistMax + 1.0
  StepImage  = DistMax / Ngrid

  #============== read table of multiGaussian approximation to scattering factors

  Type_Start = 0
  Ntypes     = len(TypesAtoms)
  Type_Final = Ntypes

  #==========================================================================

  for itype in range(Type_Start, Type_Final) :      #  INTERMEDIATE CYCLE OVER ATOMS

      AtomType = TypesAtoms[itype]

      print('')
      print('##############################################################################')
      print('')
      print(f'Calculations started for atom ',AtomType)

  #   generate file name for the coefficients for a given atomic type

      fileBCRname = 'BCR_%s_'%scattering_table + AtomType + '.temp'
      fileBCRtemp = open(fileBCRname, 'w')

      for ires in range(Nresol) :              # MAIN CYCLE OVER RESOLUTIONS

         Resolution = Resolutions[ires]

         print('')
         print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
         print('')
         print(f'Calculations started at resolution {Resolution:8.4f} for atom ',AtomType)

  ##      prepare scattering factors on a grid up to the chosen resolution limit
  #
  #       ScatFunc = SFactG(ScatAtoms[itype],Resolution,NSgrid)
  #
  ##      calculate atomic image and its maximal absolute value (expected to be at the origin)
  ##      scale Image to its max value to avoid optimization problems
  #
  #       Image = AtomImage(ScatFunc,Resolution,DistImage*Resolution,StepImage*Resolution)


         o = maptbx.atom_curves(scattering_type=AtomType, scattering_table=scattering_table)
         v = list(o.scr.as_type_gaussian_dict().values())[0]
         ff_AU_style = tuple(v.array_of_a()) + (v.c(),) + tuple(v.array_of_b()) + (0,)
         ff_AU_style = [round(_,6) for _ in ff_AU_style]
         ScatFunc = SFactG(ff_AU_style,Resolution,NSgrid)
         Image = AtomImage(ScatFunc,Resolution,DistImage*Resolution,StepImage*Resolution)

         NGridImage = len(Image)
         Distance   = [ id * StepImage for id in range(NGridImage)]
         DiffImage  = [ [ 0.0 for id in range(NGridImage)] for iprot in range(Nprotocols)]

         ImageScale = abs(Image[0])
         if Image[0] >= 0.0 :
            ImageSign = 1.0
         else :
            ImageSign = -1.0

         for id in range(NGridImage) :
             if ImageScale < abs(Image[id]) :
                ImageScale = abs(Image[id])
                if Image[id] >= 0.0 :
                   ImageSign  = 1.0
                else :
                   ImageSign  = -1.0
         ImageScale      = ImageScale * ImageSign
         ImageMaxs[ires] = ImageScale

         for id in range(NGridImage) :
             Image[id] = Image[id] / ImageScale

  #      find the value of the last extremum and define the decomposition precision

         LastMax  = GetPrecision(Image, StepImage, DistMax)
         epsc     = LastMax * Precision
         ErrorMin = 1.0

  #      define arrays

         ErrorMax = [ 0.0 for iprot in range(Nprotocols) ]
         Mterms   = [ 0   for iprot in range(Nprotocols) ]
         BCRterms = [ [ [0.0, 0.0, 0.0] for iterm in range(MaxTerms)] for iprot in range(Nprotocols) ]

         if ires == 0 or ires == int(ires/CheckProtAll) * CheckProtAll:
            NumProt = Nprotocols
         elif ires == int(ires/CheckProtFast) * CheckProtFast :
            NumProt = 4
         else :
            NumProt = 2

         for iprot in range(NumProt) :            #      INNER CYCLE OVER PROTOCOLS

             kprot = Protocols[iprot]

             print('')
             print('---------------------------------')
             print('')
             print(f'Resolution {Resolution:8.4f} Atom',AtomType,'Trying protocol',iprot,
                   f'( = ',kprot,') from',Nprotocols)

  #          try BCR from the previous atom at the same resolution

             if kprot == 0 :
                bpeak,cpeak,rpeak,mpeak,curve,curres,epsres = \
                    ext_BCR(kpres,epsc,epsp,edist,DistMax,Image,Distance,MaxTerms,
                            BCR_Save[ires],Terms_Save[ires],itype,Type_Start,fileBCRlog)

  #          try BCR from the previous resolution for the same atom

             elif kprot == 1 :

                bpeak,cpeak,rpeak,mpeak,curve,curres,epsres = \
                    ext_BCR(kpres,epsc,epsp,edist,DistMax,Image,Distance,MaxTerms,
                            BCR_Resol[ires-1],Terms_Resol[ires-1],ires,0,fileBCRlog)

  #          search for new values from scratch trying kprot > 100 protocols

             else :
                 bpeak,cpeak,rpeak,mpeak,curve,curres,epsres = \
                    get_BCR(Image,Distance,DistMax,MaxTerms,epsc,epsp,edist,kpres,kprot,fileBCRlog)

  #          save the error value, the number of terms, BCR coefficients and the difference curve;
  #          BCR subroutines count the BCR-triplets from 0 to mpeak including

             ErrorMax[iprot] = epsres
             Mterms[iprot]   = mpeak + 1

             for iterm in range(mpeak+1) :
                 BCRterms[iprot][iterm] = [ rpeak[iterm] , bpeak[iterm], cpeak[iterm] ]

             for id in range(NGridImage) :
                 DiffImage[iprot][id]  = curres[id]

  #          identify the best decomposition coefficients among the generated ones

             if ErrorMin > epsres  :
                ErrorMin = epsres
                mprot    = iprot

  #          if the multi-iteration protocol was required : test-run option ??

             if kprot > 200 :
                if fileCurves is not None:
                  print('atom, res, ImageMax, ErrorMin, ErrorRel',
                       f'{AtomType:4}{Resolution:5.2f}{ImageMax*ImageScale:12.4f}',
                       f'{ErrorMin*ImageScale:8.5f}{ErrorMin:8.5f}', file = fileCurves)
                  for id in range(NGridImage) :
                       print(f'{id:6}{Distance[id]*Resolution:9.5f}{Image[id]*ImageScale:10.5f}',
                             f'{DiffImage[nprot][id]*ImageScale:9.5f}',
                             f'{DiffImage[iprot][id]*ImageScale:9.5f}', file = fileCurves)

  #          all single-run protocols applied; check if the multi-iteration protocol is required
  #          (if it is in the list - supposed to be the last one)

             elif (Protocols[Nprotocols-1] > 200) and (iprot == (Nprotocols-2)) \
                                                   or (iprot == (Nprotocols-1)) :
                if ErrorMin <= validate:
                   break
                else :
                   print('insufficient accuracy : atom, res, ImageMax, ErrorMin, ErrorRel',
                         f'{AtomType:4}{Resolution:7.4f}{ImageMax*ImageScale:12.4f}',
                         f'{ErrorMin*ImageScale:8.5f}{ErrorMin:8.5f}')
                   nprot = mprot

  #      add the coefficients to the Table (write to the file)

         CompleteTable(AtomType,Resolution,DistMax,Ngrid,Mterms[mprot],
                   ErrorMin,BCRterms[mprot],ImageScale,Protocols[mprot],fileBCRtemp)

  #      save coefficients for a trial with the next resolution

         BCR_Save[ires]    = BCRterms[mprot]
         Terms_Save[ires]  = Mterms[mprot]

         BCR_Resol[ires]   = BCRterms[mprot]
         Terms_Resol[ires] = Mterms[mprot]
         ErrorMins[ires]   = ErrorMin

      fileBCRtemp.close()

  #==== END OF THE CYCLE OVER RESOLUTIONS =========

  #   if required refine the set of coefficients ; create a file with refined coefficients

      if RefineFlag :
         fileBCRname = 'BCR_%s_'%scattering_table + AtomType + '.table'
         fileBCR     = open(fileBCRname, 'w')

         results = {}

         ErrorMin = ErrorMins[0]
         mres     = 0
         for ires in range(Nresol) :
             if ErrorMin >= ErrorMins[ires] :
                ErrorMin  = ErrorMins[ires]
                mres = ires

         print('')
         print('---------------------------------')
         print('')
         print(f'Minimal error {ErrorMin:9.6f} for the resolution {Resolutions[mres]:8.4f}')

         Protocol = 1

         for ires in range(mres-1, -1, -1) :

  #      refine backward on the resolution starting from the best set

             Resolution = Resolutions[ires]
             ImageScale = ImageMaxs[ires]

             print('')
             print('---------------------------------')
             print('')
             print(f'Refine BCR : Resolution {Resolution:8.4f} Atom',AtomType)

             #ScatFunc = SFactG(ScatAtoms[itype],Resolution,NSgrid)
             #
             #Image = AtomImage(ScatFunc,Resolution,DistImage*Resolution,StepImage*Resolution)

             ScatFunc = SFactG(ff_AU_style,Resolution,NSgrid)
             Image = AtomImage(ScatFunc,Resolution,DistImage*Resolution,StepImage*Resolution)

             NGridImage = len(Image)
             Distance   = [ id * StepImage for id in range(NGridImage)]

             for id in range(NGridImage) :
                 Image[id] = Image[id] / ImageScale

             bpeak,cpeak,rpeak,mpeak,curve,curres,epsres = \
                  ext_BCR(kpres,epsc,epsp,edist,DistMax,Image,Distance,MaxTerms,
                          BCR_Resol[ires+1],Terms_Resol[ires+1],ires,-1,fileBCRlog)

             if epsres < ErrorMins[ires] :
                ErrorMins[ires]   = epsres
                Terms_Resol[ires] = mpeak + 1
                for iterm in range(mpeak+1) :
                    BCR_Resol[ires][iterm] = [ rpeak[iterm] , bpeak[iterm], cpeak[iterm] ]
             else :
                print('Resolution',ires,Resolution,'not improved',epsres,ErrorMins[ires])

         RefinedTable(AtomType,Resolutions,DistMax,Ngrid,Terms_Resol,
                ErrorMins,BCR_Resol,ImageMaxs,Protocol,fileBCR,
                results, scattering_table)

         fileBCR.close()
         with open("%s_%s.json"%(AtomType, scattering_table), "w") as f:
           json.dump(results, f, indent=2)

  #   ALL RESOLUTION HAVE BEEN PROCESSED FOR THE GIVEN ATOM; SWITCH TO THE NEXT ATOM

  print('finish')
