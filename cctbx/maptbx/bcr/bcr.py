from __future__ import absolute_import, division, print_function
import math
from scipy.optimize import minimize
from libtbx import adopt_init_args

from cctbx.array_family import flex
import scitbx.minimizers

# Adapted by P. Afonine from version 8.0 supplied by A. Urzhumtsev

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
def get_BCR(dens,dist,dmax,mxp,epsc,epsp=0.000,edist=1.0E-13,kpres=1,kprot=112,nfmes=None):
#
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

    #if (nfmes is not None) :
    #   print('',file=nfmes)
    #   print(30*'=',' NEW CURVE IN PROCESS ',30*'=',file=nfmes)
    #
    #   print('')
    #   print(30*'=',' NEW CURVE IN PROCESS ',30*'=')

    assert kprot in [111, 112, 122, 211, 212, 222]
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

#===========================================================



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
    if 1: # lbfgsb from cctbx
      #for it in range(1,10):
      for it in range(1,2):
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

#!!!!!!!!!!! end of the insertion

#   recover refined values

    #if (nfmes is not None) :
    #   print('',file=nfmes)
    #   print(' parameters before and after their refinement:'     ,file=nfmes)
    #   print(' Nterm    R (M)   ',6*' ','B (N)   ',6*' ','C (K)   ',5*' ',
    #                  'R (M)   ',6*' ','B (N)   ',6*' ','C (K)   ',file=nfmes)
    #
    #   print('')
    #   print(' parameters before and after their refinement:')
    #   print(' Nterm    R (M)   ',6*' ','B (N)   ',6*' ','C (K)   ',5*' ',
    #                   'R (M)   ',6*' ','B (N)   ',6*' ','C (K)   ')

    for i in range(npeak+1):
        i3 = i*3
        b0 , c0, r0 = res.x[i3] , res.x[i3+1] , res.x[i3+2]

        #if (nfmes is not None) :
        #   print(f'{i+1:6}{rpeak[i]:16.10f}{bpeak[i]:16.10f}{cpeak[i]:11.5f}   ',
        #         f'{r0:16.10f}{b0:16.10f}{c0:11.5f}',file=nfmes)
        #
        #   print(f'{i+1:6}{rpeak[i]:16.10f}{bpeak[i]:16.10f}{cpeak[i]:11.5f}   ',
        #         f'{r0:16.10f}{b0:16.10f}{c0:11.5f}')

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
         print('              shrter that the required distance  ',dmax)
    else :
       for ix in range(ndist+1) :
           if dist[ix] > dmax :
              mdist = ix - 1
              break

#   estimates for a regular grid with step = dist(1) - dist(0)

    dstep = dist[1]
    ndrop = 2.
    bmin  = 16.*dstep*dstep/ndrop
    cmin  = ceps*(bmin**1.5)/64.
    rmin  = 0.0

    #print ('',file=nfmes)
    #print (' INTERNAL DECOMPOSITION PARAMETERS ',file=nfmes)
    #print ('',file=nfmes)
    #print (' absolute max allowed error     ',f'{ceps:13.5e}' ,file=nfmes)
    #print (' drop to estimate the peak width',f'{epsp:13.5e}' ,file=nfmes)
    #if edist > 0.0 :
    #   print (' term extension limit           ',f'{edist:13.5e}',file=nfmes)
    #else :
    #   print (' term extension limit :            no limit applied',file=nfmes)
    #print (' estimated accuracy parameters :',file=nfmes)
    #print ('    bmin                        ',f'{bmin:13.5e}' ,file=nfmes)
    #print ('    cmin                        ',f'{cmin:13.5e}' ,file=nfmes)
    #print ('---------------------------',file=nfmes)

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

    #if (nfmes is not None) :
    #   space = 23*' '
    #   print('',file=nfmes)
    #   print(f' with {npeak+1:4} terms max residual peaks are',file=nfmes)
    #   print(space,' inside the interval of modeling',f'{epsres:12.7f}',file=nfmes)
    #   print(space,' in the whole range of distances',f'{totres:12.7f}',file=nfmes)
    #
    #   print('')
    #   print(f' with {npeak+1:4} terms max residual peaks are')
    #   print(space,' inside the interval of modeling',f'{epsres:12.7f}')
    #   print(space,' in the whole range of distances',f'{totres:12.7f}')

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
def OutCurves(dens,dist,curres,curve,filecurve):

    if filecurve == 'none' :
       return

    nfcurv = open(filecurve, 'w')

    ndist = len(dist) - 1

#    print(f'   dist   curve   modeled   resid    cmp=  0',
#          ''.join(f'{ip:9}' for ip in range(1,npeak+1)),file=nfcurv)
#    for i in range(ndist+1):
#        densum = dens[i]-curres[i]
#        print(f'{dist[i]:7.3f}{dens[i]:9.5f}{densum:9.5f}{curres[i]:9.5f}',
#              ''.join(f'{curve[i][ip]:9.5f}' for ip in range(npeak+1)),file=nfcurv)

    #print(f'    dist       curve         modeled         resid        max.error',file=nfcurv)

    resmax = 0.0
    for i in range(ndist+1):
        curabs = abs(curres[i])
        if curabs > resmax:
           resmax = curabs
        densum = dens[i]-curres[i]
        #print(f'{dist[i]:8.3f}{dens[i]:15.8f}{densum:15.8f}{curres[i]:15.8f}{resmax:15.8f}',
        #                                                                         file=nfcurv)
    return

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

#    exit()

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

#================================
def CheckLow(bpeak,cpeak,rpeak,npeak,ceps,nfmes):

    pi4    = 4.0 * math.pi
    pisq16 = pi4 * pi4

    bpeakc, cpeakc, rpeakc = [], [], []

    #if (nfmes is not None) :
    #
    #   print('',file=nfmes)
    #   print(' Check term values in the points r = r0 :',file=nfmes)
    #   print(' Nterm    R (M)   ',6*' ','B (N)   ',6*' ','C (K)   ',7*' ','value',file=nfmes)
    #
    #   print('')
    #   print(' Check term values in the points r = r0 :')
    #   print(' Nterm    R (M)   ',6*' ','B (N)   ',6*' ','C (K)   ',7*' ','value')

    for ip in range (npeak+1):
        b0 = bpeak[ip]
        c0 = cpeak[ip]
        r0 = rpeak[ip]
        if r0 == 0.0 :
           fval = c0 * (pi4/b0)**1.5
        else :
           fval = c0 * (1.0 - math.exp(-pisq16*r0**2/b0)) / (math.sqrt(pi4*b0) * r0**2)

        #if (nfmes is not None) :
        #   print(f'{ip+1:6}{r0:16.10f}{b0:16.10f}{c0:11.5f}   {fval:12.7f}',file=nfmes)
        #   print(f'{ip+1:6}{r0:16.10f}{b0:16.10f}{c0:11.5f}   {fval:12.7f}')

        if abs(fval) >= ceps :
           bpeakc.append(b0)
           cpeakc.append(c0)
           rpeakc.append(r0)

    return bpeakc, cpeakc, rpeakc

#============================
def FilterWeak(dens,curve,curres,dist,mdist,edist,epsres,bpeak,cpeak,rpeak,npeak,mxp,  \
                                 bmin,cmin,rmin,nfmes,ceps):

    bpeakc,cpeakc,rpeakc = CheckLow(bpeak,cpeak,rpeak,npeak,ceps,nfmes)

    mpeak = len(bpeakc) - 1

    if mpeak < npeak :

#      small contributions found; check the result after removing these terms

       bpeakc,cpeakc,rpeakc =        \
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
