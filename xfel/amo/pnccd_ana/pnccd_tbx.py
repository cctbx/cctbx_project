from __future__ import division
from psana import *
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import optimize
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


def circles(r, cent, dim, dr = 2):

    circ        = np.zeros(dim).astype(np.float64)
    y, x        = np.ogrid[0:dim[0], 0:dim[1]]
    index       = ((x-cent[0])**2 + (y-cent[1])**2 >= (r-dr)**2)&((x-cent[0])**2 + (y-cent[1])**2 <= (r+dr)**2)
    circ[index] = 1e6

    return circ

def polar2cart(r, theta, center):

    x = r  * np.cos(theta) + center[0]
    y = r  * np.sin(theta) + center[1]

    return x, y


def get_static_mask(msk,cent,r_max,msk_a,msk_w):
    """Returns static mask with specified angular sectors masked out

       @msk     Mask in cartesian coordinates                  [nx,ny]
       @cent    Beam center coordinates                        [pos.integer,pos.integer]
       @r_max   Max radial value (in pixels)                   [pos. integer]
       @mask_a  Centers of angular slice to mask               [pos.integer1 pos.integer2,...]
       @msk_w  Width of angular slice to mask                  [pos.integer1 pos.integer2,...]
    """

    dim = msk.shape
    x,y = np.ogrid[:dim[0],:dim[1]]
    cx,cy = cent

    for i in xrange(len(msk_a)) :

        tmin,tmax = np.deg2rad((int(msk_a[i])-int(msk_w[i]),int(msk_a[i])+int(msk_w[i])))

        # ensure stop angle > start angle
        if tmax < tmin:
                tmax += 2*np.pi

        # convert cartesian --> polar coordinates
        r2 = (x-cx)*(x-cx) + (y-cy)*(y-cy)
        theta = np.arctan2(x-cx,y-cy) - tmin

        # wrap angles between 0 and 2*pi
        theta %= (2*np.pi)

        # circular mask
        circmask = r2 <= r_max*r_max

        # angular mask
        anglemask = theta <= (tmax-tmin)

        ind = circmask * anglemask

        msk[ind] = 0

    return msk


def get_polar(img, msk, cent, r_max, r_min = 0, dr = 1, nPhi = None, dPhi = 1, msk_a = None, msk_w = None, plot = 0):
    """Returns cartesian images img & msk in polar coordinates pcimg[r,Phi] & pcmsk[r,Phi]

       @img     Image in cartesian coordinates                  [nx,ny]
       @msk     Mask in cartesian coordinates                   [nx,ny]
       @cent    Beam center coordinates                         [pos.integer,pos.integer]
       @r_max   Max radial value (in pixels)                    [pos. integer]
       @r_min   Min radial value (in pixels)                    [pos. integer]
       @dr      Stepsize in radius                              [pos.integer]
       @nPhi    Number of Phi points                            [pos.integer]
       @dPhi    Stepsize in Phi                                 [pos.integer]
       @mask_a  Centers of angular slice to mask                [pos.integer1 pos.integer2,...]
       @msk_w  Width of angular slice to mask                  [pos.integer1 pos.integer2,...]
       @plot    Display result of grid search                   [0/1]

    """

    if msk is None:
       msk = np.ones((img.shape))

    if nPhi is None :
       nPhi  = np.ceil(2*np.pi*r_max)

    theta , R = np.meshgrid(np.linspace(0, 2*np.pi, nPhi/dPhi,endpoint=False),np.arange(r_min, r_max, dr))

    Xcart, Ycart = polar2cart(R, theta, cent)

    Xcart = Xcart.round().astype(int)
    Ycart = Ycart.round().astype(int)

    pcimg = img[Ycart,Xcart]
    pcimg = np.reshape(pcimg,((r_max-r_min)/dr,nPhi/dPhi))

    pcmsk = msk[Ycart,Xcart]
    pcmsk = np.reshape(pcmsk,((r_max-r_min)/dr,nPhi/dPhi))

    if (msk_a is not None) and (msk_w is not None) :
       # Mask out selected angles
       for i in xrange(len(msk_a)) :
           pcmsk[:,(int(msk_a[i])-int(msk_w[i])):(int(msk_a[i])+int(msk_w[i]))] = 0


    if plot :

       plt.figure(2000)
       plt.clf()
       ax = plt.subplot(221)
       plt.imshow(img)
       plt.axis('image')
       plt.clim(0,100)
       plt.colorbar()
       ax.set_title("Cartesian image")
       ax = plt.subplot(222)
       plt.imshow(msk)
       plt.axis('image')
       plt.clim(0,1)
       plt.colorbar()
       ax.set_title("Cartesian mask")
       ax = plt.subplot(223)
       plt.imshow(pcimg)
       plt.axis('tight')
       plt.clim(0,100)
       plt.colorbar()
       plt.ylabel('r')
       plt.xlabel('Phi')
       ax.set_title("Polar image")
       ax = plt.subplot(224)
       plt.imshow(pcmsk)
       plt.axis('tight')
       plt.clim(0,1)
       plt.colorbar()
       plt.ylabel('r')
       plt.xlabel('Phi')
       ax.set_title("Polar mask")
       plt.draw()

    return pcimg, pcmsk




def get_beam(img, msk, r_max, dr = 1, cent0 = None, dx = 0, dy = 0, ang = 45, dang = 10, plot = 0 ) :
    """Returns estimated beam center coordintes (cent) refined using
       a grid search assuming Friedel symmetry in the image

       @img     Image in cartesian coordinates                  [nx,ny]
       @msk     Mask in cartesian coordinates                   [nx,ny]
       @r_max   Max radial value (in pixels) used in refinement [pos. integer]
       @dr      Stepsize in radius                              [pos.integer]
       @cent0   Initial guess of beam center coordinates        [pos.integer,pos.integer]
       @dx      Search in x-direction, x+/-dx                   [pos. integer]
       @dy      Search in y-direction, y+/-dy                   [pos. integer]
       @ang     Center of angular slice in degrees              [pos. integer]
       @dang    Size of angular slice, ang+/-dang               [pos. integer]
       @plot    Display result of grid search                   [0/1]

    """

    # Default start from center of gravity of the image
    if cent0 is None:
       cent0 = [int(round(img.shape[0]/2)) , int(round(img.shape[1]/2))]


    if (dx + dy) == 0:

       cent = cent0

    else:

       na  = np.ceil(2*np.pi*r_max)                     # Get the necessay number of phi slices

       if (na % (360/ang)):
          na  = np.ceil(na/(360/ang))*(360/ang)

       da     = round((dang/360)*na)                    # Get dang in phi slices

       # Get indices for the 4 slices
       ind_q  = np.arange(0, r_max/dr, 1).astype(int)
       ind_a1 = np.arange(((ang/360)*na-da),       ((ang/360)*na+da)+1, 1).astype(int)
       ind_a2 = np.arange((((ang+180)/360)*na-da), (((ang+180)/360)*na+da)+1, 1).astype(int)
       ind_a3 = np.arange((((ang+90)/360)*na-da),  (((ang+90)/360)*na+da)+1, 1).astype(int)
       ind_a4 = np.arange((((ang+270)/360)*na-da), (((ang+270)/360)*na+da)+1, 1).astype(int)

       C1     = np.ndarray([(2*dx+1),(2*dy+1)])
       C2     = np.ndarray([(2*dx+1),(2*dy+1)])

       x      = 0

       # Step through grid and score correlation between slices seperated by 180 deg
       for xx in xrange(cent0[0]-dx,cent0[0]+dx+1) :
           y = 0
           for yy in xrange(cent0[1]-dy,cent0[1]+dy+1) :

                  pcimg, pcmsk = get_polar(img, msk, [xx, yy], r_max, r_min = 0, dr = dr, nPhi = na, dPhi = 1 )

                  ind1 = pcmsk[np.ix_(ind_q,ind_a1)].astype(int)
                  ind2 = pcmsk[np.ix_(ind_q,ind_a2)].astype(int)
                  ind3 = pcmsk[np.ix_(ind_q,ind_a3)].astype(int)
                  ind4 = pcmsk[np.ix_(ind_q,ind_a4)].astype(int)

                  S1   = np.ma.average(pcimg[np.ix_(ind_q,ind_a1)],axis=1,weights=ind1)
                  S2   = np.ma.average(pcimg[np.ix_(ind_q,ind_a2)],axis=1,weights=ind2)
                  S3   = np.ma.average(pcimg[np.ix_(ind_q,ind_a3)],axis=1,weights=ind3)
                  S4   = np.ma.average(pcimg[np.ix_(ind_q,ind_a4)],axis=1,weights=ind4)

                  C1[x,y]   = np.corrcoef(S1,S2)[0,1]
                  C2[x,y]   = np.corrcoef(S3,S4)[0,1]

                  y += 1
           x += 1

       Ctot     = (C1/C1.max() + C2/C2.max())/2         # Equal weights for the 2 pairs for now,

       xcoord   = np.arange(cent0[0]-dx,cent0[0]+dx+1,1)
       ycoord   = np.arange(cent0[1]-dy,cent0[1]+dy+1,1)

       index    = np.unravel_index(Ctot.argmax(), Ctot.shape)
       cent     = [xcoord[index[0]] , ycoord[index[1]] ]

       if plot :

          # Generate rings from beam center

          img2 = img.copy()


          # Switch x,y for display purpose only since imshow() displays the transpose
          img2[(cent[1]-10):(cent[1]+10),cent[0]] = 1e6
          img2[cent[1],(cent[0]-10):(cent[0]+10)] = 1e6

          for rr in xrange(50,int(img2.shape[0]/2),50) :
              circ = circles(rr,cent,img2.shape)
              img2  = img2 + circ

          Y, X = np.meshgrid(ycoord, xcoord)


          fig=plt.figure(1000)
          plt.clf()
          ax = plt.subplot(211)
          plt.imshow(img2)
          plt.axis('tight')
          plt.clim(0,100)
          plt.colorbar()
          plt.xlabel('x')
          plt.ylabel('y')
          ax.set_title("Cartesian image")
          ax = fig.add_subplot(2,1,2,projection='3d')
          surf = ax.plot_surface(X, Y, Ctot, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
          ax.set_zlim(0, 1.01)
          ax.zaxis.set_major_locator(LinearLocator(10))
          ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
          fig.colorbar(surf, shrink=0.5, aspect=5)
          plt.xlabel('x')
          plt.ylabel('y')
          plt.axis('tight')
          ax.set_title("Coordinate Surface")
          plt.draw()

    return cent


def besselfit(k, *args) :

    data = args[0]
    qs   = args[1]

    theory = k[0]*(((3*((np.sin(qs*k[1]) - (qs*k[1]*np.cos(qs*k[1])))/((qs*k[1])**3))))**2);

    T      = np.corrcoef(data,theory)[0,1]
    T      = 1/T

    return T

def get_size(saxs, q, r_i, q_i = None, q_f = None, plot = 0) :
    """Spherical Besselfit of SAXS data starting from radius, r_i
       Returns optimized radius r_f in Angstrom

       @saxs    Saxs data                                       [nQ]
       @q       Q-vector                                        [nQ]
       @r_i     Initial guess for radius (in Ang)               [float]
       @q_i     Lower q-bound (Ang-1)                           [float]
       @q_f     Upper q-bound (Ang-1)                           [float]
       @plot    Display fit between data & theory               [0/1]

    """

    if q_i is None:
       q_i = q[int(np.argmax(saxs>0))]


    if q_f is None:
       q_f = q.max()

    # Get data for q-range
    ind  = (q >= q_i) & (q <= q_f)
    data = saxs[ind]
    qs   = q[ind]

    # Simplex optimization
    i_f,r_f  = sp.optimize.fmin(besselfit, x0=[data[0],r_i],args=(data,qs),disp=0)

    # Optimal solution
    theory   = i_f*(((3*((np.sin(qs*r_f) - (qs*r_f*np.cos(qs*r_f)))/((qs*r_f)**3))))**2);
    score    = np.corrcoef(data,theory)[0,1]

#    if score < 0.98 :
#       print 'Poor size estimate'

    if plot :

       scale = np.nanmean(data)/np.nanmean(theory)

       fig   = plt.figure(3000)
       plt.semilogy(qs,data,'b-',qs,scale*theory,'r-',linewidth=1)
       plt.xlabel('q (A^-1)')

    return r_f, score


def get_symmetry(c2, phi_lim = None, plot = 0) :
    """Compute pi/2 symmetry for q-ring

       @c2      C2                                              [nPhi]
       @phi_lim Truncate edge                                   [int]

    """

    n    = round(len(c2)/4)
    a    = c2[0:n]
    b    = c2[(n+1):2*n]
    brev = b[::-1]

    if phi_lim is none:
       p1 = 0
       p2 = len(brev)
    else:
       p1 = 0 + phi_lim
       p2 = len(brev)-phi_lim

    sym  =  np.corrcoef(a[p1:p2],brev[p1:p2])[0,1]

    if plot :

       plt.figure(6000)
       plt.clf()
       plt.plot(a,'bx-',brev,'rx-')
       plt.axis('tight')

    return


def get_bl(q, phi, saxs, c2, lambd, q_lim = None , phi_lim = None, l_max = 40) :
    """Determine Bl coeffecients by Legendre Polynomial
       decomposition of the angular autocorrelations (C2) after
       correcting for the curvature of the Ewald sphere.
       B0 comes from the square of the SAXS data.
       Returns even and odd coeffs and the fitted C2's.

       @q       Q-vector (Ang^-1)                               [nQ]
       @phi     Phi-vector (rad)                                [nPhi]
       @saxs    Mean intensity                                  [nQ]
       @c2      Mean normalized c2's                            [nQ x nPhi]
       @lambd   Wavelength (Ang)                                [float]
       @q_lim   Q-lim (pixels)                                  [int,int]
       @phi_lim Phi-lim (pixels)                                [int,int]
       @l_max   Maximum nr of Leg components to use             [int]

    """

def get_psana_event(file_stamps) :
    """Generate psana timestamps from time-stamps

       @file_stamps  Parameter file containing  Time, Seconds, Nanoseconds,Fiducials

    """

    timestamps = []

    for i in xrange(len(file_stamps)) :
        sec   = int(file_stamps[i][1])
        nsec  = int(file_stamps[i][2])
        fid   = int(file_stamps[i][3])

        et   = EventTime(int((sec<<32)|nsec),fid)

        timestamps.append(et)

    return timestamps
