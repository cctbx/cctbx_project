from __future__ import absolute_import, division, print_function
from six.moves import range
from psana import *
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from mpl_toolkits.mplot3d import Axes3D                 # implicit import
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

def dynamic_flatfield(pcimg, pcmsk = None):

    """Compute the dynamic flatfield image from an average image in polar coordinates
    (see Hosseinizadeh et. al, Structural Dynamics, 2015)

       @pcimg   Average Image in polar coordinates
       @pcmsk   Average Mask  in polar coordinates

    """

    # Find non-gap indicices
    if pcmsk is None :
       pcmsk       = np.zeros(pcimg.shape)
       ind         = pcimg != 0
       pcmsk[ind]  = 1
    else:
       ind         = pcmsk == 1

    # Comupte saxs data
    s           = np.ma.average(pcimg,axis=1,weights=pcmsk)
    S           = np.zeros(pcimg.shape)
    for q in range(S.shape[0]) :
        S[q,:]  = s[q]

    # Compute Flatfield
    pcflat      = np.zeros(pcimg.shape)
    pcflat[ind] = S[ind]/pcimg[ind]

    return pcflat


def analyze_quadrants(pcimg, pcmsk, qrange = None, plot = 0):

    """Analyze intensity distribution in each detector quadrant with the
       purpose to identify outlier image quadrants.

       @pcimg   Image in polar coordinates
       @pcmsk   Mask in polar coordinates
       @qrange  Q-range for analysis

    """

    # Nr of q and phi
    nQ     = pcimg.shape[0]
    nPhi   = pcimg.shape[1]

    if qrange is None:
       qrange[0] = 0
       qrange[1] = nQ

    # Quadrant angles
    ang_1  = 0
    ang_2  = int(nPhi/4)
    ang_3  = int(nPhi*(2/4))
    ang_4  = int(nPhi*(3/4))

    # Compute saxs data for each quadrant
    S1     = np.ma.average(pcimg[:,ang_1:ang_2],axis=1,weights=pcmsk[:,ang_1:ang_2])
    S2     = np.ma.average(pcimg[:,ang_2:ang_3],axis=1,weights=pcmsk[:,ang_2:ang_3])
    S3     = np.ma.average(pcimg[:,ang_3:ang_4],axis=1,weights=pcmsk[:,ang_3:ang_4])
    S4     = np.ma.average(pcimg[:,ang_4:],axis=1,weights=pcmsk[:,ang_4:])

    # Compute mean and std for each quadrant
    Means  = []
    Means.append(np.median(S1[qrange[0]:qrange[1]]))
    Means.append(np.median(S2[qrange[0]:qrange[1]]))
    Means.append(np.median(S3[qrange[0]:qrange[1]]))
    Means.append(np.median(S4[qrange[0]:qrange[1]]))
    Means  = np.array(Means)

    # Estimate variation in % between quadrants
    Var    = 100*(np.std(Means)/np.mean(Means))

    if plot :

       m1 = np.nanmean(pcimg.reshape((-1)))
       s1 = np.std(pcimg.reshape((-1)))

       q  = np.arange(nQ)
       plt.figure(10000,figsize=(10,10))
       plt.clf()
       ax = plt.subplot(211)
       plt.imshow(pcimg)
       plt.axis('tight')
       plt.clim(0,m1+3*s1)
       plt.colorbar()
       ax = plt.subplot(212)
       plt.semilogy(q[qrange[0]:qrange[1]],S1[qrange[0]:qrange[1]],'b',label='I')
       plt.semilogy(q[qrange[0]:qrange[1]],S2[qrange[0]:qrange[1]],'r',label='II')
       plt.semilogy(q[qrange[0]:qrange[1]],S3[qrange[0]:qrange[1]],'g',label='III')
       plt.semilogy(q[qrange[0]:qrange[1]],S4[qrange[0]:qrange[1]],'m',label='IV')
       plt.legend()
       plt.draw()

    return Var


def pixel_mask(imgs, thr = 1.0, plot = 1):

    """Compute pixel mask by identifying non-responsive or overloaded pixels from a set of images
       by assuming that the variance roughly follows Possion distribution

       @imgs    Assmbled image set ex: on the format 1024 x 1024 x N elements for pnCCDs
       @thr     Threshold for identifying non-changing pixels, default 1.0 * pixel std**2

    """

    # Compute average and std intensity for each pixel
    mean_image  = np.mean(imgs,axis=2)
    std_image   = np.std(imgs,axis=2)

    # Set pixels with a mean larger then thr*variance to zeros, for Poisson mean = sigma**2
    T           = abs(mean_image) < thr*(std_image**2)
    mask        = T.astype(int)

    if plot:

       m1 = np.nanmean(mean_image.reshape((-1)))
       s1 = np.std(mean_image.reshape((-1)))

       plt.figure(1200,figsize=(10,20))
       plt.clf()
       ax = plt.subplot(211)
       plt.imshow(mean_image)
       plt.axis('image')
       plt.clim(0,m1+3*s1)
       ax.set_title("Average image")
       ax = plt.subplot(212)
       plt.imshow(mask)
       plt.axis('image')
       ax.set_title("Pixel mask")
       plt.draw()

       # Display image for secs
       import time
       secs = 60
       time.sleep(secs)

    return mean_image,std_image**2, mask




def common_mode_hart(img, msk, max_int, max_com, length, orient = 0, plot = 0):

    """Apply common mode correction according to Philip Hart on raw image

       @img     Assmbled image ex: on the format 1024 x 1024 elements for pnCCDs
       @msk     Assmbled mask ex: on the format 1024 x 1024 elements for pnCCDs
       @max_int Max threshold for pixel intensity to be used in CM
       @max_com Max CM value for correction to be made
       @length  Length of consecutive pixel array, det. specific. For pnCCD 128 channel row

    """
    img2 = np.copy(img)

    dim1 = img.shape[0]
    dim2 = img.shape[1]

    # Get pixel bin indices
    bin_index = np.arange(0, dim2+length, length)

    if orient :

       # Loop over rows
       for r in range(dim1):
           # Loop over pixel bins
           for b in range(len(bin_index)-1):
               # Only use "live" pixels below max_int
               m    = msk[r,bin_index[b]:bin_index[b+1]]
               i    = img[r,bin_index[b]:bin_index[b+1]]
               m_dx = m > 0
               i_dx = i < max_int
               idx  = m_dx*i_dx
               cm   = np.median(i[idx])
               if cm < max_com:
                  img2[r,bin_index[b]:bin_index[b+1]] = img[r,bin_index[b]:bin_index[b+1]] - cm

    else:

       # Loop over rows
       for r in range(dim1):
           # Loop over pixel bins
           for b in range(len(bin_index)-1):
               # Only use "live" pixels below max_int
               m    = msk[bin_index[b]:bin_index[b+1],r]
               i    = img[bin_index[b]:bin_index[b+1],r]
               m_dx = m > 0
               i_dx = i < max_int
               idx  = m_dx*i_dx
               cm   = np.median(i[idx])
               if cm < max_com:
                  img2[bin_index[b]:bin_index[b+1],r] = img[bin_index[b]:bin_index[b+1],r] - cm


    if plot:

       m1 = np.nanmean(img.reshape((-1)))
       s1 = np.std(img.reshape((-1)))

       m2 = np.nanmean(img2.reshape((-1)))
       s2 = np.std(img2.reshape((-1)))

       plt.figure(1100,figsize=(10,20))
       plt.clf()
       ax = plt.subplot(311)
       plt.imshow(img)
       plt.axis('image')
       plt.clim(0,m1+3*s1)
       plt.colorbar()
       ax.set_title("Raw")
       ax = plt.subplot(312)
       plt.imshow(img2)
       plt.axis('image')
       plt.clim(0,m2+3*s2)
       plt.colorbar()
       ax.set_title("CM corrected")
       ax = plt.subplot(313)
       plt.imshow(img-img2)
       plt.axis('image')
       plt.clim(m2-3*s2,m2+3*s2)
       plt.colorbar()
       ax.set_title("Raw - CM")
       plt.draw()

    return img2


def common_mode(img, edge, side, plot = 0):
    """Apply quasi-common mode correction to raw image

       @img     Assmbled image ex: on the format 1024 x 1024 elements for pnCCDs
       @edge    Distance of box from edge
       @side    Box size for common mode probe area

    """
    img2 = np.copy(img)

    dim1 = img.shape[0]
    dim2 = img.shape[1]

    # Upper Left Quadrant
    Q1 = np.median(img[edge:(edge+side),edge:(edge+side)]);
    # Lower Left Quadrant
    Q2 = np.median(img[dim1-(edge+side):(dim1-edge),edge:(edge+side)]);
    # Upper Right Quadrant
    Q3 = np.median(img[edge:(edge+side),dim2-(edge+side):(dim2-edge)]);
    # Lower Right Quadrant
    Q4 = np.median(img[dim1-(edge+side):(dim1-edge),dim2-(edge+side):(dim2-edge)]);

    # Implement common mode
    img2[0:int(dim1/2),0:int(dim2/2)] = img[0:int(dim1/2),0:int(dim2/2)] - Q1
    img2[int(dim1/2):,0:int(dim2/2)]  = img[int(dim1/2):,0:int(dim2/2)] -  Q2
    img2[0:int(dim1/2),int(dim2/2):]  = img[0:int(dim1/2),int(dim2/2):] -  Q3
    img2[int(dim1/2):,int(dim2/2):]   = img[int(dim1/2):,int(dim2/2):]  -  Q4

    if plot:
       img3 = np.copy(img)

       # Display area used for cm estimate
       img3[edge:(edge+side),edge:(edge+side)]                          = 1e6
       img3[dim1-(edge+side):(dim1-edge),edge:(edge+side)]              = 1e6
       img3[edge:(edge+side),dim2-(edge+side):(dim2-edge)]              = 1e6
       img3[dim1-(edge+side):(dim1-edge),dim2-(edge+side):(dim2-edge)]  = 1e6


       m1 = np.nanmean(img.reshape((-1)))
       s1 = np.std(img.reshape((-1)))

       m2 = np.nanmean(img2.reshape((-1)))
       s2 = np.std(img2.reshape((-1)))

       plt.figure(1100,figsize=(10,20))
       plt.clf()
       ax = plt.subplot(311)
       plt.imshow(img3)
       plt.axis('image')
       plt.clim(0,m1+3*s1)
       plt.colorbar()
       ax.set_title("Raw")
       ax = plt.subplot(312)
       plt.imshow(img2)
       plt.axis('image')
       plt.clim(0,m2+3*s2)
       plt.colorbar()
       ax.set_title("CM corrected")
       ax = plt.subplot(313)
       plt.imshow(img-img2)
       plt.axis('image')
       plt.clim(m2-3*s2,m2+3*s2)
       plt.colorbar()
       ax.set_title("Raw - CM")
       plt.draw()

    return img2


def extend_image(img) :

    dim1  = img.shape[0]
    dim2  = img.shape[1]
    # Generate extended image and mask
    if dim1 > dim2:
       img2                     = np.zeros((dim1,dim1))
       edge                     = int((dim1 - dim2)/2)
       img2[:,edge:(edge+dim2)]  = img
    else:
       img2                     = np.zeros((dim2,dim2))
       edge                     = int((dim2 - dim1)/2)
       img2[edge:(edge+dim1),:]  = img

    return img2


def get_geometry(img, gap, shift, orient) :
    """Apply geometry; gap and shift to image

       @img     Assmbled image ex: on the format 1024 x 1024 elements for pnCCDs
       @gap     Gap in pixels   [pos. int]
       @shift   Shift in pixels [int]
       @orient  Orientation of panel pairs [0/1]
                0 : Gap/shift relative to upper panel pairs
                1 : Gap/shift relative to left panel pairs
    """

    # Retrieve image dimensions
    dim1      = img.shape[0]
    dim2      = img.shape[1]

    # Implement gap & shift between panel pairs
    if orient == 0 :
       # Upper-Lower orientation
       image = np.zeros((dim1+gap,dim2+abs(shift)))

       if shift >= 0 :
          image[0:int(dim1/2),0:dim2]                       =  img[0:int(dim1/2),:]
          image[int(dim1/2)+gap:,shift:(dim2+shift)]        =  img[int(dim1/2):,:]

       else:
          image[0:int(dim1/2),abs(shift):(dim2+abs(shift))] =  img[0:int(dim1/2),:]
          image[int(dim1/2)+gap:,0:dim2]                    =  img[int(dim1/2):,:]

    else :
       # Left-Right orientation
       image = np.zeros((dim1+abs(shift),dim2+gap))

       if shift >= 0 :
          image[0:dim1,0:int(dim2/2)]                       =  img[:,0:int(dim2/2)]
          image[shift:(dim1+shift),int(dim2/2)+gap:]        =  img[:,int(dim2/2):]

       else:
          image[abs(shift):(dim1+abs(shift)),0:int(dim2/2)] =  img[:,0:int(dim2/2)]
          image[0:dim1,int(dim2/2)+gap:]                    =  img[:,int(dim2/2):]

    return image


def assemble_image(img, msk = None) :
    """Assemble raw or calib format pnCCD panels

       @img     Un-assmbled image on the format 4 x 512 x 512 elements for pnCCDs
       @msk     Un-assmbled mask on  the format 4 x 512 x 512 elements for pnCCDs
    """

    if msk is None:
       msk = np.ones(img.shape)

    image  = np.zeros((1024,1024))
    mask   = np.zeros((1024,1024))

    panA   = np.zeros((512,512))
    panB   = np.zeros((512,512))
    panC   = np.zeros((512,512))
    panD   = np.zeros((512,512))

    masA   = np.zeros((512,512))
    masB   = np.zeros((512,512))
    masC   = np.zeros((512,512))
    masD   = np.zeros((512,512))

    x=0;
    for i in range(4) :
        for j in range(512) :
            for k in range(512) :
                if i == 0:
                   panA[j,k] = img[x]
                   masA[j,k] = msk[x]
                elif i == 1:
                   panB[j,k] = img[x]
                   masB[j,k] = msk[x]
                elif i == 2:
                   panC[j,k] = img[x]
                   masC[j,k] = msk[x]
                else:
                   panD[j,k] = img[x]
                   masD[j,k] = msk[x]
                x += 1

    # Assemble image, specifik to experiment
    image[0:512,0:512] = np.transpose(panC)
    image[0:512,512:]  = np.rot90(np.transpose(panD),k=2) # Rot 180
    image[512:,0:512]  = np.transpose(panB)
    image[512:,512:]   = np.rot90(np.transpose(panA),k=2) # Rot 180

    # Assemble mask, specifik to experiment
    mask[0:512,0:512] = np.transpose(masC)
    mask[0:512,512:]  = np.rot90(np.transpose(masD),k=2) # Rot 180
    mask[512:,0:512]  = np.transpose(masB)
    mask[512:,512:]   = np.rot90(np.transpose(masA),k=2) # Rot 180

    return image, mask


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

    for i in range(len(msk_a)) :

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

#        msk[ind] = 0
        msk[ind] = 1

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
       for i in range(len(msk_a)) :
           pcmsk[:,(int(msk_a[i])-int(msk_w[i])):(int(msk_a[i])+int(msk_w[i]))] = 0


    if plot :

       m1 = np.nanmean(img.reshape((-1)))
       s1 = np.std(img.reshape((-1)))

       m2 = np.nanmean(pcimg.reshape((-1)))
       s2 = np.std(pcimg.reshape((-1)))

       plt.figure(2000)
       plt.clf()
       ax = plt.subplot(221)
       plt.imshow(img)
       plt.axis('image')
       plt.clim(0,m1+3*s1)
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
       plt.clim(0,m2+3*s2)
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



def variance_norm(pcimg,pcmsk) :
    """Compute the variance normalized image in polar coordinates

       @pcimg   Image in polar coordinates
       @pcmsk   Mask in polar  coordinates

    """

    saxs_m           = np.zeros(pcimg.shape[0])
    saxs_s           = np.zeros(pcimg.shape[0])
    pcnorm           = np.zeros(pcimg.shape)

    for q in range(pcimg.shape[0]) :

        ind               = np.nonzero(pcmsk[q,:])

        if len(ind[0]) == 0 :
           saxs_m[q]     = 0
           saxs_s[q]     = 0
        else:
           saxs_m[q]     = np.nanmean(pcimg[q,ind],axis=1)
           saxs_s[q]     = np.nanstd(pcimg[q,ind],axis=1)

        if saxs_s[q] != 0 :

           pcnorm[q,ind] = (pcimg[q,ind] - saxs_m[q]) / saxs_s[q]

    pcnorm         = pcnorm*pcmsk

    return pcnorm,saxs_m,saxs_s



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
       for xx in range(cent0[0]-dx,cent0[0]+dx+1) :
           y = 0
           for yy in range(cent0[1]-dy,cent0[1]+dy+1) :

                  pcimg, pcmsk = get_polar(img, msk, [xx, yy], r_max, r_min = 0, dr = dr, nPhi = na, dPhi = 1 )

                  ind1 = pcmsk[np.ix_(ind_q,ind_a1)].astype(int)
                  ind2 = pcmsk[np.ix_(ind_q,ind_a2)].astype(int)
                  ind3 = pcmsk[np.ix_(ind_q,ind_a3)].astype(int)
                  ind4 = pcmsk[np.ix_(ind_q,ind_a4)].astype(int)

                  S1   = np.ma.average(pcimg[np.ix_(ind_q,ind_a1)],axis=1,weights=ind1)
                  S2   = np.ma.average(pcimg[np.ix_(ind_q,ind_a2)],axis=1,weights=ind2)
                  S3   = np.ma.average(pcimg[np.ix_(ind_q,ind_a3)],axis=1,weights=ind3)
                  S4   = np.ma.average(pcimg[np.ix_(ind_q,ind_a4)],axis=1,weights=ind4)


                  # Ascert that only non-zero intensities are compared
                  ind1 = (S1!=0)*(S2!=0)
                  ind2 = (S3!=0)*(S4!=0)

                  C1[x,y]   = np.corrcoef(S1[ind1],S2[ind1])[0,1]
                  C2[x,y]   = np.corrcoef(S3[ind2],S4[ind2])[0,1]

                  y += 1
           x += 1

       Ctot     = (C1/C1.max() + C2/C2.max())/2         # Equal weights for the 2 pairs for now,

       xcoord   = np.arange(cent0[0]-dx,cent0[0]+dx+1,1)
       ycoord   = np.arange(cent0[1]-dy,cent0[1]+dy+1,1)

       index    = np.unravel_index(Ctot.argmax(), Ctot.shape)
       cent     = [xcoord[index[0]] , ycoord[index[1]] ]

       if plot :

          print(cent)

          # Generate rings from beam center

          img2 = img.copy()
          img2 = img2*msk


          # Switch x,y for display purpose only since imshow() displays the transpose
          img2[(cent[1]-10):(cent[1]+10),cent[0]] = 1e6
          img2[cent[1],(cent[0]-10):(cent[0]+10)] = 1e6

          for rr in range(50,int(img2.shape[0]/2),50) :
              circ = circles(rr,cent,img2.shape)
              img2  = img2 + circ

          Y, X = np.meshgrid(ycoord, xcoord)


          m1 = np.nanmean(img.reshape((-1)))
          s1 = np.std(img.reshape((-1)))

          fig=plt.figure(1000,figsize=(5,20))
          plt.clf()
          ax = plt.subplot(211)
          plt.imshow(img2)
          plt.axis('tight')
          plt.clim(0,m1+3*s1)
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


def get_bessel(Q, R) :

    B    = 3*((np.sin(Q*R) - (Q*R*np.cos(Q*R)))/((Q*R)**3))

    return B


def bessel_fit(k, *args) :
    """Spherical Bessel fit by refining I0 and Particle Radius
   """

    data = args[0]
    qs   = args[1]

    bessel = get_bessel(qs,k[1])
    theory = k[0]*(bessel**2)

    T      = np.corrcoef(data,theory)[0,1]
    T      = 1/T

    return T


def bessel_fit2(k, *args) :
    """Spherical Bessel fit by refining I0 and Detector distance
   """

    data      = args[0]
    radius    = args[1]
    det_pix   = args[2]
    beam_l    = args[3]
    r_i       = args[4]

    # Compute q-spacing
    qs        = r_i*det_pix/k[1]*4*np.pi/beam_l/2

    bessel = get_bessel(qs,radius)
    theory = k[0]*(bessel**2)

    T      = np.corrcoef(data,theory)[0,1]
    T      = 1/T

    return T

def saxs_fit(k, *args) :
    """Saxs fit by refining I0 and Detector distance
   """

    data      = args[0]
    theory    = args[1]
    det_pix   = args[2]
    beam_l    = args[3]
    r_i       = args[4]
    q_i       = args[5]

    # Compute q-spacing
    q_s       = r_i*det_pix/k*4*np.pi/beam_l/2

    # Interpolate theory to match the experimental q-range
    theory = np.interp(q_s,q_i,theory)

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
    i_f,r_f  = sp.optimize.fmin(bessel_fit, x0=[data[0],r_i],args=(data,qs),disp=0)

    # Optimal solution
    theory   = i_f*(((3*((np.sin(qs*r_f) - (qs*r_f*np.cos(qs*r_f)))/((qs*r_f)**3))))**2);
    score    = np.corrcoef(data,theory)[0,1]

    if plot :

       scale = np.nanmean(data)/np.nanmean(theory)

       fig   = plt.figure(3000)
       plt.semilogy(qs,data,'b-',qs,scale*theory,'r-',linewidth=1)
       plt.xlabel('q (A^-1)')

    return r_f, score


def mc_fit(data, q, R, npart) :
    """Monte carlo based Bessel fit to experimental SAXS data

       @data    Saxs data                                       [nQ]
       @q       Q-vector                                        [nQ]
       @R       Pool of radius sizes
       @npart   Estimated number of particles                   [int]

    """
    # Compute Saxs Intensity for all sizes
    Iall = np.zeros((len(data),len(R)))
    dR   = R[1]-R[0]
    for r in range(len(R)) :
        Btmp      = get_bessel(q,R[r])
        Iall[:,r] = ((R[r]**6)/dR)*(Btmp**2)

    # Initialize Score
    rint  = np.random.randint(len(R))
    I0    = Iall[:,rint]

    C        = np.corrcoef(data,I0)[0,1]
    Dn       = np.zeros(len(R))
    Dn[rint] = 1
    I        = I0

    Dtmp  = np.zeros(len(R))
    # Loop over nr of particles
    for i in range(npart):
        rint        = np.random.randint(len(R))
        Dtmp[rint] += 1

        Itmp        = np.zeros(len(q))
        # Loop over nr of particle bins
        for b in range(len(Dtmp)):
            Itmp += Dtmp[b]*Iall[:,b]

        Ctmp        = np.corrcoef(data,Itmp)[0,1]

        if Ctmp > C :

           C    = Ctmp
           Dn   = Dtmp
           I    = Itmp

        elif Ctmp == C :                # if no added contribution

           C    = Ctmp
           Dn   = Dtmp
           I    = Itmp
           break


    return I,Dn,C



def get_size_distribution(saxs, q, q_i = None, q_f = None, dr = 10, npart = 100, ntrial = 5, plot = 0) :
    """Spherical Besselfit of SAXS data based on a distribution of spheres of different sizes
       Returns optimized distribution of radii r_f in Angstrom

       @saxs    Saxs data                                       [nQ]
       @q       Q-vector                                        [nQ]
       @dr      Size of radius bin                              [float, Ang]
       @npart   Estimated nr of particles                       [integer]
       @q_i     Lower q-bound (Ang-1)                           [float]
       @q_f     Upper q-bound (Ang-1)                           [float]
       @ntrial  Number of MC trials                             [int]
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

    # Determine radius range
    rmin = np.pi/qs[-1]
    rmax = np.pi/qs[0]
    R    = np.arange(rmin,rmax,dr)

    theory = np.zeros((len(qs),ntrial))
    Dn     = np.zeros((len(R),ntrial))
    score  = np.zeros(ntrial)

    # Monte-Carlo refinement
    for t in range(ntrial) :
        theory[:,t],Dn[:,t],score[t] = mc_fit(data,qs,R,npart)

    # Find Optimal solution
    index       = np.unravel_index(score.argmax(), score.shape)
    theory_opt  = theory[:,index]
    n           = int(sum(Dn[:,index]))                 # Nr of refined particles
    Dn_opt      = Dn[:,index]/n                         # Normalize
    score_opt   = score[index]

    return Dn_opt,R,n,score_opt


def get_detdistance(saxs_data, det_0, det_pix, beam_l, saxs_r = None, radius = None, q_theory = None, saxs_theory = None, plot = 0) :
    """Saxs based refinement of the detector distance using knowledge about either:
       1. The size of the calibrant (Spherical calibrants, ex: Polystyrene Latex Spheres)
       2. The theoretical scattering of the particle (Non-spherical calibrants)

       @saxs_data       Saxs data
       @saxs_r          Saxs radius in integer pixels from center to edge of image
       @det_0           Estimated detector distance (mm)
       @det_pix         Pixel size of detector (mm)
       @beam_l          Wavelength of beam (Ang)
       @radius          Radius of calibrant (Ang)
       @q_theory        Q-range of theoretical Saxs data
       @saxs_theory     Theoretical Saxs data
       @plot            Plot fit [0/1]

    """

    if (radius is None) and (saxs_theory is None):
       print("Need to provide calibrant radius or scattering curve")
       return

    if saxs_r is None:
       saxs_r = np.arange(0,len(saxs_data))

    if radius is not None:
       i0_f,det_f  = sp.optimize.fmin(bessel_fit2, x0=[saxs_data[0],det_0],args=(saxs_data,radius,det_pix,beam_l,saxs_r))

       if plot:
          qs       = saxs_r*det_pix/det_f*4*np.pi/beam_l/2
          bessel   = get_bessel(qs,radius)
          theory   = i0_f*(bessel**2)
          score    = np.corrcoef(saxs_data,theory)[0,1]
          scale    = np.nanmean(saxs_data)/np.nanmean(theory)
          fig      = plt.figure(9000)
          plt.semilogy(qs,saxs_data,'b-',qs,scale*theory,'r-',linewidth=1)
          plt.xlabel('q (A^-1)')
          plt.title("Distance= " + str(det_f) + "mm, Score=  " + str(score))

    if saxs_theory is not None:
       det_f       = sp.optimize.fmin(saxs_fit, x0=det_0,args=(saxs_data,saxs_theory,det_pix,beam_l,saxs_r,q_theory))

       if plot:
          qs       = saxs_r*det_pix/det_f*4*np.pi/beam_l/2
          theory   = np.interp(qs,q_theory,saxs_theory)
          score    = np.corrcoef(saxs_data,theory)[0,1]
          scale    = np.nanmean(saxs_data)/np.nanmean(theory)
          fig      = plt.figure(9000)
          plt.semilogy(qs,saxs_data,'b-',qs,scale*theory,'r-',linewidth=1)
          plt.xlabel('q (A^-1)')
          plt.title("Distance= " + str(det_f) + "mm, Score=  " + str(score))


    return det_f


def get_symmetry(c2, phi_lim = None, plot = 0) :
    """Compute pi/2 symmetry for s
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

    return sym


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
    # Set q-limits
    if q_lim is None:
       q_lim    = [0,len(q)]

    # Set phi-limits
    if phi_lim is None:
       phi_lim   =[0,len(phi)]

    B       = np.zeros((l_max,len(q)))
    P       = np.zeros((len(phi),l_max))
    c2recon = np.zeros(c2.shape)

    # Loop over all q
    for qq in range(q_lim[0],q_lim[1]) :
        # Compute theta
        theta = np.pi/2 - np.arcsin((q[qq]*lambd)/(4*np.pi))
        x     = np.cos(theta)**2 + np.sin(theta)**2 * np.cos(phi[phi_lim[0]:phi_lim[1]])
#        x     = min(x,1)
#        x     = max(x,-1)
        y     = 0
        # Loop over all l
        for ll in range(0,l_max):
            Plm    = np.polynomial.legendre.legval(x,ll)
            P[:,y] = (1/(4*np.pi))*Plm
            y     += 1
        # Use SVD for lsq-fitting to c2
        [U,S,V]       = np.linalg.svd(P)
        S             = np.eye(len(S),len(S))*S
        from IPython import embed; embed()
        B[:,qq]       = V*np.linalg.inv(S)*np.transpose(U[phi_lim[0]:phi_lim[1],:])*c2[phi_lim[0]:phi_lim[1],qq]

        c2recon[:,qq] = np.transpose(np.transpose(B[:,qq]) * np.transpose((U*S*np.transpose(V))))

    return B , c2recon

   ### NOT FINISHED  ###




def get_h5_event(file_stamps) :
    """Generate h5 timestamps from time-stamps

       @file_stamps  Parameter file containing  Time-stamp and h5-path to image
                      ex: /2012-07-29T01:01:13_23690/FrontPnCCDIsHit/HistData

    """

    timestamps = []
    filenr     = []

    for i in range(len(file_stamps)) :

        temp = file_stamps['f0'][i].split('/')
        timestamps.append(temp[1])
        filenr.append(file_stamps['f1'][i])

    return timestamps,filenr


def get_psana_time(evt) :
    """Generate psana timestamp from psana event (only works for xtc structures)

       @evt  Psana event

    """

    time     = np.zeros((3,))

    id       = evt.get(EventId)
    time[0]  = id.time()[0]     # Seconds
    time[1]  = id.time()[1]     # Nanon seconds
    time[2]  = id.fiducials()   # Fiducial

    return time

def get_time(file_stamps) :
    """Generate list of timestamps from parameter files (only works for xtc structures)

       @file_stamps  Time-stamp list

    """

    times    = []

    for i in range(len(file_stamps)) :
        t       = file_stamps[i]
        time    = np.zeros((3,))
        time[0] = t[1]
        time[1] = t[2]
        time[2] = t[3]
        times.append(time)

    return times



def get_psana_event(file_stamps) :
    """Generate psana timestamps from time-stamps

       @file_stamps  Parameter file containing  Time, Seconds, Nanoseconds,Fiducials

    """


    timestamps = []

    for i in range(len(file_stamps)) :


        sec   = int(file_stamps[i][1])
        nsec  = int(file_stamps[i][2])
        fid   = int(file_stamps[i][3])

        et   = EventTime(int((sec<<32)|nsec),fid)

        timestamps.append(et)

    return timestamps


def get_psana_event2(file_stamp,run = None) :
    """Generate psana timestamp from time-stamps

       @file_stamps  Parameter file containing  Time, Seconds, Nanoseconds,Fiducials

    """

    sec         = int(file_stamp[1])
    nsec        = int(file_stamp[2])
    fid         = int(file_stamp[3])

    et          = EventTime(int((sec<<32)|nsec),fid)
    timestamp   = run.event(et)


    return timestamp
