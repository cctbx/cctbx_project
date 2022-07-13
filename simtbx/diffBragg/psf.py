from __future__ import division

import numpy as np
from dials.algorithms.image.filter import convolve
import math
from scipy.signal import convolve2d
from dials.array_family import flex



def fiber2D_integ(x,y,g):
    return math.atan((x*y)/(g*math.sqrt(g*g + x*x + y*y)))/(2.0*math.pi)


def makeMoffat_integPSF(fwhm_pixel, sizeX, sizeY):
    ''' Integral form of contribution of moffat PSF to signal recorded in a square pixel'''
    g = fwhm_pixel*0.65238
    psf = np.zeros((sizeX, sizeY))
    sx = int(sizeX/2)
    sy = int(sizeY/2)
    for y in range(-sy, sy+1):
    #for y in range(sy, -(sy+1), -1):
      for x in range(-sx, sx+1):
        #for x in range(-sx, sx+1):
        # Holton 2012 paper says x,y should be pixel center; this does not seem right ?
        psf[x+sx,y+sy] = fiber2D_integ(x+1./2,y+1./2,g)-fiber2D_integ(x+1./2,y-1./2,g)-fiber2D_integ(x-1./2,y+1./2,g)+fiber2D_integ(x-1./2,y-1./2,g)
        #psf[x+sx, -y+sy] = fiber2D_integ(x+1./2,y+1./2,g)-fiber2D_integ(x+1./2,y-1./2,g)-fiber2D_integ(x-1./2,y+1./2,g)+fiber2D_integ(x-1./2,y-1./2,g)
        # Trying to get pixel center instead
        #psf[x+sx, -y+sy] = fiber2D_integ(x+1,y+1,g)-fiber2D_integ(x+1,y,g)-fiber2D_integ(x,y+1,g)+fiber2D_integ(x,y,g)
    psf = psf/psf.sum()
    psf = psf.tolist()
    psf = flex.double(psf)
    return psf


def convolve_padded_img(img, psf, sz=5):
    img = np.array(img)
    iY, iX = img.shape
    pY, pX = psf.focus()

    new_iY = iY
    if pY >= iY - sz:
        new_iY = pY + sz
    new_iX = iX
    if pX >= iX - sz:
        new_iX = pX + sz

    assert new_iX >= iX
    assert new_iY >= iY
    padX = new_iX - iX
    padY = new_iY - iY

    x = int(padX/2)
    y = int(padY/2)
    img = np.pad(img, ((y, y+1), (x, x+1)), mode='median')
    assert img.shape[0] >= pY + sz
    assert img.shape[1] >= pX + sz

    conv_img = convolve(flex.double(img), psf)
    conv_img = conv_img.as_numpy_array()[y:y+iY, x:x+iX]
    return conv_img


def convolve_with_psf(image_data, fwhm=27.0, pixel_size=177.8, psf_radius=7, sz=5, psf=None, use_scipy=True):
    ''' Given a 2D numpy array of image data, convolve with a PSF. '''
    # Currently only supporting fiber PSF i.e power law form as proposed in Holton et. al 2012, Journal of Synchotron Radiation
    if psf is None:
        xpsf=2*psf_radius+1
        ypsf=2*psf_radius+1
        fwhm_pixel=fwhm/pixel_size
        psf = makeMoffat_integPSF(fwhm_pixel, xpsf, ypsf)
    if use_scipy:
        psf_img = psf.as_numpy_array()
        med_img = (image_data[0,0] + image_data[0,-1] + image_data[-1,0] + image_data[-1,-1])*0.25 #np.median(image_data)
        convolved_image = convolve2d(image_data, psf_img, mode='same', fillvalue=med_img)

    else:
        img_shape = image_data.shape
        psf_shape = psf.focus()
        if psf_shape[0] > img_shape[0] - sz or psf_shape[1] > img_shape[1] - sz:
            convolved_image = convolve_padded_img(image_data, psf, sz)
        else:
            convolved_image = convolve(flex.double(image_data), psf)
            convolved_image = convolved_image.as_numpy_array()
    return convolved_image


if __name__=="__main__":
    psf_args = {'fwhm': 100, 'pixel_size': 80., 'psf_radius': 7}
    fwhm_pix = psf_args["fwhm"] / psf_args["pixel_size"]
    kern_size = psf_args["psf_radius"]*2 + 1
    PSF = makeMoffat_integPSF(fwhm_pix, kern_size, kern_size)

    a = np.random.random((100,100))
    a2 = convolve_with_psf(a, psf=PSF, **psf_args)
    from pylab import plt
    plt.subplot(121)
    plt.gca().set_title("without PSF")
    plt.imshow(a)
    plt.subplot(122)
    plt.gca().set_title("with PSF")
    plt.imshow(a2)
    plt.show()

