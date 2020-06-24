
from simtbx.nanoBragg import shapetype
from simtbx.nanoBragg import nanoBragg
from dxtbx.model.crystal import CrystalFactory
from scitbx.matrix import sqr

import numpy as np

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# This test simulates all the shape models on both GPU and CPU
# and verifies that the results are the same
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

def main(shape=shapetype.Tophat, cuda=False, seed=None):

    SIM = nanoBragg(verbose=10, oversample=0)
   
    # Defaults 
    cr = {'__id__': 'crystal',
          'real_space_a': (200, 0, 0),
          'real_space_b': (0, 180, 0),
          'real_space_c': (0, 0, 150),
          'space_group_hall_symbol': '-P 4 2'}
    cryst = CrystalFactory.from_dict(cr)
    SIM.Amatrix = sqr(cryst.get_A()).transpose().elems
        
    SIM.detpixels_fastslow = (64,64)
    SIM.pixel_size_mm=0.1
    SIM.wavelength_A=1.2
    SIM.verbose=10
    SIM.flux=1e12
    SIM.mosaic_spread_deg=0.02
    SIM.mosaic_domains=10
    SIM.polarization=1
    SIM.distance_mm=100
    SIM.F000=3e3  # FIXME this has to be equivalent to default_F, or else set in Fhkl, otherwise the test fails
    SIM.default_F=3e3
    SIM.progress_meter=True
    SIM.beamsize_mm=0.005
    SIM.exposure_s=1
    SIM.Ncells_abc=(15,15,15)
    SIM.show_params()

    # variable 
    SIM.xtal_shape=shape
    SIM.show_params() 

    if seed is not None:
        SIM.seed = seed
        SIM.randomize_orientation()
    if cuda:
        # NOTE: uncomment the following 4 lines and comment the add_nanoBragg_spots_cuda() call
        # in order to use current dev-mode code!
        #SIM.allocate_cuda()
        #SIM.add_nanoBragg_spots_cuda_update()
        #SIM.get_raw_pixels_cuda()
        #SIM.deallocate_cuda()
        
        SIM.add_nanoBragg_spots_cuda()
    else:
        SIM.add_nanoBragg_spots()
    img = SIM.raw_pixels.as_numpy_array()
    return img

if __name__ == "__main__":
    failures = 0
    shapes = shapetype.Tophat, shapetype.Gauss, shapetype.Square, shapetype.Round
    failed_shapes = []
    for shape in shapes:
        img_cuda = main(cuda=True)
        img = main(cuda=False)
        #np.savez("_shapes", img_cuda, img)
        
        # from the docs for np.allclose:
        #       np.allclose(a,b,rtol,atol)
        #       If the following equation is element-wise True, then allclose returns
        #       True.
        #            absolute(`a` - `b`) <= (`atol` + `rtol` * absolute(`b`))


        # check whether all values are within 0.1 photons:
        if not np.allclose(img, img_cuda, rtol=0, atol=0.1):
            failed_shapes.append( repr(shape))
    if failed_shapes:
        print ("\nThe following shape models failed the test:")
        for shape in failed_shapes:
            print ("  %s failed the test miserably" % shape)
        assert False, "This test has failed"
    else:
        print("OK")

