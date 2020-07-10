#!/usr/bin/env python
# -*- coding: utf-8 -*-



from timemory.component import WallClock


from simtbx.nanoBragg import nanoBragg, shapetype
import numpy as np
from dxtbx.model.crystal import CrystalFactory
from scitbx.matrix import sqr






def main(shape=shapetype.Tophat, cuda=False, seed=None, new_cuda=False):

    wc = WallClock("nanoBragg:py")

    wc.push()
    wc.start()
    SIM = nanoBragg(detpixels_slowfast=(2048,2048), pixel_size_mm=0.05,
                    verbose=10, oversample=0)
    wc.stop()
    wc.pop()

    # Defaults 
    cr = {'__id__': 'crystal',
          'real_space_a': (200, 0, 0),
          'real_space_b': (0, 180, 0),
          'real_space_c': (0, 0, 150),
          'space_group_hall_symbol': '-P 4 2'}
    cryst = CrystalFactory.from_dict(cr)
    SIM.Amatrix = sqr(cryst.get_A()).transpose().elems
    SIM.wavelength_A=1.2
    SIM.verbose=10
    SIM.flux=1e12
    SIM.mosaic_spread_deg=0.02
    SIM.mosaic_domains=10
    SIM.polarization=1
    SIM.distance_mm=700
    SIM.F000=3e3  # FIXME this has to be equivalent to default_F, or else set in Fhkl, otherwise the test fails
    SIM.default_F=3e3
    SIM.progress_meter=True
    SIM.beamsize_mm=0.005
    SIM.exposure_s=1
    SIM.Ncells_abc=(15,15,15)
    SIM.show_params()
    # variable 
    SIM.xtal_shape=shape
    # boost up the signal (number of mosaic blocks in crystal)
    SIM.spot_scale = 1e3
    SIM.show_params() 
    if seed is not None:
        SIM.seed = seed
        SIM.randomize_orientation()
    if cuda:
        # NOTE: uncomment the following 4 lines and comment the add_nanoBragg_spots_cuda() call
        # in order to use current dev-mode code!
        if new_cuda:
            # TODO: `allocate_cuda` doesn't seem to be defined anywhere
            SIM.allocate_cuda()
            SIM.add_nanoBragg_spots_cuda_update()
            SIM.get_raw_pixels_cuda()
            SIM.deallocate_cuda()
        else:
            SIM.add_nanoBragg_spots_cuda_nvtx()
    else:
        SIM.add_nanoBragg_spots_nvtx()
    SIM.add_noise_nvtx()
    #SIM.raw_pixels += 
    SIM.to_smv_format_nvtx(fileout="intimage_001.img")
    img = SIM.raw_pixels.as_numpy_array()
    return img


# img = main(cuda=True, new_cuda=False)
# Call twice -- to give cudaMalloc a chance
# img = main(cuda=True, new_cuda=False)
# TODO: new_cuda is broken
#img = main(cuda=True, new_cuda=True)
img = main(cuda=False)

#import pylab as plt
#plt.imshow( img)
#plt.show()


