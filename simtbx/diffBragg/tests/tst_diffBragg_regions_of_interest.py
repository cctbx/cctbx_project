from __future__ import division

def main():

    import numpy as np
    from simtbx.diffBragg.utils import get_diffBragg_instance

    D = get_diffBragg_instance()
    D.nopolar = False
    D.interpolate = 0
    D.spot_scale = 1e5
    rois = (20, 100, 30, 100), (10, 90, 40, 80)
    #px = 39,52
    for x1, x2, y1, y2 in rois:
        D.raw_pixels *= 0
        D.region_of_interest = ((x1, x2), (y1, y2))
        #D.printout_pixel_fastslow = x1+px[0], y1+px[1]

        D.add_nanoBragg_spots()
        nano_pixels = D.raw_pixels.as_numpy_array()
        nano_roi_pixels = nano_pixels[y1:y2, x1:x2]

        D.raw_pixels *= 0
        D.vectorize_umats()
        D.add_diffBragg_spots()
        npix_sim = (x2-x1)*(y2-y1)
        diff_roi_pixels = D.raw_pixels_roi.as_numpy_array()[:npix_sim].reshape((y2-y1, x2-x1))

        assert np.allclose(diff_roi_pixels, nano_roi_pixels, atol=1e-9)


if __name__ == "__main__":
    main()
    print ("OK!")
