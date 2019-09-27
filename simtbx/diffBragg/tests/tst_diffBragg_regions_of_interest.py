
def main():

    import numpy as np
    from simtbx.diffBragg.utils import get_diffBragg_instance

    D = get_diffBragg_instance()
    D.spot_scale = 1e5
    rois = (20, 100, 30, 100), (10, 90, 40, 80)
    for x1, x2, y1, y2 in rois:
        D.raw_pixels *= 0
        D.region_of_interest = ((x1, x2), (y1, y2))
        D.add_nanoBragg_spots()
        nano_pixels = D.raw_pixels.as_numpy_array()
        nano_roi_pixels = nano_pixels[y1:y2+1, x1:x2+1]

        D.raw_pixels *= 0
        D.region_of_interest = ((x1, x2), (y1, y2))
        D.vectorize_umats()
        D.add_diffBragg_spots()
        diff_roi_pixels = D.raw_pixels_roi.as_numpy_array()

        assert np.allclose(diff_roi_pixels, nano_roi_pixels, atol=1e-9)


if __name__ == "__main__":
    main()
    print ("OK!")


