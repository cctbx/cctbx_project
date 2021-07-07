from __future__ import division


def main():
    from simtbx.nanoBragg.sim_data import SimData

    S = SimData(use_default_crystal=True)
    S.instantiate_diffBragg()
    S.D.nopolar = True
    S.D.oversample = 3

    S.D.add_diffBragg_spots()
    S._add_background()

    diff_img = S.D.raw_pixels.as_numpy_array()

    S.D.raw_pixels *= 0
    S.D.add_nanoBragg_spots()
    S._add_background()

    nano_img = S.D.raw_pixels.as_numpy_array()

    import numpy as np

    assert np.allclose(diff_img, nano_img, atol=1e-9)


if __name__ == "__main__":
    main()
    print ("OK!")
