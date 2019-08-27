
def main():

    from dxtbx.model.crystal import CrystalFactory
    from dxtbx.model.detector import DetectorFactory
    from dxtbx.model.beam import BeamFactory
    from simtbx.nanoBragg.tst_nanoBragg_basic import fcalc_from_pdb
    from simtbx.nanoBragg import shapetype
    from simtbx.diffBragg import diffBragg
    import numpy as np

    n_trials = 10
    wavelen = 1.24
    flux = 1e12
    SHAPE = shapetype.Gauss
    np.random.seed( n_trials)
    angles_XYZ = np.random.random( (n_trials, 3))
    NCELLS_ABC=(15,15,15)
    
    beam_descr = {'direction': (0.0, 0.0, 1.0),
                 'divergence': 0.0,
                 'flux': 5e11,
                 'polarization_fraction': 1.,
                 'polarization_normal': (0.0, 1.0, 0.0),
                 'sigma_divergence': 0.0,
                 'transmission': 1.0,
                 'wavelength': wavelen}

    cryst_descr = {'__id__': 'crystal',
                   'real_space_a': (79, 0, 0),
                   'real_space_b': (0, 79, 0),
                   'real_space_c': (0, 0, 38),
                   'space_group_hall_symbol': '-P 4 2'}

    det_descr = {'panels':
                   [{'fast_axis': (-1.0, 0.0, 0.0),
                     'gain': 1.0,
                     'identifier': '',
                     'image_size': (196, 196),
                     'mask': [],
                     'material': '',
                     'mu': 0.0,
                     'name': 'Panel',
                     'origin': (19.6, -19.6, -550),
                     'pedestal': 0.0,
                     'pixel_size': (0.1, 0.1),
                     'px_mm_strategy': {'type': 'SimplePxMmStrategy'},
                     'raw_image_offset': (0, 0),
                     'slow_axis': (0.0, 1.0, 0.0),
                     'thickness': 0.0,
                     'trusted_range': (0.0, 65536.0),
                     'type': ''}]}

    DET = DetectorFactory.from_dict(det_descr)
    BEAM = BeamFactory.from_dict(beam_descr)
    crystal = CrystalFactory.from_dict(cryst_descr)

    Fhkl = fcalc_from_pdb(resolution=4, algorithm="fft", wavelength=wavelen)

    D = diffBragg(DET, BEAM, verbose=0, panel_id=0)
    D.xtal_shape = SHAPE
    D.Ncells_abc = NCELLS_ABC
    D.wavelength_A = wavelen
    D.flux = flux
    D.mosaic_spread_deg = 0.01
    D.mosaic_domains = 10
    D.Fhkl = Fhkl

    rotX,rotY,rotZ = 0,1,2
    D.refine(rotX)  # rotX
    D.refine(rotY)  # rotY
    D.refine(rotZ)  # rotZ

    D.initialize_managers()
    D.vectorize_umats()

    from scitbx.matrix import col, sqr
    x = col((1,0,0))
    y = col((0,-1,0))
    z = col((0,0,1))

    for i_ang, (thetaX, thetaY, thetaZ) in enumerate(angles_XYZ):

        RX = x.axis_and_angle_as_r3_rotation_matrix(thetaX, deg=True)
        RY = y.axis_and_angle_as_r3_rotation_matrix(thetaY, deg=True)
        RZ = z.axis_and_angle_as_r3_rotation_matrix(thetaZ, deg=True)

        Arecip_orig = sqr(crystal.get_A())
        Areal = Arecip_orig.inverse()
        Areal = RX*RY*RZ*Areal
        Arecip = Areal.inverse()

        D.thetaX=0
        D.thetaY=0
        D.thetaZ=0
        D.Amatrix = Arecip.transpose().elems
        D.add_diffBragg_spots()
        imgA = D.raw_pixels.as_numpy_array()

        D.raw_pixels *= 0
        D.set_value(rotX, thetaX)
        D.set_value(rotY, thetaY)
        D.set_value(rotZ, thetaZ)
        D.Amatrix = Arecip_orig.transpose().elems
        D.add_diffBragg_spots()
        imgB = D.raw_pixels.as_numpy_array()

        assert(np.allclose(imgA, imgB, atol=1e-4))
        print("OK (%d / %d)" % (i_ang+1, len(angles_XYZ)))


if __name__=="__main__":
    main()
    print("OK")
