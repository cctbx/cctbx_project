from __future__ import division

from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
from dxtbx.model import Crystal
from simtbx.diffBragg.utils import  fcalc_from_pdb

Fcomplex = fcalc_from_pdb(resolution=2, wavelength=1, algorithm='fft', ucell=(79,79,38.90,90,90), symbol="P43212",as_amplitudes=False)
Famp = fcalc_from_pdb(resolution=2, wavelength=1, algorithm='fft', ucell=(79,79,38.90,90,90), symbol="P43212",as_amplitudes=True)


def sim_spots(F):

    a_real = (79,0,0)
    b_real = (0,79,0)
    c_real = (0,0,38)
    C = Crystal(a_real, b_real, c_real, 'P43212')

    nbcryst = NBcrystal(init_defaults=True)
    nbcryst.dxtbx_crystal = C  # simulate ground truth
    nbcryst.thick_mm = 0.1
    nbcryst.Ncells_abc = 10,10,10

    nbcryst.miller_array = F
    print("Ground truth ncells = %f" % (nbcryst.Ncells_abc[0]))

    # ground truth detector
    from simtbx.nanoBragg.sim_data import SimData
    DET_gt = SimData.simple_detector(150, 0.177, (600, 600))

    # initialize the simulator
    SIM = SimData(use_default_crystal=True)
    SIM.detector = DET_gt
    SIM.crystal = nbcryst
    SIM.instantiate_diffBragg(oversample=0)
    SIM.D.default_F = 0
    SIM.D.progress_meter = False
    SIM.D.add_diffBragg_spots()
    SIM.D.F000 = 0
    SPOTS = SIM.D.raw_pixels.as_numpy_array()
    SIM.D.free_all()
    SIM.D.free_Fhkl2()
    return SPOTS

A = sim_spots(Fcomplex)
B = sim_spots(Famp)
import numpy as np
assert np.allclose(A,B)
print("OK!")
