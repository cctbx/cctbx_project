
import numpy as np
from simtbx.diffBragg.sim_data2 import  SimData
from simtbx.diffBragg.process_simdata import process_simdata
from simtbx.diffBragg.refiners import RefineRot
from dials.algorithms.indexing import compare_orientation_matrices
from copy import deepcopy
from simtbx.diffBragg import utils

angles = np.array([.4, .3, .2])
# NOTE: make this more obvious (whats going on in process simdata:
spot_rois, spot_hkl, Amat_init, Amat_known, abc_init, img, misset = \
    process_simdata(plot=True, angles=angles)
# NOTE: misset is an RX*RY*RZ matrix where the angles of rotation are .4,.3,.2 in degrees respectively
S = SimData()
S.crystal.missetting_matrix = misset

# Step 1:
RR = RefineRot(spot_rois=spot_rois,
               abc_init=abc_init, img=img, SimData_instance=S)
angles_refined = np.array(RR.x[-4:-1]) * 180 / np.pi
assert np.all(np.round(angles_refined, 1) == -angles)

# Step 2:
crystal = RR.S.crystal
C0 = deepcopy(crystal.dxtbx_crystal)  # NOTE: ground truth
C = crystal.dxtbx_crystal_with_missetting()
#C.set_A(sqr(crystal.Amatrix_realspace).inverse())
angles_rad = RR.x[-4:-1]
C2 = utils.refine_model_from_angles(C, angles=angles_rad)
ang = compare_orientation_matrices.difference_rotation_matrix_axis_angle(C, C0)[2]
ang2 = compare_orientation_matrices.difference_rotation_matrix_axis_angle(C2, C0)[2]
assert ang > 0.5
assert ang2 < 0.002
print("OK!")

