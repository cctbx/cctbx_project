
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
    process_simdata(plot=False, angles=angles)
# NOTE: misset is an RX*RY*RZ matrix where the angles of rotation are .4,.3,.2 in degrees respectively
S = SimData()
S.crystal.missetting_matrix = misset
S.instantiate_diffBragg()
# Step 1:
RR = RefineRot(spot_rois=spot_rois,
               abc_init=abc_init, img=img, SimData_instance=S)
RR.run()
angles_refined = np.array(RR.x[-4:-1]) * 180 / np.pi
assert np.all(np.round(angles_refined, 1) == -angles)

# Step 2:
crystal = RR.S.crystal
C0 = deepcopy(crystal.dxtbx_crystal)  # NOTE: ground truth

# rotate the ground truth crystal with the .4,.3,.2 degree perturbations
C = deepcopy(C0)
q = misset.r3_rotation_matrix_as_unit_quaternion()
rot_ang, rot_ax = q.unit_quaternion_as_axis_and_angle(deg=True)
C.rotate_around_origin(rot_ax, rot_ang)
initial_ang_offset = compare_orientation_matrices.difference_rotation_matrix_axis_angle(C, C0)[2]

# get the correction missetting from refiner
correction_ang, correction_ax = RR.get_correction_misset(as_axis_angle_deg=True)
C.rotate_around_origin(correction_ax, correction_ang)
final_ang_offset = compare_orientation_matrices.difference_rotation_matrix_axis_angle(C, C0)[2]

print("Initial missorientation=%1.2g degrees" % initial_ang_offset)
print("Final missorientation=%1.2g degrees" % final_ang_offset)

assert initial_ang_offset > 0.5
assert final_ang_offset < 0.005
print("OK!")

