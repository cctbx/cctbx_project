
import numpy as np
from simtbx.diffBragg.sim_data import SimData
from simtbx.diffBragg.process_simdata import process_simdata
from simtbx.diffBragg.refiners import RefineRot
from copy import deepcopy

angles = np.array([.4, .3, .2])
# NOTE: make this more obvious (whats going on in process simdata:
spot_rois, spot_hkl, Amat_init, Amat_known, abc_init, img, misset = \
    process_simdata(plot=False, angles=angles)
# NOTE: misset is an RX*RY*RZ matrix where the angles of rotation are .4,.3,.2 in degrees respectively
S = SimData()
ang, ax = misset.r3_rotation_matrix_as_unit_quaternion().unit_quaternion_as_axis_and_angle(deg=True)
C0 = deepcopy(S.crystal.dxtbx_crystal)  # ground truth
S.crystal.dxtbx_crystal.rotate_around_origin(ax, ang)

S.instantiate_diffBragg()
# Step 1:
RR = RefineRot(spot_rois=spot_rois,
               abc_init=abc_init, img=img, SimData_instance=S)
RR.run()
angles_refined = np.array(RR.x[-4:-1]) * 180 / np.pi

# most important test result here is that angles_refined is approx -.4,-.3,-.2
assert np.all(np.round(angles_refined, 1) == -angles)

# Step 2:
crystal = RR.S.crystal

# rotate the ground truth crystal with the .4,.3,.2 degree perturbations
C = S.crystal.dxtbx_crystal
initial_miss = np.diff(zip(C.get_U(), C0.get_U()))[:, 0].sum()

# get the correction missetting from refiner
correction_ang, correction_ax = RR.get_correction_misset(as_axis_angle_deg=True)
C.rotate_around_origin(correction_ax, correction_ang)
final_miss = np.diff(zip(C.get_U(), C0.get_U()))[:, 0].sum()
assert final_miss < 5e-2*initial_miss   # this is volatile criterion, no worries if it starts failing
print("OK!")
