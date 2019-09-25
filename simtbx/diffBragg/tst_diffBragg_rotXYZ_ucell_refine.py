from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
args = parser.parse_args()

from dxtbx.model.crystal import Crystal
from cctbx import uctbx
from scitbx.matrix import sqr, rec, col
import numpy as np
from scipy.spatial.transform import Rotation
import pylab as plt

from simtbx.diffBragg.nanoBragg_crystal import nanoBragg_crystal
from simtbx.diffBragg.sim_data2 import SimData
from simtbx.diffBragg import utils
from simtbx.diffBragg.refiners import RefineMissetAndUcell
from simtbx.diffBragg.refiners.crystal_systems import MonoclinicManager

ucell = (55, 65, 75, 90, 95, 90)
ucell2 = (55.2, 66, 74, 90, 94.3, 90)
symbol = "P121"

# generate a random raotation
rotation = Rotation.random(num=1, random_state=100)[0]
Q = rec(rotation.as_quat(), n=(4, 1))
rot_ang, rot_axis = Q.unit_quaternion_as_axis_and_angle()

# generate a small perturbation rotation
np.random.seed(1)
perturb_rot_axis = np.random.random(3)
perturb_rot_axis /= np.linalg.norm(perturb_rot_axis)
perturb_rot_ang = 0.1  # 0.1 degree random perturbtation

# make the ground truth crystal:
a_real, b_real, c_real = sqr(uctbx.unit_cell(ucell).orthogonalization_matrix()).transpose().as_list_of_lists()
C = Crystal(a_real, b_real, c_real, symbol)
C.rotate_around_origin(rot_axis, rot_ang)

a2_real, b2_real, c2_real = sqr(uctbx.unit_cell(ucell2).orthogonalization_matrix()).transpose().as_list_of_lists()
C2 = Crystal(a2_real, b2_real, c2_real, symbol)
C2.rotate_around_origin(rot_axis, rot_ang)
assert np.allclose(C2.get_U(), C.get_U())
C2.rotate_around_origin(col(perturb_rot_axis), perturb_rot_ang)

# Setup the simulation and create a realistic image
# with background and noise
# <><><><><><><><><><><><><><><><><><><><><><><><><>
nbcryst = nanoBragg_crystal()
nbcryst.dxtbx_crystal = C   # simulate ground truth
nbcryst.thick_mm = 0.1
nbcryst.Ncells_abc = 12, 12, 12

SIM = SimData()
SIM.detector = SimData.simple_detector(150, 0.1, (512, 512))
SIM.crystal = nbcryst
SIM.instantiate_diffBragg(oversample=0)
SIM.D.progress_meter = False
SIM.water_path_mm = 0.005
SIM.air_path_mm = 0.1
SIM.add_air = True
SIM.add_Water = True
SIM.include_noise = True
SIM.D.add_diffBragg_spots()
spots = SIM.D.raw_pixels.as_numpy_array()
SIM._add_background()
SIM._add_noise()
# This is the ground truth image:
img = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels *= 0

# Simulate the perturbed image for comparison
SIM.D.Bmatrix = C2.get_B()
SIM.D.Umatrix = C2.get_U()
SIM.D.add_diffBragg_spots()
SIM._add_background()
SIM._add_noise()
# Perturbed image:
img_pet = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels *= 0

# NOTE: NEED TO DO SPOT FINDING AND WHAT NOT
# spot_rois, abc_init, img, SimData_instance,
# locate the strong spots and fit background planes
# <><><><><><><><><><><><><><><><><><><><><><><><><>
spot_data = utils.get_spot_data(spots, thresh=19)
ss_spot, fs_spot = map(np.array, zip(*spot_data["maxIpos"]))  # slow/fast  scan coords of strong spots
num_spots = len(ss_spot)
if args.plot:
    plt.imshow(img, vmax=200)
    plt.plot(fs_spot, ss_spot, 'o', mfc='none', mec='r')
    plt.title("Simulated image with strong spots marked")
    plt.show()

is_bg_pixel = np.ones(img.shape, bool)
for bb_ss, bb_fs in spot_data["bboxes"]:
    is_bg_pixel[bb_ss, bb_fs] = False

# now fit tilting planes
shoebox_sz = 20
tilt_abc = np.zeros((num_spots, 3))
spot_roi = np.zeros((num_spots, 4), int)
if args.plot:
    patches = []
img_shape = SIM.detector[0].get_image_size()
for i_spot, (x_com, y_com) in enumerate(zip(fs_spot, ss_spot)):
    i1 = int(max(x_com - shoebox_sz / 2., 0))
    i2 = int(min(x_com + shoebox_sz / 2., img_shape[0]-1))
    j1 = int(max(y_com - shoebox_sz / 2., 0))
    j2 = int(min(y_com + shoebox_sz / 2., img_shape[1]-1))

    shoebox_img = img[j1:j2, i1:i2]
    shoebox_mask = is_bg_pixel[j1:j2, i1:i2]

    tilt, bgmask, coeff = utils.tilting_plane(
        shoebox_img,
        mask=shoebox_mask,  # mask specifies which spots are bg pixels...
        zscore=2)

    tilt_abc[i_spot] = coeff[1], coeff[2], coeff[0]  # store as fast-scan coeff, slow-scan coeff, offset coeff

    spot_roi[i_spot] = i1, i2, j1, j2
    if args.plot:
        R = plt.Rectangle(xy=(x_com-shoebox_sz/2, y_com-shoebox_sz/2.),
                      width=shoebox_sz,
                      height=shoebox_sz,
                      fc='none', ec='r')
        patches.append(R)

if args.plot:
    patch_coll = plt.mpl.collections.PatchCollection(patches, match_original=True)
    plt.imshow(img, vmin=0, vmax=200)
    plt.gca().add_collection(patch_coll)
    plt.show()

UcellMan = MonoclinicManager(
    a=ucell2[0],
    b=ucell2[1],
    c=ucell2[2],
    beta=ucell2[4]*np.pi/180.)

nbcryst.dxtbx_crystal = C2
SIM.crystal = nbcryst

init_Umat_norm = np.abs(np.array(C2.get_U()) - np.array(C.get_U())).sum()
init_Bmat_norm = np.abs(np.array(C2.get_B()) - np.array(C.get_B())).sum()

RUC = RefineMissetAndUcell(
    spot_rois=spot_roi,
    abc_init=tilt_abc,
    img=img,
    SimData_instance=SIM,
    plot_images=args.plot,
    ucell_manager=UcellMan)
RUC.trad_conv = True
RUC.trad_conv_eps = 1e-5
RUC.run()

ang, ax = RUC.get_correction_misset(as_axis_angle_deg=True)
C2.rotate_around_origin(ax,ang)
C2.set_B(RUC.get_refined_Bmatrix())

final_Umat_norm = np.abs(np.array(C2.get_U()) - np.array(C.get_U())).sum()
final_Bmat_norm = np.abs(np.array(C2.get_B()) - np.array(C.get_B())).sum()

print("Results!")
print("Before refinement: Umatrix distance=%2.7g, Bmatrix distance=%2.7g" % (init_Umat_norm, init_Bmat_norm))
print("After refinement: Umatrix distance=%2.7g, Bmatrix distance=%2.7g" % (final_Umat_norm, final_Bmat_norm))
print("")
print("ground truth unit cell: %2.7g,%2.7g,%2.7g,%2.7g,%2.7g,%2.7g" % ucell)
print("unit cell passed to refinement: %2.7g,%2.7g,%2.7g,%2.7g,%2.7g,%2.7g" % ucell2)
print("refined unit cell: %2.7g,%2.7g,%2.7g,%2.7g,%2.7g,%2.7g" % C2.get_unit_cell().parameters())
print("")
print("Perturbation axis =%+2.7g,%+2.7g,%+2.7g and angle=%+2.7g deg"
      % (perturb_rot_axis[0], perturb_rot_axis[1], perturb_rot_axis[2], perturb_rot_ang))
print("Misset applied during refinement: axis=%+2.7g,%+2.7g,%+2.7g and angle=%+2.7g deg"
      % (ax[0], ax[1], ax[2], ang))

print("OK")
