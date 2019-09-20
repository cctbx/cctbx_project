from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
args = parser.parse_args()

from IPython import embed

from cctbx import sgtbx
from rstbx.symmetry.constraints import parameter_reduction
import numpy as np
from scipy.spatial.transform import Rotation

from scitbx.matrix import sqr
from scitbx.array_family import flex
from scitbx.matrix import rec
from simtbx.diffBragg.nanoBragg_crystal import nanoBragg_crystal
from cctbx import uctbx
from simtbx.diffBragg.sim_data2 import SimData
from dxtbx.model.crystal import Crystal
from itertools import cycle
import pylab as plt



def toggle_image(imgs, titles, **kwargs):
    imgs = cycle(imgs)
    titles = cycle(titles)
    icount = 0
    while icount < 6:
        plt.cla()
        plt.imshow(imgs.next(), **kwargs)
        plt.title(titles.next() + "\nContinuing in %d" % (6 - icount))
        plt.draw()
        plt.pause(.7)
        icount += 1
        plt.close()

ucell = (55, 66, 77, 90, 95, 90)
symbol = "P1211"

#a_real, b_real, c_real = sqr(uctbx.unit_cell(ucell).orthogonalization_matrix()).transpose().as_list_of_lists()
a_real, b_real, c_real = np.reshape(uctbx.unit_cell(ucell).orthogonalization_matrix(), (3,3)).T   #transpose().as_list_of_lists()
C = Crystal(a_real, b_real, c_real, symbol)
point_group = sgtbx.space_group_info(symbol=symbol).group().build_derived_point_group()
S = parameter_reduction.symmetrize_reduce_enlarge(point_group)

#rotation = Rotation.random(num=1, random_state=0)[0]
#Q = rec(rotation[0].as_quat(), n=(4, 1))
#rot_ang, rot_axis = Q.unit_quaternion_as_axis_and_angle()

# get a random rotation matrix
#rot = Rotation.random(num=1, random_state=1)[0].as_dcm()

# parameters for refinement
#S.set_orientation(np.hstack((a_real, b_real, c_real))*1e-10, length_unit=1)
S.set_orientation(np.hstack((a_real, b_real, c_real)), length_unit=1e-10)
X = S.forward_independent_parameters()
dX = S.forward_gradients()  # gradients in meters

#ar = np.dot(rot, a_real)
#br = np.dot(rot, b_real)
#cr = np.dot(rot, c_real)
#S.set_orientation(np.hstack((ar, br, cr)))
#dX2 = S.forward_gradients()  # same as dX
#
#_ = S.forward_independent_parameters()
#dX3 = S.forward_gradients()  # now its rotated or something...
#embed()
#
#B = S.backward_orientation(independent=X)
#dX3 = S.forward_gradients()
#
#
#B = S.backward_orientation(independent=X)
#embed()

n_ucell_params = len(X)
B = S.backward_orientation(independent=X)
print("Number of parameters = %d" % n_ucell_params)

assert np.allclose(C.get_B(), B.reciprocal_matrix())
#Bstar = B.reciprocal_matrix()

print ("OK!")

#C.rotate_around_origin(rot_axis, rot_ang)

SIM = SimData()
nbcryst = nanoBragg_crystal()
nbcryst.dxtbx_crystal = C
nbcryst.thick_mm = 0.01
nbcryst.Ncells_abc = (20, 20, 20)
SIM.crystal = nbcryst
#SIM.add_water = True
#SIM.water_path_mm = 0.1
#SIM.add_air = False
#SIM.include_noise = True
SIM.instantiate_diffBragg()
SIM.D.add_diffBragg_spots()
img = SIM.D.raw_pixels.as_numpy_array()
SIM.D.free_all()

#ucell2 = (55.2, 66.1, 77.8, 90, 96.1, 90)
#C2 = nanoBragg_crystal.dxtbx_crystal_from_ucell_and_symbol(ucell2, symbol)
#C2.rotate_around_origin(rot_axis, rot_ang)
#a2, b2, c2 = C2.get_real_space_vectors()
#S.set_orientation(a2+b2+c2)
#X2 = S.forward_independent_parameters()
#dX2 = S.forward_gradients()
##n_ucell_params2 = len(X2)

np.random.seed(1)
X_peturb = np.random.normal(X, 1e-10)


# muck around with the X parameters, toggle one at a time and compare finite differences to analytical derivs
for i_param in range(n_ucell_params):

    X2 = list(X)
    X2[i_param] += 1e-8 #X_peturb[i_param]
    B2 = S.backward_orientation(independent=X2).direct_matrix()
    a2_real = B2[0], B2[1], B2[2]
    b2_real = B2[3], B2[4], B2[5]
    c2_real = B2[6], B2[7], B2[8]
    C2 = Crystal(a2_real, b2_real, c2_real, symbol)
    print("Peturbing parameter %d" % i_param)
    print("Ground truth unit cell:")
    print C.get_unit_cell()
    print("Peturbed unit cell:")
    print C2.get_unit_cell()

    # simulate the peturbed crystal image
    nbcryst.dxtbx_crystal = C2
    SIM.crystal = nbcryst
    SIM.instantiate_diffBragg()
    SIM.D.add_diffBragg_spots()
    img2 = SIM.D.raw_pixels.as_numpy_array()
    SIM.D.free_all()
    # done.

    gt_ucell = ",".join(map(lambda x: "%1.3f" % x, C.get_unit_cell().parameters()))
    perturbed_ucell = ",".join(map(lambda x: "%1.3f" % x, C2.get_unit_cell().parameters()))

    if args.plot:
        imgs = cycle((img, img2))
        titles = cycle(("GT: %s" % gt_ucell,"perturbed: %s" % perturbed_ucell))
        icount = 0
        while icount < 6:
            plt.cla()
            plt.imshow(imgs.next(), vmax=100)
            plt.title(titles.next() + "\nContinuing in %d" % (6-icount))
            plt.draw()
            plt.pause(.7)
            icount += 1
            plt.close()

    # Simulate the derivative
    roi = ((310, 320), (155, 161))

    # reset the ground truth crystal to compute the analytical derivative
    nbcryst.dxtbx_crystal = C
    SIM.crystal = nbcryst
    SIM.instantiate_diffBragg()
    D = SIM.D
    D.region_of_interest = roi
    D.printout_pixel_fastslow = (roi[0][0], roi[1][0])
    D.refine(3+i_param)
    D.initialize_managers()

    # TODO: verify if we need the transpose or not
    #z = dX[i_param]
    #z2 = sqr(z).transpose()
    #dX_transpose = flex.double(z2.elems)

    D.set_ucell_derivative_matrix(3+i_param, dX[i_param]/1e10)
    D.add_diffBragg_spots()
    analy_deriv = D.get_derivative_pixels(3+i_param).as_numpy_array()
    SIM.D.free_all()

    delta_param = X2[i_param] - X[i_param]
    (x1, x2), (y1, y2) = roi
    finite_deriv = (img2[y1:y2+1, x1:x2+1] - img[y1:y2+1, x1:x2+1]) / delta_param


    #titles = cycle(("GT: %s" % gt_ucell, "perturbed: %s" % perturbed_ucell))
    embed()
    toggle_image()
    print("OK")
