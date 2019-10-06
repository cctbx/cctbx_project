"""
We need a way to update geometric properties
via dxtbx objects, without re-instantiating the simulator
These unit tests check our functon update_dxtbx_geoms which allows
one to update those models for beam/detector post-instantiation
"""

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
args = parser.parse_args()

from scipy.stats import linregress
import numpy as np
import pylab as plt
from simtbx.diffBragg import sim_data
from dxtbx.model import Panel
from IPython import embed
from scitbx.matrix import sqr, col

# part 1:
# simple k_diffracted derivative
fs = np.array([1, 0, 0])
ss = np.array([0, -1, 0])
k_incident = np.array([1, 0, 0])

# FIXME:
k_incident = np.array([0, 0, -1])

i, j = 10, 10
pixsize_mm = 0.11
detdist_mm = 150
camera_size_mm = 50, 50
origin_z = -detdist_mm
origin_mm = np.array([-camera_size_mm[1]/2., camera_size_mm[0]/2., origin_z])

unit = 1e-3

# diffracted beam
k = origin_mm + fs*i*pixsize_mm + ss*j*pixsize_mm
k *= unit  # convert to meters like the code

# FIXME mod1
k = np.array([.01685, .02055, -.15])


# derivative of k vector w.r.t. detdist
h = .0115 * unit
origin2_mm = np.array([-camera_size_mm[1]/2., camera_size_mm[0]/2., origin_z+h])
k2 = origin2_mm + fs*i*pixsize_mm + ss*j*pixsize_mm
k2 *= unit  # convert to meters

# FIXME mod1
k2 = np.array([.01685, .02055, -.15+h])
fdiff = (k2-k) / h  #/ unit
dk = (0, 0, 1)

print ("\ndiffracted vector:")
print ("finite deriv", fdiff)
print ("deriv", dk)

# take derivative of unit diffracted vector
air_path = np.linalg.norm(k)
unit_k = k / air_path

air_path2 = np.linalg.norm(k2)
unit_k2 = k2 / air_path2

fdiff = (unit_k2 - unit_k) / h #/ unit

deriv = dk/air_path - k/air_path/air_path/air_path*(np.dot(k, dk))
d_unitk = deriv
print ("\nunit diffracted vector:")
print ("finite deriv", fdiff)
print ("deriv", deriv)
print("error = %f" % np.abs(deriv-fdiff).mean())
print("h=%f" % h)


# now take derivative of q_vecor
wavelen = 1.8
q = 1/wavelen * (unit_k - k_incident)
q2 = 1/wavelen * (unit_k2 - k_incident)
fdiff = (q2-q) / h #/ unit

deriv = 1/wavelen * d_unitk
dq = deriv

print ("\nunit momentum transfer vector:")
print ("finite deriv", fdiff)
print ("deriv", deriv)

# take derivative of V-vector
# V = N*(U*B*q- h0)
SIM = sim_data.SimData()
C = SIM.crystal.dxtbx_crystal
B = sqr(C.get_B()).inverse().transpose()
U = sqr(C.get_U())
O = SIM.crystal.Omatrix
UBO = (U*B*O).transpose()

nn = SIM.crystal.Ncells_abc[0]
N = sqr((nn, 0, 0, 0, nn, 0, 0, 0, nn))
h0 = col((4, 3, -1))
V = N*(UBO*col(q) - h0)
V2 = N*(UBO*col(q2)-h0)
fdiff = (V2-V)/h #/unit
deriv = N*(UBO*col(dq))
dV = deriv

print ("\nV-vector")
print("Values", V.elems)
print ("finite deriv", fdiff.elems)
print ("deriv", deriv.elems)

hrad_sqr = np.dot(V.elems, V.elems)
hrad_sqr2 = np.dot(V2.elems, V2.elems)
fd_hrad = (hrad_sqr2 - hrad_sqr) / h

Fl = nn**3*np.exp(-(hrad_sqr / 0.63 * 1))
Fl2 = nn**3*np.exp(-(hrad_sqr2 / 0.63 * 1))
fdiff_Fl = (Fl2-Fl)/h
dH = 2 * np.dot(V.elems, dV.elems)

print ("\nHrad_sqr")
print(" finite diff=%f" % fd_hrad)
print(" derivative=%f" % dH)

deriv_Fl = -1/.63*Fl * dH

print("\ndFlatt")
print ("finite deriv", fdiff_Fl)
print("deriv: %f" % deriv_Fl)

print("\n---------------------------------\n")
# make a simple detector
#det = sim_data.SimData.simple_detector(detector_distance_mm=200, pixelsize_mm=0.025, image_shape=(2400, 2400))
det = sim_data.SimData.simple_detector(detector_distance_mm=150, pixelsize_mm=0.1, image_shape=(512, 512))
orig_idx = 10  # id of origin coordinate in diffBragg

B = SIM.beam.nanoBragg_constructor_beam

# set the detector
SIM.detector = det
SIM.instantiate_diffBragg()
D = SIM.D
D.refine(orig_idx)
D.initialize_managers()
#D.region_of_interest = ((2070, 2130), (70, 150))
D.add_diffBragg_spots()

# get the simulated image and derivative
img = D.raw_pixels_roi.as_numpy_array()
deriv = D.get_derivative_pixels(orig_idx).as_numpy_array()
D.raw_pixels_roi *= 0
D.raw_pixels *= 0

node = det[0]
node_d = node.to_dict()
O = node_d["origin"][0], node_d["origin"][1], node_d["origin"][2]
distance = O[2]

# update the detector distance
percs = [np.power(1e-1*i, 2) for i in range(1, 9)]
all_shifts = []
all_errors = []
for i_shift, p in enumerate(percs):
    # update the detector model
    delta_shift = abs(p*1e-2*distance)
    O2 = O[0], O[1], (distance + delta_shift)
    node_d["origin"] = O2
    det[0] = Panel.from_dict(node_d)
    D.update_dxtbx_geoms(det, B, 0)
    D.add_diffBragg_spots()
    img2 = D.raw_pixels_roi.as_numpy_array()
    D.raw_pixels_roi *= 0
    D.raw_pixels *= 0
    delta_shift_meters = delta_shift*1e-3
    fdiff = (img2-img)/delta_shift_meters
    bragg = img > 1e-2
    error = np.abs(fdiff[bragg] - deriv[bragg]).mean()
    all_errors.append(error)
    all_shifts.append(delta_shift_meters)
    print("Error=%f shift=%f mm, distance=%f mm" % (error, delta_shift, D.distance_mm))
    if args.plot:
        y = slice(40, 65, 1)
        x = slice(415, 437, 1)
        plt.subplot(121)
        plt.imshow(fdiff[y, x])
        plt.title("finite diff")
        plt.subplot(122)
        plt.imshow(deriv[y, x])
        plt.title("analytical")
        plt.draw()
        plt.suptitle("Detdist=%f mm, delta=%f mm, Shift %d / %d"
                     % (distance+delta_shift, delta_shift, i_shift + 1, len(percs)))
        plt.pause(0.8)

if args.plot:
    plt.close()

    plt.plot(all_shifts, all_errors, 'o')
    plt.show()

l = linregress(all_shifts, all_errors)
assert l.rvalue > .99
assert l.slope > 0
assert l.pvalue < 1e-6
print("OK!")
