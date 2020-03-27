"""
We need a way to update geometric properties
via dxtbx objects, without re-instantiating the simulator
These unit tests check our functon update_dxtbx_geoms which allows
one to update those models for beam/detector post-instantiation
"""

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plotimages", action='store_true')
parser.add_argument("--plotlines", action='store_true')
parser.add_argument("--curvatures",action="store_true")
parser.add_argument("--nopolar", action="store_true")
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
#h_vals = [1e-5*2**i for i in range(20)]

h_vals = [(2*i*(1e-2))*1e-3 for i in range(1,30,2)]

all_q_err = []
all_q_err2 = []
all_uk_err = []
for h in h_vals:
    # derivative of k vector w.r.t. detdist
    #origin2_mm = np.array([-camera_size_mm[1]/2., camera_size_mm[0]/2., origin_z+h])
    #k2 = origin2_mm + fs*i*pixsize_mm + ss*j*pixsize_mm
    #k2 *= unit  # convert to meters

    # FIXME mod1
    k2 = np.array([.01685, .02055, -.15+h])
    k0 = np.array([.01685, .02055, -.15-h])
    fdiff = (k2-k) / h  #/ unit
    dk =np.array( (0, 0, 1))


    #print ("\ndiffracted vector:")
    #print ("finite deriv", fdiff)
    #print ("deriv", dk)

    # take derivative of unit diffracted vector
    air_path = np.linalg.norm(k)
    unit_k = k / air_path

    air_path2 = np.linalg.norm(k2)
    unit_k2 = k2 / air_path2

    air_path0 = np.linalg.norm(k0)
    unit_k0 = k0 / air_path0

    fdiff = (unit_k2 - unit_k) / h  # / unit

    G = np.dot(k, dk)
    per_k = 1./air_path
    per_k3 = np.power(per_k, 3)
    per_k5 = np.power(per_k, 5)
    deriv =dk*per_k - k*(per_k3*G)
    d_unitk = deriv

    #print ("\nunit diffracted vector:")
    #print ("finite deriv", fdiff)
    #print ("deriv", d_unitk)
    err = np.abs(d_unitk-fdiff).mean()
    #print("error = %2.7g, h=%2.7g" % (err, h))
    all_uk_err.append(err)
    #print("h=%2.7g" % h)
    # now take derivative of q_vecor

    wavelen = 1.8
    q = 1 / wavelen * (unit_k - k_incident)
    q2 = 1 / wavelen * (unit_k2 - k_incident)
    q0 = 1 / wavelen * (unit_k0 - k_incident)
    fdiff = (q2 - q) / h  # / unit

    deriv = 1 / wavelen * d_unitk

    dq = 1/wavelen*d_unitk
    err = np.abs(dq-fdiff).mean()
    #print("error = %2.7g, h=%2.7g" % (err, h))
    all_q_err.append(err)

    # second derivative
    fdiff2 = (q2+q0 - 2.*q) / h/h  # / unit

    # analytical
    d_unitk2 = (3*per_k5*G*G)*k - per_k3*(np.dot(dk, dk))*k - (2*per_k3*G)*dk
    dq2 = 1 / wavelen * d_unitk2

    err2 = np.abs(dq2-fdiff2).mean()
    #print("\terror2 = %2.7g, h=%2.7g" % (err, h))
    all_q_err2.append(err2)

    print ("\nmomentum transfer vector:")
    print ("finite deriv2", fdiff2)
    print ("deriv2", dq2)

luk = linregress(h_vals, all_uk_err)
lq = linregress(h_vals, all_q_err)
lq2 = linregress(np.array(h_vals)**2, all_q_err2)

print (luk.rvalue, luk.slope)
print (lq.rvalue, lq.slope)
print (lq2.rvalue, lq2.slope)
print ("\nmomentum transfer vector:")
print ("finite deriv", fdiff)
#print ("deriv", deriv)
print ("deriv", dq)

print ("\nmomentum transfer vector:")
print ("finite deriv2", fdiff2)
print ("deriv", dq2)

assert luk.rvalue > 0.99
assert lq.rvalue > 0.99
assert lq2.rvalue > 0.99


wavelen = 1.8
h = h_vals[0]
k2 = np.array([.01685, .02055, -.15 + h])
k0 = np.array([.01685, .02055, -.15 - h])
dk = np.array((0, 0, 1))
air_path = np.linalg.norm(k)
unit_k = k / air_path
air_path2 = np.linalg.norm(k2)
unit_k2 = k2 / air_path2
air_path0 = np.linalg.norm(k0)
unit_k0 = k0 / air_path0

q = 1/wavelen * (unit_k - k_incident)
q2 = 1/wavelen * (unit_k2 - k_incident)
q0 = 1/wavelen * (unit_k0 - k_incident)
fdiff = (q2-q) / h

per_k = 1/air_path
per_k3 = np.power(per_k, 3)
per_k5 = np.power(per_k, 5)
G = np.dot(dk, k)
d_unitk = dk * per_k - k * (per_k3 * G)
dq = 1/wavelen * d_unitk

print ("\nunit momentum transfer vector:")
print ("finite deriv", fdiff)
print ("deriv", dq)


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
V0 = N*(UBO*col(q0)-h0)
fdiff = (V2-V)/h #/unit
fdiff2 = (V2+V0-2*V)/h/h #/unit


deriv = N*(UBO*col(dq))
dV = deriv

d_unitk2 = (3 * per_k5 * G * G) * k - per_k3 * (np.dot(dk, dk)) * k - (2 * per_k3 * G) * dk
dq2 = 1 / wavelen * d_unitk2
dV2 = N*(UBO*col(dq2))

print ("\nV-vector")
print("Values", V.elems)
print ("finite deriv", fdiff.elems)
print ("deriv", dV.elems)
print ("finite deriv2", fdiff2.elems)
print ("deriv2", dV2.elems)

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

#deriv_Fl2 = -2/.63*((dV2*V)*Fl + (dV*dV)*Fl + (dV*V)*deriv_Fl)

print("\n---------------------------------\n")

# make a simple detector
#det = sim_data.SimData.simple_detector(detector_distance_mm=200, pixelsize_mm=0.025, image_shape=(2400, 2400))

##############################


#   BEGIN FINITE DIFF TEST


##############################

det = sim_data.SimData.simple_detector(detector_distance_mm=150, pixelsize_mm=0.1, image_shape=(512, 512))
orig_idx = 10  # id of origin coordinate in diffBragg

B = SIM.beam.nanoBragg_constructor_beam

# set the detector
SIM.detector = det
SIM.instantiate_diffBragg()
D = SIM.D

D.oversample_omega = False
D.nopolar = args.nopolar

D.refine(orig_idx)
D.initialize_managers()
#D.region_of_interest = ((2070, 2130), (70, 150))
D.add_diffBragg_spots()

# get the simulated image and derivative
img = D.raw_pixels_roi.as_numpy_array()
deriv = D.get_derivative_pixels(orig_idx).as_numpy_array()
deriv2 = D.get_second_derivative_pixels(orig_idx).as_numpy_array()
D.raw_pixels_roi *= 0
D.raw_pixels *= 0

node = det[0]
node_d = node.to_dict()
O = node_d["origin"][0], node_d["origin"][1], node_d["origin"][2]
distance = O[2]

# update the detector distance
percs = [np.power(1e-3*i, 2) for i in range(1, 30, 2)]
all_shifts = []
all_errors = []
all_errors2 = []
shifts_mm = [2*i*(1e-2) for i in range(1,30,2)]
print (shifts_mm)

for i_shift, p in enumerate(percs):
    # update the detector model
    #delta_shift = abs(p*1e-2*distance)
    delta_shift = shifts_mm[i_shift]  # + distance

    O2 = O[0], O[1], (distance + delta_shift)
    node_d["origin"] = O2

    det[0] = Panel.from_dict(node_d)
    D.update_dxtbx_geoms(det, B, 0)
    D.add_diffBragg_spots()
    img_forward = D.raw_pixels_roi.as_numpy_array()
    D.raw_pixels_roi *= 0
    D.raw_pixels *= 0
    delta_shift_meters = delta_shift*1e-3
    fdiff = (img_forward-img)/delta_shift_meters

    if args.curvatures:
        O2 = O[0], O[1], (distance - delta_shift)
        node_d["origin"] = O2
        det[0] = Panel.from_dict(node_d)
        D.update_dxtbx_geoms(det, B, 0)

        D.add_diffBragg_spots()
        img_backward = D.raw_pixels_roi.as_numpy_array()
        D.raw_pixels_roi *= 0
        D.raw_pixels *= 0

        fdiff2 = (img_forward - 2*img + img_backward) / delta_shift_meters / delta_shift_meters

    bragg = img > 1e-2

    error = (np.abs(fdiff[bragg] - deriv[bragg]) ).mean()
    all_errors.append(error)

    if args.curvatures:
        error2 = (np.abs(fdiff2[bragg] - deriv2[bragg]) / 1).mean()
        all_errors2.append(error2)

    all_shifts.append(delta_shift_meters)

    print("Error=%2.7g shift=%2.7g mm, originZ=%2.7g mm" % (error, delta_shift, D.distance_mm))
    if args.plotimages:
        plt.clf()
        y = slice(40, 65, 1)
        x = slice(415, 437, 1)
        plt.subplot(121)
        plt.imshow(fdiff[y, x])
        plt.title("finite diff")
        plt.subplot(122)
        plt.imshow(deriv[y, x])
        plt.title("analytical")
        plt.draw()
        plt.suptitle("originZ=%f mm, delta=%f mm, Shift %d / %d"
                     % (distance+delta_shift, delta_shift, i_shift + 1, len(percs)))
        plt.pause(0.3)
        if args.curvatures:
            plt.clf()
            plt.subplot(121)
            plt.imshow(fdiff2[y, x])
            plt.title("finite 2nd diff")
            plt.subplot(122)
            plt.imshow(deriv2[y, x])
            plt.title("analytical")
            plt.draw()
            plt.suptitle("originZ=%f mm, delta=%f mm, Shift %d / %d"
                         % (distance + delta_shift, delta_shift, i_shift + 1, len(percs)))
            plt.pause(0.3)

print( "\n\nSecond derivative\nerror          (shift)^2\n-----------------------------")
for err2, shift in zip(all_errors2, all_shifts):
    print ("%10.7g  %10.7g" %(err2, shift**2))

if args.plotlines:
    plt.close()

    plt.plot(all_shifts, all_errors, 'o')
    plt.suptitle("finite 1st difference error", fontsize=16)
    plt.xlabel("h", fontsize=16)
    plt.show()

l = linregress(all_shifts, all_errors)
print("finite diff l.rvalue=%10.7g" % l.rvalue)
assert l.rvalue > .99
assert l.slope > 0
assert l.pvalue < 1e-6

if args.curvatures:
    if args.plotlines:
        plt.close()

        plt.plot(np.array(all_shifts)**2, all_errors2, 'o')
        plt.suptitle("finite 2nd difference error", fontsize=16)
        plt.xlabel("$h^2$", fontsize=16)
        plt.show()

    l = linregress(np.array(all_shifts)**2, all_errors2)
    print("finite 2nd diff l.rvalue=%10.7g" % l.rvalue)
    assert l.rvalue > .99
    assert l.slope > 0
    assert l.pvalue < 1e-6
print("OK!")
