from __future__ import division
import numpy as np
from simtbx.nanoBragg.anisotropic_mosaicity import AnisoUmats
from dxtbx.model import Crystal

cryst_dict = dict([('__id__', 'crystal'),
                   ('real_space_a', (-48.93914505851325, -61.4985726090971, 0.23980318971727585)),
                   ('real_space_b', (-27.63556200961052, 72.26768337463876, 13.81410546001183)),
                   ('real_space_c', (-42.92524538136074, 33.14788397044063, -259.2845460893375)),
                   ('space_group_hall_symbol', '-P 6 2')])

cryst = Crystal.from_dict(cryst_dict)

# mosaicities in degrees
eta_a, eta_b, eta_c = 0.02, 0.1, 0.03

# Generate mosaic models with randomized sampling, no derivatives
Nmos = 1000
A = AnisoUmats(num_random_samples=Nmos)

# spherical cap model
U, Uprime, Udblprime = A.generate_Umats(eta=[eta_a,eta_a,eta_a], crystal=cryst,
                                        how=2, verbose=True, compute_derivs=True)
from simtbx.nanoBragg.tst_gaussian_mosaicity import plotter2
U_by2, _,_ =A.generate_Umats(eta_a/2., crystal=cryst,
                             how=2, verbose=True, compute_derivs=False)
plotter2(U, U_by2, False)

# isotropic case should not depend on the crystal model
Unoxtal, _, _ =A.generate_Umats(eta_a, crystal=None, how=2, compute_derivs=False)
for i in range(len(U)):
    assert np.allclose(U[i].elems, Unoxtal[i].elems)

eps = 0.00001  # degrees
Ushift, _,_ = A.generate_Umats(eta_a+eps, how=2)
Ushift_minus, _,_ =A.generate_Umats(eta_a-eps,how=2)
from simtbx.nanoBragg.tst_gaussian_mosaicity2 import check_finite
check_finite(U, Ushift, Uprime, eps)
print("Finite differences check!")

# anisotropic cap model
U,_,_ =A.generate_Umats((eta_a,eta_b,eta_c), cryst,
                 how=1, compute_derivs=False)
Ucryst,_,_ =A.generate_Umats((eta_a,eta_b,eta_c), cryst, transform_eta=True,
                        how=1, compute_derivs=False)
from scitbx.matrix import col
ABC = map(lambda x:col(x).normalize(), cryst.get_real_space_vectors() )
sigs_xyz = []
sigs_abc = []
for v in ABC:
    for ii, Umat_list in enumerate([U,Ucryst]):
        angs = []
        for u in Umat_list:
            vrot = u*v
            angs.append(np.arccos(vrot.dot(v)))
        angs = np.array(angs)
        sig = angs[~np.isnan(angs)].std()*180/np.pi*2
        if ii==0:
            sigs_xyz.append(sig)
        else:
            sigs_abc.append(sig)

# in the anisotropic model, sigmas should inversely correlate with eta_abc
from scipy.stats import pearsonr
eta_abc = eta_a, eta_b, eta_c
cor_abc = pearsonr(sigs_abc, eta_abc)[0]
cor_xyz = pearsonr(sigs_xyz, eta_abc)[0]
print("eta_abc", eta_abc)
print("sigs_abc", sigs_abc, " Eta corr=", cor_abc)
print("sigs_xyz", sigs_xyz, " Eta corr=", cor_xyz)
assert cor_abc < cor_xyz
assert cor_abc < -0.99
print("OK")
