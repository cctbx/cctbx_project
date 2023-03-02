"""Tests for kramers_kroning.py"""

from __future__ import division

import numpy as np
import torch

import cctbx.dispersion.kramers_kronig.kramers_kronig_helper as kramers_kronig_helper
import cctbx.dispersion.kramers_kronig.kramers_kronig as kramers_kronig

# Get constants
Fe3 = kramers_kronig_helper.parse_data(kramers_kronig_helper.SAMPLE_DATA_PATH + "/pf-rd-ox_fftkk.out")[6064:6165,:]
Fe2 = kramers_kronig_helper.parse_data(kramers_kronig_helper.SAMPLE_DATA_PATH + "/pf-rd-red_fftkk.out")[6064:6165,:]

path = kramers_kronig_helper.SAMPLE_DATA_PATH + "/Fe_fake.dat"
sf = kramers_kronig_helper.parse_data(path)
ind_0 = np.argmin(np.abs(sf[:,0]-7070))
ind_1 = np.argmin(np.abs(sf[:,0]-7170))
sf = sf[ind_0:ind_1,:]
Fe0 = sf


def get_f_p_get_f_dp(sf, padn=10):
    """Helper for test that finding f_p and then finding f_dp yields the input f_dp."""

    energy = torch.Tensor(sf[:,0])
    f_dp = torch.Tensor(sf[:,2])

    energy_padded,f_p_pred,f_p_pred_padded,f_dp_padded = kramers_kronig.get_f_p(energy,
                                                                                f_dp,
                                                                                padn=padn,
                                                                                )


    energy_padded,f_dp_pred,f_dp_pred_padded,_ = kramers_kronig.get_f_dp(energy_padded,
                                                                         f_p_pred_padded,
                                                                         padn=0,
                                                                         )


    # add back DC term
    F_dp_pred = np.fft.fft(f_dp_pred)
    F_dp_pred[0] = np.fft.fft(f_dp_padded)[0]
    f_dp_pred = np.fft.ifft(F_dp_pred).real

    if padn != 0:
        f_dp_padded = f_dp_padded[padn:-padn]
        f_dp_pred = f_dp_pred[padn:-padn]

    return(f_dp_padded, f_dp_pred)


def test_get_f_p_get_f_dp_Fe3(Fe3):
    """Test that finding f_p and then finding f_dp yields the input f_dp."""
    f_dp, f_dp_pred = get_f_p_get_f_dp(Fe3)

    np.testing.assert_allclose(f_dp, f_dp_pred, atol=1e-4, rtol=1e-4)


def test_get_f_p_get_f_dp_Fe2(Fe2):
    """Test that finding f_p and then finding f_dp yields the input f_dp."""
    f_dp, f_dp_pred = get_f_p_get_f_dp(Fe2)

    np.testing.assert_allclose(f_dp, f_dp_pred, atol=1e-4, rtol=1e-4)


def get_penalty(sf, padn=10):
    """Helper for test that finding f_p and then calculating the penalty yields 0 penalty"""

    energy = torch.Tensor(sf[:,0])
    f_dp = torch.Tensor(sf[:,2])

    energy_padded,f_p_pred,f_p_pred_padded,f_dp_padded  = kramers_kronig.get_f_p(energy,
                                                                                 f_dp,
                                                                                 padn=padn,
                                                                                 )

    mse = kramers_kronig.get_penalty(energy_padded, f_p_pred_padded,
                                     f_dp_padded,
                                     padn=0, trim=0)
    return(mse)


def test_get_penalty_Fe3(Fe3):
    """Test that finding f_p and then calculating the penalty yields 0 penalty"""
    mse = get_penalty(Fe3)
    np.testing.assert_allclose(0,mse,atol=1e-8, rtol=1e-8)


def test_get_penalty_Fe2(Fe2):
    """Test that finding f_p and then calculating the penalty yields 0 penalty"""
    mse = get_penalty(Fe2)
    np.testing.assert_allclose(0,mse,atol=1e-8, rtol=1e-8)


def test_get_penalty_known_response(Fe2):
    """Test that the penalty of Fe2 is 0 when the known response is the same as the input function"""

    energy = torch.Tensor(Fe2[:,0])
    f_p = torch.Tensor(Fe2[:,1])
    f_dp = torch.Tensor(Fe2[:,2])


    mse = kramers_kronig.get_penalty(energy, f_p,
                                     f_dp,
                                     padn=0, trim=0,
                                     known_response_energy=energy,
                                     known_response_f_p=f_p,
                                     known_response_f_dp=f_dp,
                                     )

    np.testing.assert_allclose(0,mse,atol=1e-8, rtol=1e-8)

def test_cos_wave():
    """Test that finding f_p when f_dp is cos(energy) yields sin(energy)"""
    energy = torch.Tensor(np.arange(-np.pi,np.pi,.05))

    u = torch.Tensor(np.cos(energy))
    h_u = np.sin(energy)

    _,f_p_pred,_,_ = kramers_kronig.get_f_p(energy, u,
                                            padn=0, trim=0,
                                            )

    np.testing.assert_allclose(h_u, f_p_pred, rtol=1e-2, atol=1e-2)


def test_get_f_p_cos_wave_0():
    """Test that finding f_p with an added cos wave is same as subtracting the cos wave,
    transforming, and adding the known response"""

    padn=0
    energy = torch.Tensor(np.arange(7070,7170))

    T = energy[-1]-energy[0]
    w = 2*np.pi/T
    u = torch.Tensor(np.cos(w*energy))

    _,h_u,_,_ = kramers_kronig.get_f_p(energy, u, padn=padn,
                                            known_response_energy=None,
                                            known_response_f_p=None,
                                            known_response_f_dp=None,)

    _,f_p_pred_subtract,_,_ = kramers_kronig.get_f_p(energy, u, padn=padn,
                                                      known_response_energy=energy,
                                                      known_response_f_p=h_u,
                                                      known_response_f_dp=u,
                                                      )

    np.testing.assert_allclose(f_p_pred_subtract, h_u, rtol=1e-4, atol=1e-4)



def test_get_f_p_cos_wave(Fe3):
    """Test that finding f_p with an added cos wave is same as subtracting the cos wave,
    transforming, and adding the known response"""

    padn=0
    energy = torch.Tensor(Fe3[:,0])
    f_dp = torch.Tensor(Fe3[:,2])

    T = energy[-1]-energy[0]
    w = 2*np.pi/T
    u = torch.Tensor(np.cos(w*energy))

    _,h_u,_,_ = kramers_kronig.get_f_p(energy, u, padn=padn,
                                            known_response_energy=None,
                                            known_response_f_p=None,
                                            known_response_f_dp=None,)

    _,f_p_pred_subtract,_,_ = kramers_kronig.get_f_p(energy, f_dp+u, padn=padn,
                                                     known_response_energy=energy,
                                                     known_response_f_p=h_u,
                                                     known_response_f_dp=u,
                                                     )



    _,f_p_pred,_,_ = kramers_kronig.get_f_p(energy, f_dp+u, padn=padn,
                                            known_response_energy=None,
                                            known_response_f_p=None,
                                            known_response_f_dp=None,
                                            )

    _,f_p_pred_subtract,_,_ = kramers_kronig.get_f_p(energy, u, padn=padn,
                                                      known_response_energy=energy,
                                                      known_response_f_p=h_u,
                                                      known_response_f_dp=u,
                                                      )

    _,f_p_pred,_,_ = kramers_kronig.get_f_p(energy, u, padn=padn,
                                            known_response_energy=None,
                                            known_response_f_p=None,
                                            known_response_f_dp=None,
                                            )

    np.testing.assert_allclose(f_p_pred_subtract, f_p_pred, rtol=1e-4, atol=1e-4)

def test_get_f_dp_cos_wave(Fe3):
    """Test that finding f_dp with an added cos wave is same as subtracting the cos wave,
    transforming, and adding the known response"""

    padn=0
    energy = torch.Tensor(Fe3[:,0])
    f_p = torch.Tensor(Fe3[:,1])

    T = energy[-1]-energy[0]
    w = 2*np.pi/T
    u = torch.Tensor(np.cos(w*energy))
    _,h_u,_,_ = kramers_kronig.get_f_p(energy, u, padn=padn)

    _,f_dp_pred_subtract,_,_ = kramers_kronig.get_f_dp(energy, f_p+h_u, padn=padn,
                                                       known_response_energy=energy,
                                                       known_response_f_p=h_u,
                                                       known_response_f_dp=u,
                                                       )

    _,f_dp_pred,_,_ = kramers_kronig.get_f_dp(energy, f_p+h_u, padn=padn,
                                              known_response_energy=None,
                                              known_response_f_p=None,
                                              known_response_f_dp=None,
                                              )

    np.testing.assert_allclose(f_dp_pred_subtract, f_dp_pred, rtol=1e-2, atol=1e-2)

def test_get_f_p_nonuniform(Fe0, padn=100):
    """Test that finding f_p for Fe_nonuniform is the same as finding it for Fe_uniform"""

    # remove some values to make non-uniform
    energy_0 = torch.Tensor(np.concatenate((Fe0[0:2,0],Fe0[2::2,0],Fe0[-1:,0]),axis=0))
    f_dp_0 = torch.Tensor(np.concatenate((Fe0[0:2,2],Fe0[2::2,2], Fe0[-1:,2]),axis=0))

    # original vector
    energy_1 = torch.Tensor(Fe0[:-1,0])
    f_dp_1 = torch.Tensor(Fe0[:-1,2])

    energy_padded_0,f_p_pred_0,_,_ = kramers_kronig.get_f_p(energy_0,
                                              f_dp_0,
                                              padn=padn,
                                              )

    energy_padded_1,f_p_pred_1,_,_ = kramers_kronig.get_f_p(energy_1,
                                              f_dp_1,
                                              padn=padn,
                                              )

    np.testing.assert_allclose(energy_padded_0, energy_padded_1)
    np.testing.assert_allclose(f_p_pred_0, f_p_pred_1, rtol=1e-1, atol=1e-1)

def run():
    """Run all tests"""
    test_get_f_p_get_f_dp_Fe3(Fe3)
    test_get_f_p_get_f_dp_Fe2(Fe2)
    test_get_penalty_Fe3(Fe3)
    test_get_penalty_Fe2(Fe2)
    test_get_penalty_known_response(Fe2)
    test_cos_wave()
    test_get_f_p_cos_wave_0()
    test_get_f_p_cos_wave(Fe3)
    test_get_f_dp_cos_wave(Fe3)
    test_get_f_p_nonuniform(Fe0, padn=100)
    print("OK")

if __name__ == '__main__':
    run()
