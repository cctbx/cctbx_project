"""Tests for kramkron core_functions"""

import pytest
import numpy as np

import core_functions

@pytest.fixture
def Fe3(): # example 0
    path = "sample_data/pf-rd-ox_fftkk.out"
    return(core_functions.parse_data(path))

@pytest.fixture
def Fe2(): #example 1
    path = "sample_data/pf-rd-red_fftkk.out"
    return(core_functions.parse_data(path))

@pytest.fixture
def Fe0():
    path = "sample_data/fe.nff"
    return(core_functions.parse_data(path,remove_first_line=True))


def test_parse_Fe3_beginning(Fe3):
    """Test that input is parsed properly."""
    np.testing.assert_equal(Fe3[:4], np.array([[1006.0,-1.95442043311, 5.92009170594],
                                                    [1007.0, -2.81223119888, 8.52503033764],
                                                    [1008.0, -3.33235759037, 10.1156545543],
                                                    [1009.0, -3.56113395273, 10.8295300422]]))

    
def test_parse_Fe3_end(Fe3):
    """Test that input is parsed properly."""
    np.testing.assert_equal(Fe3[-4:], np.array([[24902.0,0.244888423888, 0.418575531827],
                                                     [24903.0,0.237396106135, 0.405759333696],
                                                     [24904.0,0.220802089543, 0.377388435958],
                                                     [24905.0,0.18454034239, 0.315405661584]]))


def test_parse_Fe0_beginning(Fe0):
    """Test that input is parsed properly."""
    np.testing.assert_equal(Fe0[:4], np.array([[10.0000, -9999.00, 1.37852], 
                                                    [10.1617, -9999.00, 1.42961],
                                                    [10.3261, -9999.00, 1.48259],
                                                    [10.4931, -9999.00, 1.53754]]))

    
def test_parse_Fe0_end(Fe0):
    """Test that input is parsed properly."""
    np.testing.assert_equal(Fe0[-4:], np.array([[28590.2, 26.2151, 0.333497], 
                                                     [29052.6, 26.2100, 0.323310],
                                                     [29522.5, 26.2050, 0.313422],
                                                     [30000.0, 26.2000, 0.303827]]))


def test_get_f_p_Fe3(Fe3):
    """Test that the Hilbert transform and the approach taken by the Sherrell thesis match."""
    padn=5000
    Z=26
    include_Z_term=False
    energy = Fe3[:,0]
    f_dp = Fe3[:,2]
    energy_interp,f_p_pred,_,_,_ = core_functions.get_f_p(energy, f_dp, padn=padn,
                                     Z = Z, # atomic number
                                     include_Z_term=include_Z_term,
                                     hilbert_transform_func=core_functions.hilbert_transform)
    
    energy_interp,f_p_pred_sherrell,_,_,_ = core_functions.get_f_p(energy, f_dp, padn=padn,
                                              Z = Z, # atomic number
                                              include_Z_term=include_Z_term,
                                              hilbert_transform_func=core_functions.hilbert_transform_sherrell)
    
    np.testing.assert_allclose(f_p_pred, f_p_pred_sherrell)


def test_get_f_p_get_f_dp_Fe3(Fe3):
    """Test that finding f_p and then finding f_dp yields the input f_dp."""
    padn=5000
    crop=500
    Z=26
    include_Z_term=False
    energy = Fe3[:,0]
    f_dp = Fe3[:,2]
    
    energy_interp,f_p_pred, energy_interp_pad, f_p_pred_pad, f_dp_pad = core_functions.get_f_p(energy, f_dp, padn=padn,
                                     Z = Z, # atomic number
                                     include_Z_term=include_Z_term,
                                     hilbert_transform_func=core_functions.hilbert_transform)
    
    
    f_dp_interp = core_functions.INTERP_FUNC(energy, f_dp)(energy_interp)
    
    energy_interp_pad,f_dp_pred,_,_ = core_functions.get_f_dp(energy_interp_pad, f_p_pred_pad, padn=0,
                                                      Z = Z, # atomic number
                                                      include_Z_term=include_Z_term,
                                                      hilbert_transform_func=core_functions.hilbert_transform)

    
    # add back DC term
    F_dp_pred = np.fft.fft(f_dp_pred)
    F_dp_pred[0] = np.fft.fft(f_dp_pad)[0]
    f_dp_pred = np.fft.ifft(F_dp_pred).real
    
    if padn != 0:
        f_dp_pred = f_dp_pred[padn:-padn]

    np.testing.assert_allclose(f_dp_interp[crop:-crop], f_dp_pred[crop:-crop], atol=1e-4, rtol=1e-4)


def test_penalty_Fe3(Fe3):
    """Test that finding f_p and then calculating the penalty yields 0 penalty"""
    padn=5000
    Z=26
    include_Z_term=False
    energy = Fe3[:,0]
    f_dp = Fe3[:,2]
    
    energy_interp,f_p_pred,energy_interp_pad,f_p_pred_pad,f_dp_pad = core_functions.get_f_p(energy, f_dp, padn=padn,
                                                                            Z = Z, # atomic number
                                                                            include_Z_term=include_Z_term,
                                                                            hilbert_transform_func=core_functions.hilbert_transform)
    
    mse = core_functions.penalty(energy_interp_pad, f_p_pred_pad, 
                                 f_dp_pad, 
                                 padn=0, trim=0, Z=Z, include_Z_term=include_Z_term,
                                 hilbert_transform_func=core_functions.hilbert_transform)
    np.testing.assert_allclose(0,mse,atol=1e-8, rtol=1e-8)
    
    
def test_get_f_p_get_f_dp_Fe0(Fe0):
    """Test that finding f_p and then finding f_dp yields the input f_dp."""
    padn=5000
    crop=500
    Z=26
    include_Z_term=False
    energy = Fe0[:,0]
    f_dp = Fe0[:,2]
    
    
    energy_interp,f_p_pred, energy_interp_pad, f_p_pred_pad, f_dp_pad = core_functions.get_f_p(energy, f_dp, padn=padn,
                                                                                               Z = Z, # atomic number
                                                                                               include_Z_term=include_Z_term,
                                                                                               hilbert_transform_func=core_functions.hilbert_transform)
    
    
    f_dp_interp = core_functions.INTERP_FUNC(energy, f_dp)(energy_interp)
    

    
    energy_interp_pad,f_dp_pred,_,_ = core_functions.get_f_dp(energy_interp_pad, f_p_pred_pad, padn=0,
                                                              Z = Z, # atomic number
                                                              include_Z_term=include_Z_term,
                                                              hilbert_transform_func=core_functions.hilbert_transform)

    
    # add back DC term
    F_dp_pred = np.fft.fft(f_dp_pred)
    F_dp_pred[0] = np.fft.fft(f_dp_pad)[0]
    f_dp_pred = np.fft.ifft(F_dp_pred).real
    
    if padn != 0:
        f_dp_pred = f_dp_pred[padn:-padn]
    np.testing.assert_allclose(f_dp_interp[crop:-crop], f_dp_pred[crop:-crop], atol=1e-4, rtol=1e-4)


def test_penalty_Fe0(Fe0):
    """Test that finding f_p and then calculating the penalty yields 0 penalty"""

    padn=5000
    Z=26
    include_Z_term=False
    energy = Fe0[:,0]
    f_dp = Fe0[:,2]
    
    energy_interp,f_p_pred,energy_interp_pad,f_p_pred_pad,f_dp_pad = core_functions.get_f_p(energy, f_dp, padn=padn,
                                                                            Z = Z, # atomic number
                                                                            include_Z_term=include_Z_term,
                                                                            hilbert_transform_func=core_functions.hilbert_transform)
    
    
    # f_dp_interp = core_functions.INTERP_FUNC(energy, f_dp)(energy_interp)
    
    mse = core_functions.penalty(energy_interp_pad, f_p_pred_pad, 
                                 f_dp_pad, 
                                 padn=0, trim=0, Z=Z, include_Z_term=include_Z_term,
                                 hilbert_transform_func=core_functions.hilbert_transform)
    np.testing.assert_allclose(0,mse,atol=1e-8, rtol=1e-8)
    
    
def test_cos_wave():
    """Test that finding f_p when f_dp is cos(energy) yields sin(energy)"""
    energy = np.arange(-np.pi,np.pi,.001)
    
    u = np.cos(energy)
    h_u = np.sin(energy)
    
    energy_interp,f_p_pred,_,_,_ = core_functions.get_f_p(energy, u, padn=0,
                                                    trim=0,
                                                    Z = 26, # atomic number
                                                    include_Z_term=False,
                                                    hilbert_transform_func=core_functions.hilbert_transform,
                                                    window_type='cosine',
                                                    )
    np.testing.assert_allclose(h_u, f_p_pred, rtol=1e-3, atol=1e-3)


def test_filter_out_dc():
    """Test filtering out the DC component"""
    energy = np.arange(-np.pi,np.pi,.001)
    
    u = np.cos(energy)
    u_dc = u+1
    u_filtered = core_functions.filter_out_dc(u_dc)
    

    np.testing.assert_allclose(u_filtered, u, rtol=1e-3, atol=1e-3)
    
    
def test_get_f_p_cos_wave(Fe3):
    """Test that finding f_p with an added cos wave is same as subtracting the cos wave, 
    transforming, and adding the known response"""
    
    padn=0
    energy = Fe3[:,0]
    f_dp = Fe3[:,2]
    
    T = energy[-1]-energy[0]
    w = 2*np.pi/T
    u = np.cos(w*energy)
    h_u = np.sin(w*energy)
    
    energy_interp,f_p_pred_subtract, \
    energy_interp_pad, f_p_pred_pad, \
    f_dp_pad = core_functions.get_f_p(energy, f_dp+u, padn=padn,
                                      hilbert_transform_func=core_functions.hilbert_transform,
                                      known_response_energy=energy,
                                      known_response_f_p=h_u,
                                      known_response_f_dp=u,
                                      )
    

    energy_interp,f_p_pred, \
    energy_interp_pad, f_p_pred_pad, \
    f_dp_pad = core_functions.get_f_p(energy, f_dp+u, padn=padn,
                                      hilbert_transform_func=core_functions.hilbert_transform,
                                      known_response_energy=None,
                                      known_response_f_p=None,
                                      known_response_f_dp=None,
                                      )    

    np.testing.assert_allclose(f_p_pred_subtract, f_p_pred, rtol=1e-4, atol=1e-4)


def test_get_f_dp_cos_wave(Fe3):
    """Test that finding f_dp with an added cos wave is same as subtracting the cos wave, 
    transforming, and adding the known response"""
    
    padn=0
    energy = Fe3[:,0]
    f_p = Fe3[:,1]
    
    T = energy[-1]-energy[0]
    w = 2*np.pi/T
    u = np.cos(w*energy)
    h_u = np.sin(w*energy)
    
    energy_interp,f_dp_pred_subtract, \
    energy_interp_pad, f_dp_pred_pad = core_functions.get_f_dp(energy, f_p+h_u, padn=padn,
                                                               hilbert_transform_func=core_functions.hilbert_transform,
                                                               known_response_energy=energy,
                                                               known_response_f_p=h_u,
                                                               known_response_f_dp=u,
                                                               )
    

    energy_interp,f_dp_pred, \
    energy_interp_pad, f_dp_pred_pad = core_functions.get_f_dp(energy, f_p+h_u, padn=padn,
                                                               hilbert_transform_func=core_functions.hilbert_transform,
                                                               known_response_energy=None,
                                                               known_response_f_p=None,
                                                               known_response_f_dp=None,
                                                               )    

    np.testing.assert_allclose(f_dp_pred_subtract, f_dp_pred, rtol=1e-3, atol=1e-3)