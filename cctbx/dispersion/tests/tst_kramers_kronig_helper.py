"""Tests for kramers_kronig_helper.py"""

from __future__ import division

import numpy as np
import torch

import cctbx.dispersion.kramers_kronig.kramers_kronig_helper as kramers_kronig_helper

# Get constants
Fe3 = kramers_kronig_helper.parse_data(kramers_kronig_helper.SAMPLE_DATA_PATH + "/pf-rd-ox_fftkk.out")
Fe2 = kramers_kronig_helper.parse_data(kramers_kronig_helper.SAMPLE_DATA_PATH + "/pf-rd-red_fftkk.out")

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

def test_interpolate_scipy(Fe3):
    """Test that interpolate with scipy mode results in a uniform x array."""
    Fe3 = Fe3[np.array([0,100,101,102,103,104]),:]
    x_new, y_new = kramers_kronig_helper.interpolate(Fe3[:,0],Fe3[:,1], mode="scipy")
    dx_new = x_new[1:]-x_new[:-1]
    np.testing.assert_equal(False,np.any(dx_new-dx_new[0]))

def test_interpolate_torch_uniformity(Fe3):
    """Test that interpolate with torch mode results in a uniform x array."""
    Fe3 = Fe3[np.array([0,100,101,102,103,104]),:]
    x_new, y_new = kramers_kronig_helper.interpolate(torch.tensor(Fe3[:,0]),torch.tensor(Fe3[:,1]), mode="torch")
    dx_new = x_new[1:]-x_new[:-1]
    np.testing.assert_equal(False,np.any(np.array(dx_new-dx_new[0])))

def test_interpolate_torch_0(Fe3):
    """Test that interpolate_torch yields correct result on linear example."""
    x_original = torch.tensor([1,5,10])
    y_original = torch.tensor([2,10,20])
    x_new = torch.tensor([4,7])
    y_new = kramers_kronig_helper.interpolate_torch(x_original, y_original, x_new)
    np.testing.assert_equal(np.array(y_new),np.array([8,14]))

def test_interpolate_torch_1(Fe3):
    """Test that interpolate_torch yields correct result on linear example."""
    x_original = torch.tensor([1,2,3,4,5,6,7,8,9,10])
    y_original = 2*torch.tensor([1,2,3,4,5,6,7,8,9,10])
    x_new = torch.tensor([4,7])
    y_new = kramers_kronig_helper.interpolate_torch(x_original, y_original, x_new)
    np.testing.assert_equal(np.array(y_new),np.array([8,14]))

def run():
    """Run all tests"""
    test_parse_Fe3_beginning(Fe3)
    test_parse_Fe3_end(Fe3)
    test_interpolate_scipy(Fe3)
    test_interpolate_torch_uniformity(Fe3)
    test_interpolate_torch_0(Fe3)
    test_interpolate_torch_1(Fe3)
    print("OK")

if __name__ == '__main__':
    run()
