"""Tests for optimizations using the kramers_kronig_optimize.py"""

from __future__ import division

import numpy as np
import cctbx.dispersion.kramers_kronig.kramers_kronig_optimize as kramers_kronig_optimize

def test_run_example_opt_params_0():
    """Test optimization with Kramers-Kronig penalty yields decreasing actual loss"""
    _,_,_,_,_,_,_,_,_,_,_,\
    actual_loss_vec  = kramers_kronig_optimize.run_example_opt(width=5,
                                                               padn=10,
                                                               trim=30,
                                                               spacing=20,
                                                               noise_loc=[0,0],
                                                               noise_scale=[1e-3,1e-3],
                                                               learning_rate=1e-1,
                                                               num_iter=100,
                                                               )
    np.testing.assert_array_less(actual_loss_vec[-1].detach().numpy(),actual_loss_vec[0].detach().numpy())

def test_run_example_opt_params_1():
    """Test optimization with Kramers-Kronig penalty yields decreasing actual loss"""
    _,_,_,_,_,_,_,_,_,_,_,\
    actual_loss_vec  = kramers_kronig_optimize.run_example_opt(width=5,
                                                               padn=10,
                                                               trim=30,
                                                               spacing=20,
                                                               noise_loc=[0,0],
                                                               noise_scale=[1e-2,1e-2],
                                                               learning_rate=1e-1,
                                                               num_iter=100,
                                                               )
    np.testing.assert_array_less(actual_loss_vec[-1].detach().numpy(),actual_loss_vec[0].detach().numpy())

def test_run_example_opt_params_2():
    """Test optimization with Kramers-Kronig penalty yields decreasing actual loss"""
    _,_,_,_,_,_,_,_,_,_,_,\
    actual_loss_vec  = kramers_kronig_optimize.run_example_opt(width=5,
                                                               padn=10,
                                                               trim=30,
                                                               spacing=5,
                                                               noise_loc=[0,0],
                                                               noise_scale=[1e-3,1e-3],
                                                               learning_rate=1e-1,
                                                               num_iter=100,
                                                               )
    np.testing.assert_array_less(actual_loss_vec[-1].detach().numpy(),actual_loss_vec[0].detach().numpy())

def test_run_example_opt_params_3():
    """Test optimization with Kramers-Kronig penalty yields decreasing actual loss"""
    _,_,_,_,_,_,_,_,_,_,_,\
    actual_loss_vec  = kramers_kronig_optimize.run_example_opt(width=5,
                                                               padn=10,
                                                               trim=30,
                                                               spacing=5,
                                                               noise_loc=[0,0],
                                                               noise_scale=[1e-3,1e-3],
                                                               learning_rate=1e-1,
                                                               num_iter=100,
                                                               uniform_energy=False,
                                                               )
    np.testing.assert_array_less(actual_loss_vec[-1].detach().numpy(),actual_loss_vec[0].detach().numpy())

def run():
    """Run all tests"""
    test_run_example_opt_params_0()
    test_run_example_opt_params_1()
    test_run_example_opt_params_2()
    test_run_example_opt_params_3()
    print("OK")

if __name__ == '__main__':
    run()
