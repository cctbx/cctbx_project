"""Tests for optimizations using the kramkron API"""

import numpy as np

import example_opt

def test_run_example_opt_params_0():
    """Test optimization with Kramers-Kronig penalty yields decreasing actual loss"""
    _,_,_,_,_,_,_,_,_,_,_,\
    actual_loss_vec  = example_opt.run_example_opt(width=5,
                                                   padn=100,
                                                   trim=30,
                                                   spacing=20,
                                                   noise_loc=[0,0],
                                                   noise_scale=[1e-3,1e-3],
                                                   learning_rate=1e-1,
                                                   num_iter=1000,
                                                   )
    np.testing.assert_array_less(actual_loss_vec[-1].detach().numpy(),actual_loss_vec[0].detach().numpy())
    
def test_run_example_opt_params_1():
    """Test optimization with Kramers-Kronig penalty yields decreasing actual loss"""
    _,_,_,_,_,_,_,_,_,_,_,\
    actual_loss_vec  = example_opt.run_example_opt(width=5,
                                                   padn=100,
                                                   trim=30,
                                                   spacing=20,
                                                   noise_loc=[0,0],
                                                   noise_scale=[1e-2,1e-2],
                                                   learning_rate=1e-1,
                                                   num_iter=1000,
                                                   )
    np.testing.assert_array_less(actual_loss_vec[-1].detach().numpy(),actual_loss_vec[0].detach().numpy())
    
def test_run_example_opt_params_2():
    """Test optimization with Kramers-Kronig penalty yields decreasing actual loss"""
    _,_,_,_,_,_,_,_,_,_,_,\
    actual_loss_vec  = example_opt.run_example_opt(width=5,
                                                   padn=100,
                                                   trim=30,
                                                   spacing=5,
                                                   noise_loc=[0,0],
                                                   noise_scale=[1e-3,1e-3],
                                                   learning_rate=1e-1,
                                                   num_iter=1000,
                                                   )
    np.testing.assert_array_less(actual_loss_vec[-1].detach().numpy(),actual_loss_vec[0].detach().numpy())