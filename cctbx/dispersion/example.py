"""
Example use of kramers_kronig API:
1. Simulate f" based on a very simple model of the K-edge.
2. Use the dispersion relations to calculate f′.
3. Sample both of these curves with Gaussian noise to simulate experimental measurement of the two curves.
4. Use restraint to optimize the parameters. Use automatic differentiation for first-derivatives.
5. Compare the optimized model to the initial ground truth. 
6. Show result in matplotlib.
"""

import kramers_kronig.kramers_kronig_opt as kramers_kronig_opt

energy,\
f_p,\
f_dp,\
energy_ss,\
f_p_noisy_ss,\
f_dp_noisy_ss,\
f_p_pred_0,\
f_dp_pred_0,\
f_p_opt,\
f_dp_opt,\
loss_vec,\
actual_loss_vec  = kramers_kronig_opt.run_example_opt(width=5,
                                                      padn=100,
                                                      trim=30,
                                                      spacing=10,
                                                      noise_loc=[0,0],
                                                      noise_scale=[1e-3,1e-3],
                                                      learning_rate=1e-1,
                                                      num_iter=10000,
                                                      )

kramers_kronig_opt.visualize(energy,
                             f_p,
                             f_dp,
                             energy_ss,
                             f_p_noisy_ss,
                             f_dp_noisy_ss,
                             f_p_pred_0,
                             f_dp_pred_0,
                             f_p_opt,
                             f_dp_opt,
                             loss_vec,
                             actual_loss_vec,
                             )