"""
This test checks the setter and getter for Ncells parameter
"""

from IPython import embed
from simtbx.diffBragg import sim_data

S = sim_data.SimData()
S.instantiate_diffBragg()
S.D.spot_scale = 1000

Ncells_id = 9

S.D.Ncells_abc = 12
i = S.D.get_value(Ncells_id)
S.D.set_value(Ncells_id, 14)
i2 = S.D.get_value(Ncells_id)
S.D.Ncells_abc = 80
i3 = S.D.get_value(Ncells_id)

assert i == 12
assert i2 == 14
assert i3 == 80

print("OK")

S.D.refine(9)
S.D.initialize_managers()

S.D.region_of_interest = ((0, 511), (0, 511))
S.D.printout_pixel_fastslow = 10, 10
S.D.add_diffBragg_spots()
img = S.D.raw_pixels.as_numpy_array()
deriv = S.D.get_derivative_pixels(9).as_numpy_array()

N = S.D.get_value(9)
delta_N = 0.001
S.D.set_value(9, N+delta_N)

S.D.raw_pixels *= 0
S.D.region_of_interest = ((0, 511), (0, 511))
S.D.add_diffBragg_spots()
img2 = S.D.raw_pixels.as_numpy_array()

fdiff = (img2 - img) / delta_N

import numpy as np
bragg = img > 1
error = np.abs(fdiff[bragg] - deriv[bragg]).mean()
print ("error=%f"%error)
embed()
