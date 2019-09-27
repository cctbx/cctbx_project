from simtbx.diffBragg import sim_data
from scitbx.matrix import sqr, col
from cctbx import uctbx
from dxtbx.model import Crystal
from simtbx.diffBragg import nanoBragg_crystal
import numpy as np

ucell = (70, 60, 50, 90.0, 110, 90.0)
symbol = "C121"
a_real, b_real, c_real = sqr(uctbx.unit_cell(ucell).orthogonalization_matrix()).transpose().as_list_of_lists()
C = Crystal(a_real, b_real, c_real, symbol)

nbr = nanoBragg_crystal.nanoBragg_crystal()
nbr.dxtbx_crystal = C

S = sim_data.SimData()
S.crystal = nbr
S.instantiate_diffBragg()
S.D.add_diffBragg_spots()
img = S.D.raw_pixels.as_numpy_array()

# simulate the primitive cell directly
to_p1 = C.get_space_group().info().change_of_basis_op_to_primitive_setting()
Cp1 = C.change_basis(to_p1)
nbr2 = nanoBragg_crystal.nanoBragg_crystal()
nbr2.dxtbx_crystal = Cp1

S2 = sim_data.SimData()
S2.crystal = nbr2
S2.instantiate_diffBragg()
S2.D.add_diffBragg_spots()
img2 = S2.D.raw_pixels.as_numpy_array()

# rescale because currently volume is computed incorrectly
img2 = img2 * S.D.spot_scale / S2.D.spot_scale

assert S.D.Omatrix == tuple(to_p1.c_inv().r().transpose().as_double())
assert S2.D.Omatrix == (1, 0, 0, 0, 1, 0, 0, 0, 1)
assert np.allclose(img, img2)

print("OK")



