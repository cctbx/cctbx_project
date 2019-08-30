

from tst_diffBragg_rotXYZ import get_diffBragg_instance

D = get_diffBragg_instance()
x1, x2, y1, y2 = 20, 100, 50, 100
#D.region_of_interest = ((x1, x2), (y1, y2))
D.spot_scale = 1e5
D.add_nanoBragg_spots()
nano_pixels = D.raw_pixels.as_numpy_array()

D.raw_pixels *= 0
#D.region_of_interest = ((x1, x2), (y1, y2))
D.vectorize_umats()
D.add_diffBragg_spots()
roi_pixels = D.raw_pixels.as_numpy_array()
#roi_pixels = D.raw_pixels_roi.as_numpy_array()
print nano_pixels.max()
print roi_pixels.max()




