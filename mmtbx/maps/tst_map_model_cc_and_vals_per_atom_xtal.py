from __future__ import division
import iotbx.pdb
import mmtbx.model
import mmtbx.f_model
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from mmtbx.maps import correlation

pdb_str = """
CRYST1  120.670  120.670  136.370  90.00  90.00  90.00 I 4 2 2
ATOM      1  O   HOH A1001       0.385 140.717  50.480  0.50 10.00           O
ATOM      9  O   HOH A1002       6.262 141.644  51.299  0.20 10.00           O
ATOM     24  O   HOH A1004       9.580 133.515  50.443  1.00 10.00           O
TER
ATOM    400  O   HOH A1055      26.182 123.726  20.298  1.00 10.00           O
ATOM    401  O   HOH A1055      25.854 122.538  21.075  1.00 10.00           O
ATOM    402  O   HOH A1055      26.638 121.342  20.536  1.00 10.00           O
ATOM    403  O   HOH A1055      26.375 120.862  19.435  1.00 10.00           O
ATOM    404  O   HOH A1055      24.353 122.353  21.000  1.00 10.00           O
ATOM    405  O   HOH A1055      23.698 121.398  21.929  1.00 10.00           O
ATOM    406  O   HOH A1055      23.055 121.735  23.085  1.00 10.00           O
ATOM    407  O   HOH A1055      23.560 119.980  21.783  1.00 10.00           O
ATOM    408  O   HOH A1055      22.543 120.613  23.678  1.00 10.00           O
ATOM    409  O   HOH A1055      22.835 119.525  22.907  1.00 10.00           O
ATOM    410  O   HOH A1055      23.969 119.044  20.827  1.00 10.00           O
ATOM    411  O   HOH A1055      22.515 118.180  23.107  1.00 10.00           O
ATOM    412  O   HOH A1055      23.658 117.716  21.026  1.00 10.00           O
ATOM    413  O   HOH A1055      22.940 117.293  22.156  1.00 10.00           O
ATOM    414  O   HOH A1056      27.660 120.961  21.291  1.00 10.00           O
ATOM    415  O   HOH A1056      28.617 119.926  20.939  1.00 10.00           O
ATOM    416  O   HOH A1056      28.936 119.084  22.170  1.00 10.00           O
ATOM    417  O   HOH A1056      29.791 119.479  22.957  1.00 10.00           O
ATOM    418  O   HOH A1056      29.969 120.503  20.457  1.00 10.00           O
ATOM    419  O   HOH A1056      29.789 121.456  19.403  1.00 10.00           O
ATOM    420  O   HOH A1056      30.835 119.383  19.884  1.00 10.00           O
ATOM    421  O   HOH A1057      28.254 117.978  22.432  1.00 10.00           O
ATOM    422  O   HOH A1057      28.579 117.228  23.654  1.00 10.00           O
ATOM    423  O   HOH A1057      30.035 116.807  23.734  1.00 10.00           O
ATOM    424  O   HOH A1057      30.644 116.451  22.736  1.00 10.00           O
ATOM    425  O   HOH A1057      27.655 116.014  23.546  1.00 10.00           O
ATOM    426  O   HOH A1057      26.471 116.555  22.790  1.00 10.00           O
ATOM    427  O   HOH A1057      27.090 117.423  21.727  1.00 10.00           O
ATOM    428  O   HOH A1058      30.589 116.855  24.946  1.00 10.00           O
ATOM    429  O   HOH A1058      31.990 116.512  25.116  1.00 10.00           O
ATOM    430  O   HOH A1058      32.223 115.225  25.897  1.00 10.00           O
ATOM    431  O   HOH A1058      33.376 114.801  26.021  1.00 10.00           O
ATOM    432  O   HOH A1058      32.729 117.660  25.819  1.00 10.00           O
ATOM    433  O   HOH A1058      32.742 118.931  24.991  1.00 10.00           O
ATOM    434  O   HOH A1058      32.085 120.062  25.441  1.00 10.00           O
ATOM    435  O   HOH A1058      33.403 118.955  23.779  1.00 10.00           O
ATOM    436  O   HOH A1058      32.118 121.215  24.691  1.00 10.00           O
ATOM    437  O   HOH A1058      33.444 120.109  23.019  1.00 10.00           O
ATOM    438  O   HOH A1058      32.789 121.228  23.486  1.00 10.00           O
ATOM    439  O   HOH A1059      31.173 114.619  26.441  1.00 10.00           O
ATOM    440  O   HOH A1059      31.363 113.477  27.312  1.00 10.00           O
ATOM    441  O   HOH A1059      30.460 112.314  26.941  1.00 10.00           O
ATOM    442  O   HOH A1059      29.379 112.486  26.390  1.00 10.00           O
ATOM    443  O   HOH A1059      31.084 113.855  28.770  1.00 10.00           O
ATOM    444  O   HOH A1059      32.012 114.945  29.291  1.00 10.00           O
ATOM    445  O   HOH A1059      31.573 115.430  30.966  1.00 10.00           O
ATOM    446  O   HOH A1059      32.850 114.619  31.914  1.00 10.00           O
ATOM    447  O   HOH A1060      30.890 111.112  27.331  1.00 10.00           O
ATOM    448  O   HOH A1060      29.940 110.001  27.198  1.00 10.00           O
ATOM    449  O   HOH A1060      28.873 110.111  28.287  1.00 10.00           O
ATOM    450  O   HOH A1060      28.983 110.902  29.213  1.00 10.00           O
ATOM    451  O   HOH A1060      30.638 108.642  27.328  1.00 10.00           O
ATOM    452  O   HOH A1060      31.329 108.669  28.563  1.00 10.00           O
ATOM    453  O   HOH A1060      31.768 108.465  26.324  1.00 10.00           O
END
"""

def run():
  pdb_inp = iotbx.pdb.input(source_info = None, lines = pdb_str)
  model = mmtbx.model.manager(model_input = pdb_inp)
  xrs_good = model.get_xray_structure()
  f_obs = abs(xrs_good.structure_factors(d_min=1.5).f_calc())
  #
  xrs_poor = xrs_good.deep_copy_scatterers()
  sites_cart = xrs_poor.sites_cart()
  sites_cart[2] = (sites_cart[2][0]+5,
                   sites_cart[2][1]+2,
                   sites_cart[2][2]+3)
  xrs_poor = xrs_poor.replace_sites_cart(new_sites = sites_cart)
  #
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    xray_structure = xrs_poor)
  fmodel.update_all_scales()
  print ("r_work=%6.4f r_free=%6.4f"%(fmodel.r_work(), fmodel.r_free()))
  #
  r = correlation.map_model_cc_and_vals_per_atom_xtal(fmodel = fmodel)
  for i, a in enumerate(model.get_hierarchy().atoms()):
    cc = r.ccs[i]
    mv = r.vals[i]
    #print(a.format_atom_record(), cc, mv)
  mean = flex.mean(r.vals[3:])
  a0 = round(mean/r.vals[0],0)
  a1 = round(mean/r.vals[1],0)
  assert approx_equal(a0, 2)
  assert approx_equal(a1, 5)

if(__name__ == "__main__"):
  run()
