from __future__ import absolute_import, division, print_function
from cctbx import french_wilson
from cctbx.development import random_structure
from scitbx.array_family import flex
import boost_adaptbx.boost.python as bp
from six.moves import zip
fw_ext = bp.import_ext("cctbx_french_wilson_ext")
from libtbx.utils import null_out, Sorry
from libtbx.test_utils import Exception_expected
import random

def exercise_00():
  x = flex.random_double(1000)
  y = flex.random_double(1000)
  xa = flex.double()
  ya = flex.double()
  ba = flex.bool()
  for x_, y_ in zip(x,y):
    scale1 = random.choice([1.e-6, 1.e-3, 0.1, 1, 1.e+3, 1.e+6])
    scale2 = random.choice([1.e-6, 1.e-3, 0.1, 1, 1.e+3, 1.e+6])
    b = random.choice([True, False])
    x_ = x_*scale1
    y_ = y_*scale2
    v1 = fw_ext.expectEFW(eosq=x_, sigesq=y_, centric=b)
    v2 = fw_ext.expectEsqFW(eosq=x_, sigesq=y_, centric=b)
    assert type(v1) == type(1.)
    assert type(v2) == type(1.)
    xa.append(x_)
    ya.append(y_)
    ba.append(b)
  fw_ext.is_FrenchWilson(F=xa, SIGF=ya, is_centric=ba, eps=0.001)

def exercise_01():
  """
  Sanity check - don't crash when mean intensity for a bin is zero.
  """
  xrs = random_structure.xray_structure(
    unit_cell=(50,50,50,90,90,90),
    space_group_symbol="P1",
    n_scatterers=1200,
    elements="random")
  fc = abs(xrs.structure_factors(d_min=1.5).f_calc())
  fc = fc.set_observation_type_xray_amplitude()
  cs = fc.complete_set(d_min=1.4)
  ls = cs.lone_set(other=fc)
  f_zero = ls.array(data=flex.double(ls.size(), 0))
  f_zero.set_observation_type_xray_amplitude()
  fc = fc.concatenate(other=f_zero)
  sigf = flex.double(fc.size(), 0.1) + (fc.data() * 0.03)
  fc = fc.customized_copy(sigmas=sigf)
  try :
    fc_fc = french_wilson.french_wilson_scale(miller_array=fc, log=null_out())
  except Sorry :
    pass
  else :
    raise Exception_expected
  ic = fc.f_as_f_sq().set_observation_type_xray_intensity()
  fc_fc = french_wilson.french_wilson_scale(miller_array=ic, log=null_out())

# (centric, h = I/sigI, sigI, <I>, F, SIGF): exact French-Wilson posterior
# moments F = <sqrt(J)>, SIGF = sd(sqrt(J)), computed with mpmath (30
# digits, adaptive quadrature over u = J/sigI with weight
# u^p exp(-(u - h)^2/2), p = -1/2 centric, 0 acentric; then scaled by
# sqrt(sigI)). <I> is chosen so that the prior term sigI/<I> is O(1).
french_wilson_reference = [
  (False, 15.4724887407, 2962.11391564, 19077.3821844, 213.970176194, 6.93093844489),
  (False, 16.383699352, 17.9040398534, 138.748221383, 17.1190042807, 0.52354424971),
  (False, 23.1682367932, 9.57944554571, 32.2405431498, 14.8941372725, 0.321772680024),
  (True, 23.7182696341, 601.569005297, 711.27746253, 119.369586483, 2.52231709079),
  (False, 28.8370067771, 457.865166965, 5371.30171853, 114.889021708, 1.99339251616),
  (True, -2.50497568066, 192.359572332, 244.098251375, 4.59344735172, 3.35966686543),
  (True, 4.98411945095, 6.15485423483, 39.3690085536, 5.44748136753, 0.581619222368),
  (False, 24.8000952411, 180.443959637, 2336.91529412, 66.8820289838, 1.34966058514),
  (False, 18.8608345524, 117.74965657, 358.243318068, 47.1093637548, 1.25085482032),
  (False, 29.8578241747, 1658.12265381, 28132.4725791, 222.472554497, 3.72788856501),
  (True, 4.57897475306, 36.8205139565, 48.7618932879, 12.7260389089, 1.50148956433),
  (False, 10.2131935623, 1732.6900042, 8131.36855568, 132.866549674, 6.54051787761),
  (False, 24.961222519, 5.01885692834, 11.6123810774, 11.1904632229, 0.224360057399),
  (False, 12.5095801085, 4365.62902154, 21401.4748304, 233.504697648, 9.36707772586),
  (True, 17.7720121037, 1082.6861479, 3185.30318056, 138.548123081, 3.9143272794),
  (True, 7.97532564029, 3901.20393964, 80918.8884617, 175.318100878, 11.2331806565),
  (True, 5.13080231347, 10.04867803, 12.768927726, 7.06951648873, 0.730101028086),
  (False, 2.86337823055, 238.152655901, 1425.97545652, 25.6974419272, 4.81054177886),
  (True, 21.1525091186, 12.3558109706, 162.225439742, 16.1529161629, 0.382950041795),
  (True, 10.8849353691, 21.7553975263, 64.01015852, 15.3390273724, 0.712669073932),
  (False, 23.5125796017, 40.8700768134, 1407.90579572, 30.9923488595, 0.659732201545),
  (True, 10.0110630234, 1828.52418721, 23827.7725155, 134.781905727, 6.82337952499),
  (True, 29.6469560178, 21.812232089, 61.2878212639, 25.4187619849, 0.429333870214),
  (False, 7.85552904308, 38.7207973754, 51.9338643209, 17.4046481138, 1.11828994),
  (True, 16.2302483395, 26.7922239521, 296.855969116, 20.823030919, 0.64473190882),
  (True, 11.9558674559, 3770.28362445, 26102.9507536, 211.749207835, 8.93902713882),
  (False, 25.5953470375, 17.6788150207, 32.750175045, 21.2678747886, 0.415821648935),
  (False, 23.9874643322, 28.0198433932, 59.8664121403, 25.9197322888, 0.540806862716),
  (False, 28.0333615711, 19.4419306899, 869.54972366, 23.3420028463, 0.416624282876),
  (False, 16.9166290879, 91.9055368348, 139.228991917, 39.4128084823, 1.16722118375),
]

def exercise_invert_french_wilson():
  """ (I, sigI, <I>) -> exact French-Wilson (F, SIGF) -> back. """
  for i, (c, h, s, mi, f_fw, sigf_fw) in enumerate(french_wilson_reference):
    i_true = s * (h + (0.5 if c else 1.0) * s / mi)
    valid, i_obs, sig_i_obs, h_inv, prior_dominated = \
      fw_ext.invert_french_wilson(F=f_fw, SIGF=sigf_fw, mean_intensity=mi,
        centric=c)
    assert valid and not prior_dominated, i
    assert abs(sig_i_obs / s - 1) < 1e-5, i
    assert abs(i_obs - i_true) < 1e-5 * s * max(1, abs(h)), i
  # SIGF/F at or beyond the large-sigma limit: prior-dominated
  for sigf, c in [(0.6, False), (0.8, True)]:
    r = fw_ext.invert_french_wilson(F=1.0, SIGF=sigf, mean_intensity=1.0,
      centric=c)
    assert r[4], (sigf, c)

if (__name__ == "__main__"):
  exercise_00()
  exercise_01()
  exercise_invert_french_wilson()
  print("OK")
