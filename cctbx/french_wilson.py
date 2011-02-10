import sys, math
from libtbx.utils import Sorry
import libtbx.phil

#acentric tables from French-Wilson supplement, 1978
ac_zj =    [ 0.226,0.230,0.235,0.240,0.246,0.251,0.257,0.263,0.270,
             0.276,0.283,0.290,0.298,0.306,0.314,0.323,0.332,0.341,0.351,
             0.362,0.373,0.385,0.397,0.410,0.424,0.439,0.454,0.470,0.487,
             0.505,0.525,0.545,0.567,0.590,0.615,0.641,0.668,0.698,0.729,
             0.762,0.798,0.835,0.875,0.917,0.962,1.009,1.059,1.112,1.167,
             1.226,1.287,1.352,1.419,1.490,1.563,1.639,1.717,1.798,1.882,
             1.967,2.055,2.145,2.236,2.329,2.422,2.518,2.614,2.710,2.808,
             2.906,3.004 ]
ac_zj_sd = [ 0.217,0.221,0.226,0.230,0.235,0.240,0.245,0.250,0.255,
             0.261,0.267,0.273,0.279,0.286,0.292,0.299,0.307,0.314,0.322,
             0.330,0.339,0.348,0.357,0.367,0.377,0.387,0.398,0.409,0.421,
             0.433,0.446,0.459,0.473,0.488,0.503,0.518,0.535,0.551,0.568,
             0.586,0.604,0.622,0.641,0.660,0.679,0.698,0.718,0.737,0.757,
             0.776,0.795,0.813,0.831,0.848,0.865,0.881,0.895,0.909,0.921,
             0.933,0.943,0.953,0.961,0.968,0.974,0.980,0.984,0.988,0.991,
             0.994,0.996 ]
ac_zf =    [ 0.423,0.428,0.432,0.437,0.442,0.447,0.453,0.458,0.464,
             0.469,0.475,0.482,0.488,0.495,0.502,0.509,0.516,0.524,0.532,
             0.540,0.549,0.557,0.567,0.576,0.586,0.597,0.608,0.619,0.631,
             0.643,0.656,0.670,0.684,0.699,0.714,0.730,0.747,0.765,0.783,
             0.802,0.822,0.843,0.865,0.887,0.911,0.935,0.960,0.987,1.014,
             1.042,1.070,1.100,1.130,1.161,1.192,1.224,1.257,1.289,1.322,
             1.355,1.388,1.421,1.454,1.487,1.519,1.551,1.583,1.615,1.646,
             1.676,1.706 ]
ac_zf_sd = [ 0.216,0.218,0.220,0.222,0.224,0.226,0.229,0.231,0.234,
             0.236,0.239,0.241,0.244,0.247,0.250,0.253,0.256,0.259,0.262,
             0.266,0.269,0.272,0.276,0.279,0.283,0.287,0.291,0.295,0.298,
             0.302,0.307,0.311,0.315,0.319,0.324,0.328,0.332,0.337,0.341,
             0.345,0.349,0.353,0.357,0.360,0.364,0.367,0.369,0.372,0.374,
             0.375,0.376,0.377,0.377,0.377,0.376,0.374,0.372,0.369,0.366,
             0.362,0.358,0.353,0.348,0.343,0.338,0.332,0.327,0.321,0.315,
             0.310,0.304 ]

#centric tables from French-Wilson supplement, 1978
c_zj =     [ 0.114,0.116,0.119,0.122,0.124,0.127,0.130,0.134,0.137,
             0.141,0.145,0.148,0.153,0.157,0.162,0.166,0.172,0.177,0.183,
             0.189,0.195,0.202,0.209,0.217,0.225,0.234,0.243,0.253,0.263,
             0.275,0.287,0.300,0.314,0.329,0.345,0.363,0.382,0.402,0.425,
             0.449,0.475,0.503,0.534,0.567,0.603,0.642,0.684,0.730,0.779,
             0.833,0.890,0.952,1.018,1.089,1.164,1.244,1.327,1.416,1.508,
             1.603,1.703,1.805,1.909,2.015,2.123,2.233,2.343,2.453,2.564,
             2.674,2.784,2.894,3.003,3.112,3.220,3.328,3.435,3.541,3.647,
             3.753,3.962 ]
c_zj_sd =  [ 0.158,0.161,0.165,0.168,0.172,0.176,0.179,0.184,0.188,
             0.192,0.197,0.202,0.207,0.212,0.218,0.224,0.230,0.236,0.243,
             0.250,0.257,0.265,0.273,0.282,0.291,0.300,0.310,0.321,0.332,
             0.343,0.355,0.368,0.382,0.397,0.412,0.428,0.445,0.463,0.481,
             0.501,0.521,0.543,0.565,0.589,0.613,0.638,0.664,0.691,0.718,
             0.745,0.773,0.801,0.828,0.855,0.881,0.906,0.929,0.951,0.971,
             0.989,1.004,1.018,1.029,1.038,1.044,1.049,1.052,1.054,1.054,
             1.053,1.051,1.049,1.047,1.044,1.041,1.039,1.036,1.034,1.031,
             1.029,1.028 ]
c_zf =     [ 0.269,0.272,0.276,0.279,0.282,0.286,0.289,0.293,0.297,
             0.301,0.305,0.309,0.314,0.318,0.323,0.328,0.333,0.339,0.344,
             0.350,0.356,0.363,0.370,0.377,0.384,0.392,0.400,0.409,0.418,
             0.427,0.438,0.448,0.460,0.471,0.484,0.498,0.512,0.527,0.543,
             0.560,0.578,0.597,0.618,0.639,0.662,0.687,0.713,0.740,0.769,
             0.800,0.832,0.866,0.901,0.938,0.976,1.016,1.057,1.098,1.140,
             1.183,1.227,1.270,1.313,1.356,1.398,1.439,1.480,1.519,1.558,
             1.595,1.632,1.667,1.701,1.735,1.767,1.799,1.829,1.859,1.889,
             1.917,1.945 ]
c_zf_sd =  [ 0.203,0.205,0.207,0.209,0.211,0.214,0.216,0.219,0.222,
             0.224,0.227,0.230,0.233,0.236,0.239,0.243,0.246,0.250,0.253,
             0.257,0.261,0.265,0.269,0.273,0.278,0.283,0.288,0.293,0.298,
             0.303,0.309,0.314,0.320,0.327,0.333,0.340,0.346,0.353,0.361,
             0.368,0.375,0.383,0.390,0.398,0.405,0.413,0.420,0.427,0.433,
             0.440,0.445,0.450,0.454,0.457,0.459,0.460,0.460,0.458,0.455,
             0.451,0.445,0.438,0.431,0.422,0.412,0.402,0.392,0.381,0.370,
             0.360,0.349,0.339,0.330,0.321,0.312,0.304,0.297,0.290,0.284,
             0.278,0.272 ]

# tables from CTRUNCATE
# Copyright (C) 2006-2008 Norman Stein

# acentric
# look up tables calculated using quadpack
# tables give values from h = -4.0 to h = 3.0 in steps of 0.1
#ac_zj = [
#  0.22561, 0.23037, 0.23531, 0.24046, 0.24581, 0.25139, 0.25720, 0.26327, 0.26959, 0.27620,
#  0.28310, 0.29032, 0.29787, 0.30577, 0.31406, 0.32274, 0.33186, 0.34143, 0.35150, 0.36208,
#  0.37322, 0.38495, 0.39731, 0.41036, 0.42413, 0.43868, 0.45406, 0.47033, 0.48755, 0.50580,
#  0.52514, 0.54564, 0.56740, 0.59050, 0.61503, 0.64108, 0.66876, 0.69817, 0.72942, 0.76262,
#  0.79788, 0.83533, 0.87507, 0.91722, 0.96188, 1.00916, 1.05915, 1.11192, 1.16756, 1.22611,
#  1.28760, 1.35205, 1.41944, 1.48974, 1.56288, 1.63879, 1.71735, 1.79844, 1.88189, 1.96756,
#  2.05525, 2.14478, 2.23597, 2.32863, 2.42258, 2.51764, 2.61365, 2.71046, 2.80794, 2.90596,
#  3.00444 ]
#
#ac_zj_sd = [
#  0.21604, 0.22024, 0.22459, 0.22910, 0.23377, 0.23861, 0.24362, 0.24882, 0.25422, 0.25982,
#  0.26563, 0.27167, 0.27794, 0.28446, 0.29124, 0.29828, 0.30562, 0.31324, 0.32118, 0.32945,
#  0.33805, 0.34701, 0.35634, 0.36606, 0.37618, 0.38671, 0.39768, 0.40910, 0.42099, 0.43335,
#  0.44620, 0.45956, 0.47343, 0.48781, 0.50272, 0.51815, 0.53410, 0.55056, 0.56751, 0.58494,
#  0.60281, 0.62109, 0.63974, 0.65869, 0.67789, 0.69726, 0.71673, 0.73619, 0.75555, 0.77470,
#  0.79353, 0.81192, 0.82977, 0.84696, 0.86339, 0.87895, 0.89357, 0.90718, 0.91972, 0.93117,
#  0.94152, 0.95076, 0.95894, 0.96609, 0.97226, 0.97755, 0.98200, 0.98573, 0.98880, 0.99130,
#  0.99331 ]
#
#ac_zf = [
#  0.42331, 0.42784, 0.43250, 0.43731, 0.44226, 0.44737, 0.45264, 0.45808, 0.46369, 0.46949,
#  0.47548, 0.48168, 0.48809, 0.49473, 0.50160, 0.50873, 0.51611, 0.52378, 0.53173, 0.53999,
#  0.54857, 0.55749, 0.56677, 0.57643, 0.58649, 0.59696, 0.60788, 0.61927, 0.63115, 0.64354,
#  0.65648, 0.66999, 0.68410, 0.69884, 0.71424, 0.73033, 0.74715, 0.76471, 0.78305, 0.80220,
#  0.82218, 0.84301, 0.86472, 0.88732, 0.91082, 0.93522, 0.96052, 0.98672, 1.01379, 1.04171,
#  1.07043, 1.09992, 1.13013, 1.16097, 1.19239, 1.22431, 1.25663, 1.28928, 1.32215, 1.35516,
#  1.38822, 1.42125, 1.45416, 1.48689, 1.51936, 1.55152, 1.58334, 1.61476, 1.64576, 1.67632,
#  1.70644 ]
#
#ac_zf_sd = [
#  0.21545, 0.21753, 0.21966, 0.22185, 0.22409, 0.22639, 0.22874, 0.23116, 0.23363, 0.23618,
#  0.23878, 0.24146, 0.24420, 0.24702, 0.24990, 0.25287, 0.25591, 0.25903, 0.26222, 0.26550,
#  0.26886, 0.27230, 0.27583, 0.27944, 0.28313, 0.28690, 0.29075, 0.29468, 0.29868, 0.30274,
#  0.30687, 0.31106, 0.31529, 0.31956, 0.32386, 0.32816, 0.33246, 0.33673, 0.34095, 0.34510,
#  0.34915, 0.35307, 0.35683, 0.36039, 0.36372, 0.36678, 0.36951, 0.37190, 0.37389, 0.37544,
#  0.37653, 0.37711, 0.37716, 0.37667, 0.37561, 0.37398, 0.37179, 0.36906, 0.36581, 0.36207,
#  0.35789, 0.35332, 0.34841, 0.34323, 0.33783, 0.33228, 0.32663, 0.32095, 0.31529, 0.30968,
#  0.30416 ]
#
# centric
# look up tables calculated using quadpack
# tables give values from h = -4.0 to h = 4.0 in steps of 0.1
#c_zj = [
#  0.11544, 0.11798, 0.12063, 0.12340, 0.12628, 0.12929, 0.13244, 0.13573, 0.13917, 0.14279,
#  0.14657, 0.15055, 0.15473, 0.15912, 0.16374, 0.16862, 0.17376, 0.17918, 0.18492, 0.19099,
#  0.19742, 0.20424, 0.21148, 0.21917, 0.22736, 0.23609, 0.24540, 0.25535, 0.26599, 0.27739,
#  0.28960, 0.30272, 0.31681, 0.33198, 0.34833, 0.36596, 0.38499, 0.40556, 0.42781, 0.45190,
#  0.47799, 0.50627, 0.53692, 0.57016, 0.60619, 0.64523, 0.68752, 0.73327, 0.78271, 0.83605,
#  0.89346, 0.95513, 1.02117, 1.09167, 1.16666, 1.24610, 1.32989, 1.41787, 1.50981, 1.60539,
#  1.70427, 1.80605, 1.91030, 2.01659, 2.12447, 2.23353, 2.34338, 2.45369, 2.56417, 2.67456,
#  2.78470, 2.89443, 3.00367, 3.11236, 3.22047, 3.32800, 3.43496, 3.54139, 3.64732, 3.75278,
#  3.85783, 3.96249, 4.06681, 4.17083, 4.27457, 4.37807, 4.48134, 4.58442, 4.68732, 4.79006,
#  4.89265 ]
#
#c_zj_sd = [
#  0.15779, 0.16105, 0.16444, 0.16796, 0.17162, 0.17543, 0.17939, 0.18351, 0.18781, 0.19229,
#  0.19697, 0.20185, 0.20694, 0.21227, 0.21784, 0.22367, 0.22977, 0.23617, 0.24287, 0.24990,
#  0.25728, 0.26503, 0.27317, 0.28173, 0.29073, 0.30020, 0.31018, 0.32068, 0.33175, 0.34341,
#  0.35571, 0.36867, 0.38233, 0.39673, 0.41191, 0.42790, 0.44473, 0.46244, 0.48106, 0.50060,
#  0.52108, 0.54251, 0.56489, 0.58819, 0.61238, 0.63741, 0.66320, 0.68964, 0.71661, 0.74395,
#  0.77148, 0.79898, 0.82620, 0.85289, 0.87877, 0.90354, 0.92694, 0.94869, 0.96857, 0.98639,
#  1.00200, 1.01532, 1.02636, 1.03515, 1.04181, 1.04651, 1.04946, 1.05089, 1.05106, 1.05021,
#  1.04859, 1.04642, 1.04389, 1.04115, 1.03835, 1.03558, 1.03291, 1.03039, 1.02805, 1.02591,
#  1.02396, 1.02219, 1.02061, 1.01919, 1.01791, 1.01677, 1.01574, 1.01482, 1.01398, 1.01322,
#  1.01253 ]
#
#c_zf = [
#  0.27272, 0.27577, 0.27892, 0.28217, 0.28553, 0.28900, 0.29259, 0.29630, 0.30015, 0.30413,
#  0.30826, 0.31255, 0.31700, 0.32163, 0.32644, 0.33145, 0.33666, 0.34210, 0.34777, 0.35369,
#  0.35988, 0.36635, 0.37312, 0.38022, 0.38767, 0.39549, 0.40370, 0.41234, 0.42144, 0.43103,
#  0.44115, 0.45183, 0.46311, 0.47505, 0.48769, 0.50108, 0.51528, 0.53035, 0.54634, 0.56332,
#  0.58137, 0.60054, 0.62092, 0.64256, 0.66554, 0.68993, 0.71578, 0.74314, 0.77207, 0.80257,
#  0.83468, 0.86836, 0.90359, 0.94031, 0.97842, 1.01780, 1.05830, 1.09974, 1.14193, 1.18465,
#  1.22766, 1.27075, 1.31368, 1.35625, 1.39826, 1.43957, 1.48003, 1.51954, 1.55804, 1.59549,
#  1.63188, 1.66722, 1.70153, 1.73485, 1.76725, 1.79876, 1.82946, 1.85939, 1.88862, 1.91719,
#  1.94516, 1.97257, 1.99946, 2.02587, 2.05182, 2.07736, 2.10249, 2.12726, 2.15168, 2.17576,
#  2.19952 ]
#
#c_zf_sd = [
#  0.20265, 0.20478, 0.20697, 0.20923, 0.21155, 0.21394, 0.21640, 0.21894, 0.22156, 0.22425,
#  0.22704, 0.22991, 0.23288, 0.23595, 0.23912, 0.24240, 0.24579, 0.24930, 0.25293, 0.25669,
#  0.26059, 0.26462, 0.26880, 0.27314, 0.27763, 0.28228, 0.28710, 0.29210, 0.29728, 0.30265,
#  0.30821, 0.31396, 0.31991, 0.32605, 0.33240, 0.33893, 0.34565, 0.35255, 0.35962, 0.36683,
#  0.37417, 0.38159, 0.38908, 0.39658, 0.40403, 0.41138, 0.41855, 0.42546, 0.43200, 0.43809,
#  0.44360, 0.44842, 0.45243, 0.45551, 0.45756, 0.45846, 0.45814, 0.45655, 0.45364, 0.44944,
#  0.44397, 0.43732, 0.42959, 0.42092, 0.41149, 0.40146, 0.39103, 0.38038, 0.36969, 0.35912,
#  0.34880, 0.33886, 0.32936, 0.32037, 0.31193, 0.30405, 0.29671, 0.28990, 0.28360, 0.27776,
#  0.27234, 0.26731, 0.26263, 0.25826, 0.25417, 0.25032, 0.24670, 0.24327, 0.24002, 0.23693,
#  0.23398 ]

master_phil = libtbx.phil.parse("""
  max_bins = 60
    .type = int
    .short_caption = Max. resolution bins
    .help = '''Maximum number of resolution bins'''
  min_bin_size = 40
    .type = int
    .short_caption = Minimum bin size
    .help = '''Minimum number of reflections per bin'''
""")

def fw_acentric(I, sigma_I, mean_intensity) :
  h = (I/sigma_I) - (sigma_I/mean_intensity)
  if (I/sigma_I) < -3.7 or h < -4.0:
    return -1.0, -1.0, -1.0, -1.0
  else:
    if h < 3.0:
      point = 10.0*(h+4.0)
      pt_1 = int(point)
      pt_2 = pt_1 + 1
      delta = point - pt_1
      J = interpolate(pt_1=ac_zj[pt_1],
                      pt_2=ac_zj[pt_2],
                      delta=delta) * sigma_I
      sigma_J = interpolate(pt_1=ac_zj_sd[pt_1],
                            pt_2=ac_zj_sd[pt_2],
                            delta=delta) * sigma_I
      F = interpolate(pt_1=ac_zf[pt_1],
                      pt_2=ac_zf[pt_2],
                      delta=delta) * math.sqrt(sigma_I)
      sigma_F = interpolate(pt_1=ac_zf_sd[pt_1],
                            pt_2=ac_zf_sd[pt_2],
                            delta=delta) * math.sqrt(sigma_I)
    else:
      J = h*sigma_I
      sigma_J = sigma_I
      F = math.sqrt(J)
      sigma_F = 0.5*(sigma_I/F)
    return J, sigma_J, F, sigma_F

def fw_centric(I, sigma_I, mean_intensity) :
  h = (I/sigma_I) - ( sigma_I/(2.0*mean_intensity) )
  if (I/sigma_I) < -3.7 or h < -4.0:
    return -1.0, -1.0, -1.0, -1.0
  else:
    if h < 4.0: #4.0 for truncate tables, 5.0 for ctruncate tables
      point = 10.0*(h+4.0)
      pt_1 = int(point)
      pt_2 = pt_1 + 1
      delta = point - pt_1
      J = interpolate(pt_1=c_zj[pt_1],
                      pt_2=c_zj[pt_2],
                      delta=delta) * sigma_I
      sigma_J = interpolate(pt_1=c_zj_sd[pt_1],
                            pt_2=c_zj_sd[pt_2],
                            delta=delta) * sigma_I
      F = interpolate(pt_1=c_zf[pt_1],
                      pt_2=c_zf[pt_2],
                      delta=delta) * math.sqrt(sigma_I)
      sigma_F = interpolate(pt_1=c_zf_sd[pt_1],
                            pt_2=c_zf_sd[pt_2],
                            delta=delta) * math.sqrt(sigma_I)
    else:
      #adapted from French-Wilson w/ added x^6 term in the expansion
      h_2 = 1.0 / (h*h)
      h_4 = h_2 * h_2
      h_6 = h_2 * h_4
      #posterier of F
      post_F = math.sqrt(h) * (1.0 - (3.0/8.0)*h_2 - (87.0/128.0)*h_4 - (2889.0/1024.0)*h_6)
      #posterier of sigma_F
      post_sig_F = math.sqrt( h * ((1.0/4.0)*h_2 + (15.0/32.0)*h_4 + (273.0/128.0)*h_6) )
      J = h*sigma_I*(1.0 - (1.0/2.0)*h_2 - (3.0/4.0)*h_4 - 3.0*h_6)
      sigma_J = 2.0*sigma_I*post_F*post_sig_F
      F = post_F*math.sqrt(sigma_I)
      sigma_F = post_sig_F*math.sqrt(sigma_I)
  return J, sigma_J, F, sigma_F

def get_mean_intensity(miller_array):
  sum = 0.0
  for d in miller_array.data():
    sum += d
  return (sum / len(miller_array.data()))

# default number of bins is 60, but require that each bin has at least 40 reflections
# if not try again with less bins until condition is satisfied
# function runs recursively
def f_w_binning(miller_array, max_bins=60, min_bin_size=40, log=sys.stderr):
  miller_array.setup_binner(n_bins=max_bins)
  for i_bin in miller_array.binner().range_all():
    sel = miller_array.binner().selection(i_bin)
    bin = miller_array.select(sel)
    if bin.size() > 0:
      if bin.size() < min_bin_size:
        new_max_bins = max_bins - 1
        print >> log, "bin too small, trying %d bins" % new_max_bins
        f_w_binning(miller_array, max_bins=new_max_bins)
        break

def get_bin_centers(miller_array):
  from cctbx.array_family import flex
  centers = flex.double()
  for i_bin in miller_array.binner().range_all():
    sel = miller_array.binner().selection(i_bin)
    bin = miller_array.select(sel)
    bin_center = (bin.d_max_min()[0]+bin.d_max_min()[1])/2
    centers.append(bin_center)
  return centers

def interpolate(pt_1, pt_2, delta):
  return ( ((1.0-delta)*pt_1) + (delta*pt_2) )

def calculate_mean_intensities(miller_array, log=sys.stdout):
  print >> log, "** Calculating bin mean intensity values for each intensity **"
  bin_mean_intensities = miller_array.mean(use_binning=True).data
  bin_centers = get_bin_centers(miller_array=miller_array)
  d_mean_intensities = dict()
  for i_bin in miller_array.binner().range_all():
    sel = miller_array.binner().selection(i_bin)
    bin = miller_array.select(sel)
    if bin.size() > 0:
      bin_center = bin_centers[i_bin]
      for index, d in bin.d_spacings():
        # d is between bin_center[i-1] and bin_center[i]
        if d > bin_center:
          d_1 = bin_centers[i_bin-1]
          d_2 = bin_centers[i_bin]
          m_1 = bin_mean_intensities[i_bin-1]
          m_2 = bin_mean_intensities[i_bin]
          # there is no bin[i-1]
          if m_1 == None:
            mean_i = bin_mean_intensities[i_bin]
            #TO-DO deal with tail
            #d_1 = d_2
            #d_2 = bin_centers[i_bin+1]
            #m_1 = m_2
            #m_2 = bin_mean_intensities[i_bin+1]
            #slope = (m_2-m_1) / (d_2-d_1)
            #width = bin.d_max_min()[0] - bin.d_max_min()[1]
            #d_2 = d_1
            #d_1 = d_1 - width
            #m_2 = m_1
            #m_1 = -1 * (slope*(d_2-d_1)-m_2)
            #delta = d - d_1
            #mean_i = interpolate(pt_1=m_1,
            #                     pt_2=m_2,
            #                     delta=delta)
            #d_mean_intensities[index] = mean_i
          else:
            delta = (d_1 - d) / (d_1 - d_2)
            mean_i = interpolate(pt_1=m_1,
                                 pt_2=m_2,
                                 delta=delta)
            assert (d_1 > d and d > d_2)
            assert ((m_1 > mean_i and mean_i > m_2) or (m_2 > mean_i and mean_i > m_1))
          d_mean_intensities[index] = mean_i
        # d is between bin_center[i] and bin_center[i+1]
        elif d < bin_center:
          d_1 = bin_centers[i_bin]
          d_2 = bin_centers[i_bin+1]
          m_1 = bin_mean_intensities[i_bin]
          m_2 = bin_mean_intensities[i_bin+1]
          # there is no bin[i+1]
          if m_2 == None:
            mean_i = bin_mean_intensities[i_bin]
            #TO-DO deal with tail
            #d_2 = d_1
            #d_1 = bin_centers[i_bin-1]
            #m_2 = m_1
            #m_1 = bin_mean_intensities[i_bin-1]
            #slope = (m_2-m_1) / (d_2-d_1)
            #width = bin.d_max_min()[0] - bin.d_max_min()[1]
            #d_1 = d_2
            #d_2 = d_1 + width
            #m_1 = m_2
            #m_2 = slope*(d_2-d_1)+m_1
            #delta = d - d_1
            #mean_i = interpolate(pt_1=m_1,
            #                     pt_2=m_2,
            #                     delta=delta)
            #d_mean_intensities[index] = mean_i
          else:
            delta = (d_1 - d) / (d_1 - d_2)
            mean_i = interpolate(pt_1=m_1,
                                 pt_2=m_2,
                                 delta=delta)
            assert (d_1 > d and d > d_2)
            assert ((m_1 > mean_i and mean_i > m_2) or (m_2 > mean_i and mean_i > m_1))
            #print d_1, d, d_2, m_1, mean_i, m_2
          d_mean_intensities[index] = mean_i
        # d = the current bin center
        else:
          mean_i = bin_mean_intensities[i_bin]
          d_mean_intensities[index] = mean_i
  return d_mean_intensities

def french_wilson_scale(miller_array, params=None, log=None):
  from cctbx.array_family import flex
  if not miller_array.is_xray_intensity_array():
    raise Sorry("Input array appears to be amplitudes. This method is only appropriate for input intensities.")
  if params == None:
    params = master_phil.extract()
  if log == None:
    log = sys.stdout
  rejected = []
  print >> log, "** Scaling input intensities via French-Wilson Method **"
  print >> log, "Trying %d bins..." % params.max_bins
  f_w_binning(miller_array=miller_array,
              max_bins=params.max_bins,
              min_bin_size=params.min_bin_size,
              log=log)
  print >> log, "Number of bins = %d" % miller_array.binner().n_bins_used()
  new_I = flex.double()
  new_sigma_I = flex.double()
  new_F = flex.double()
  new_sigma_F = flex.double()
  new_indices = flex.miller_index()
  bin_mean_intensities = miller_array.mean(use_binning=True).data
  d_mean_intensities = calculate_mean_intensities(miller_array=miller_array)
  assert len(d_mean_intensities) == miller_array.data().size()
  for i_bin in miller_array.binner().range_all():
    sel = miller_array.binner().selection(i_bin)
    bin = miller_array.select(sel)
    if bin.size() > 0:
      bin_mean_intensity = bin_mean_intensities[i_bin]
      cen = bin.select_centric()
      acen = bin.select_acentric()
      for I, sigma_I, index in zip(cen.data(),
                                   cen.sigmas(),
                                   cen.indices()):
        if (sigma_I <= 0) :
          rejected.append( (index, I, sigma_I, mean_intensity) )
          continue
        mean_intensity = d_mean_intensities[index]
        J, sigma_J, F, sigma_F = fw_centric(
                                   I=I,
                                   sigma_I=sigma_I,
                                   mean_intensity=mean_intensity)
        if J >= 0:
          assert sigma_J >= 0 and F >= 0 and sigma_F >= 0
          new_I.append(J)
          new_indices.append(index)
          new_sigma_I.append(sigma_J)
          new_F.append(F)
          new_sigma_F.append(sigma_F)
        else:
          rejected.append( (index, I, sigma_I, mean_intensity) )
      for I, sigma_I, index in zip(acen.data(),
                                   acen.sigmas(),
                                   acen.indices()):
        if (sigma_I <= 0) :
          rejected.append( (index, I, sigma_I, mean_intensity) )
          continue
        mean_intensity = d_mean_intensities[index]
        J, sigma_J, F, sigma_F = fw_acentric(
                                   I=I,
                                   sigma_I=sigma_I,
                                   mean_intensity=mean_intensity)
        if J >= 0:
          assert sigma_J >= 0 and F >= 0 and sigma_F >= 0
          new_I.append(J)
          new_indices.append(index)
          new_sigma_I.append(sigma_J)
          new_F.append(F)
          new_sigma_F.append(sigma_F)
        else:
          rejected.append( (index, I, sigma_I, mean_intensity) )
  f_obs = miller_array.customized_copy(indices=new_indices,
                                       data=new_F,
                                       sigmas=new_sigma_F)
  f_obs.set_observation_type_xray_amplitude()
  show_rejected_summary(rejected=rejected, log=log)
  return f_obs

def show_rejected_summary(rejected, log=sys.stderr):
  print >> log, "** Summary or rejected intensities **"
  print >> log, "-----------------------------------------------------------------"
  print >> log, "Miller Index  :  Intesity  :  Sigma  :  Bin Mean Intensity"
  for rej in rejected:
    print >> log, "%s    %.3f      %.3f    %.3f" % \
                  (str(rej[0]),rej[1],rej[2],rej[3])
  print >> log, "-----------------------------------------------------------------"
  print >> log, "** Total # rejected intensities: %d **" % len(rejected)
