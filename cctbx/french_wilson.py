import sys, math
from libtbx.str_utils import make_sub_header
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

def fw_acentric(
      I,
      sigma_I,
      mean_intensity,
      sigma_iobs_rejection_criterion) :
  h = (I/sigma_I) - (sigma_I/mean_intensity)
  h_min = sigma_iobs_rejection_criterion
  i_sig_min = h_min+0.3
  if (I/sigma_I) < i_sig_min or h < h_min:
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

def fw_centric(
      I,
      sigma_I,
      mean_intensity,
      sigma_iobs_rejection_criterion) :
  h = (I/sigma_I) - ( sigma_I/(2.0*mean_intensity) )
  h_min = sigma_iobs_rejection_criterion
  i_sig_min = h_min+0.3
  if (I/sigma_I) < i_sig_min or h < h_min:
    return -1.0, -1.0, -1.0, -1.0
  else:
    if h < 4.0:
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
def f_w_binning(miller_array, max_bins=60, min_bin_size=40, log=None):
  if log == None:
    log = sys.stdout
  bin_success = False
  while not bin_success:
    miller_array.setup_binner(n_bins=max_bins)
    bin_success = True
    for i_bin in miller_array.binner().range_all():
      sel = miller_array.binner().selection(i_bin)
      bin = miller_array.select(sel)
      if bin.size() > 0:
        if bin.size() < min_bin_size:
          max_bins = max_bins - 1
          if max_bins == 0:
            print >> log, "not enough reflections for accurate binning\n"+ \
                          "** skipping French-Wilson scaling **"
            return False
          print >> log, "bin too small, trying %d bins" % max_bins
          bin_success = False
          break
          #f_w_binning(miller_array, max_bins=new_max_bins, log=log)
  return True

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

def calculate_mean_intensities(miller_array, log=None):
  if log == None:
    log = sys.stdout
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

def french_wilson_scale(
      miller_array,
      params=None,
      sigma_iobs_rejection_criterion=None,
      log=None):
  from cctbx.array_family import flex
  if not miller_array.is_xray_intensity_array():
    raise Sorry("Input array appears to be amplitudes. This method is only appropriate for input intensities.")
  if miller_array.unit_cell() is None:
    raise Sorry("No unit cell information found. Please supply unit cell data.")
  if miller_array.crystal_symmetry() is None:
    raise Sorry("No crystal symmetry information found. Please supply "+
                "crystal symmetry data.")
  if params == None:
    params = master_phil.extract()
  if (params.max_bins is None) : # XXX reset in case of user error
    params.max_bins = 60
  if log == None:
    log = sys.stdout
  if (sigma_iobs_rejection_criterion is None) :
    sigma_iobs_rejection_criterion = -4.0
  elif ((sigma_iobs_rejection_criterion < -4.0) or
        (sigma_iobs_rejection_criterion > -1.0)) :
    raise Sorry(
      "For French and Wilson scaling, sigma_iobs_rejection_criterion " +
      "must be a value between -4.0 and -1.0, or None.")
  rejected = []
  make_sub_header("Scaling input intensities via French-Wilson Method",
    out=log)
  print >> log, "Trying %d bins..." % params.max_bins
  if not f_w_binning(miller_array=miller_array,
              max_bins=params.max_bins,
              min_bin_size=params.min_bin_size,
              log=log):
    return None
  print >> log, "Number of bins = %d" % miller_array.binner().n_bins_used()
  new_I = flex.double()
  new_sigma_I = flex.double()
  new_F = flex.double()
  new_sigma_F = flex.double()
  new_indices = flex.miller_index()
  bin_mean_intensities = miller_array.mean(use_binning=True).data
  d_mean_intensities = \
    calculate_mean_intensities(miller_array=miller_array, log=log)
  assert len(d_mean_intensities) == miller_array.data().size()
  for i_bin in miller_array.binner().range_all():
    sel = miller_array.binner().selection(i_bin)
    bin = miller_array.select(sel)
    if bin.size() > 0:
      #bin_mean_intensity = bin_mean_intensities[i_bin]
      cen = bin.select_centric()
      acen = bin.select_acentric()
      for I, sigma_I, index in zip(cen.data(),
                                   cen.sigmas(),
                                   cen.indices()):
        mean_intensity = d_mean_intensities[index]
        if (sigma_I <= 0) :
          if I <= 0 or sigma_I < 0 :
            rejected.append( (index, I, sigma_I, mean_intensity) )
            continue
          else:
            J = I
            sigma_J = sigma_I
            F = math.sqrt(I)
            sigma_F = sigma_I
        else :
          J, sigma_J, F, sigma_F = fw_centric(
                                     I=I,
                                     sigma_I=sigma_I,
                                     mean_intensity=mean_intensity,
                                     sigma_iobs_rejection_criterion=\
                                     sigma_iobs_rejection_criterion)
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
          if I <= 0 or sigma_I < 0 :
            rejected.append( (index, I, sigma_I, mean_intensity) )
            continue
          else:
            J = I
            sigma_J = sigma_I
            F = math.sqrt(I)
            sigma_F = sigma_I
        else :
          mean_intensity = d_mean_intensities[index]
          J, sigma_J, F, sigma_F = fw_acentric(
                                     I=I,
                                     sigma_I=sigma_I,
                                     mean_intensity=mean_intensity,
                                     sigma_iobs_rejection_criterion=\
                                     sigma_iobs_rejection_criterion)
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

def show_rejected_summary(rejected, log=None):
  if log == None:
    log = sys.stdout
  print >> log, "** Total # rejected intensities: %d **" % len(rejected)
  if len(rejected) > 0:
    print >> log, "** Summary or rejected intensities **"
    print >> log, "-----------------------------------------------------------------"
    print >> log, "Miller Index  :  Intesity  :  Sigma  :  Bin Mean Intensity"
    for rej in rejected:
      print >> log, "%s    %.3f      %.3f    %.3f" % \
                    (str(rej[0]),rej[1],rej[2],rej[3])
    print >> log, "-----------------------------------------------------------------"
