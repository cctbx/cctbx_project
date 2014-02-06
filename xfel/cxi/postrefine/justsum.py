"""
Read pickle files output from integration process to determine the polarity and
refine the rotation matrix
Usage:
libtbx.python ~mu238/XfelProject/justsum.py hklin=1jw8-sf-asu.mtz pickle_dir=/net/viper/raid1/mu238/L614/results/myo_stills/015/integration/ frame_start=0 frame_end=757 flag_polar=1

L650-myoglobin 
(wavelength = 1.748, pickles=380) unit_cell=90.309,90.309,45.204,90,90,120 
libtbx.python ~mu238/XfelProject/justsum.py hklin=1jw8-sf-asu.mtz pickle_dir=/net/viper/raid1/mu238/L650/pickles/pickles/ frame_start=0 frame_end=380

L650-hydrogenase
(pickles=178) unit_cell=111.7610,111.7610,103.7700,90,90,90 space_group=P4(2)2(1)2 
libtbx.python ~mu238/XfelProject/justsum.py hklin=3c8y-sf-asu.mtz pickle_dir=/net/viper/raid1/mu238/L650-hydrogenase/results/myo_stills/002/integration frame_start=0 frame_end=178
"""
def get_overall_correlation (data_a, data_b) :
  """
  Correlate the averaged intensities to the intensities from the
  reference data set.
  """

  assert len(data_a) == len(data_b)
  corr = 0
  slope = 0
  try:
    sum_xx = 0
    sum_xy = 0
    sum_yy = 0
    sum_x  = 0
    sum_y  = 0
    N      = 0
    for i in xrange(len(data_a)):
      I_r       = data_a[i]
      I_o       = data_b[i]
      N      += 1
      sum_xx += I_r**2
      sum_yy += I_o**2
      sum_xy += I_r * I_o
      sum_x  += I_r
      sum_y  += I_o
    slope = (N * sum_xy - sum_x * sum_y) / (N * sum_xx - sum_x**2)
    corr  = (N * sum_xy - sum_x * sum_y) / (math.sqrt(N * sum_xx - sum_x**2) *
               math.sqrt(N * sum_yy - sum_y**2))
  except ZeroDivisionError:
    pass
    
  return corr, slope

def get_observations (dir_name,data_subset):
  file_names = []
  for file_name in os.listdir(dir_name):
    if (file_name.endswith("_00000.pickle")):
      if data_subset==0 or \
        (data_subset==1 and (int(os.path.basename(file_name).split("_00000.pickle")[0][-1])%2==1)) or \
        (data_subset==2 and (int(os.path.basename(file_name).split("_00000.pickle")[0][-1])%2==0)):
        file_names.append(os.path.join(dir_name, file_name))
    elif (file_name.endswith(".pickle")):
      if data_subset==0 or \
        (data_subset==1 and (int(os.path.basename(file_name).split(".pickle")[0][-1])%2==1)) or \
        (data_subset==2 and (int(os.path.basename(file_name).split(".pickle")[0][-1])%2==0)):
        file_names.append(os.path.join(dir_name, file_name))
  print "Number of pickle files found:", len(file_names)
  return file_names

def organize_input(observations, anomalous_flag=False):
  """Given the supplied observations and Miller array, return a tuple of
  the original observations, the ASU-mapped observations, and the
  observations reindexed using (k, h, -l).  Merge Bijvoet mates unless
  anomalous_flag is True.
  """
  
  observations = observations.resolution_filter(d_min=d_min, d_max=d_max)
  
  miller_set = symmetry(
      unit_cell=target_unit_cell,
      space_group_symbol=target_space_group
    ).build_miller_set(
      anomalous_flag=anomalous_flag,
      d_min=d_min)
  
  i_I_positive = (observations.data() > 0)
  miller_indices_positive = observations.indices().select(i_I_positive)
  I_positive = observations.data().select(i_I_positive)
  sigI_positive = observations.sigmas().select(i_I_positive)
  
  observations = observations.customized_copy(indices=miller_indices_positive,
      data=I_positive,
      sigmas=sigI_positive,
      anomalous_flag=anomalous_flag,
      crystal_symmetry=miller_set.crystal_symmetry()
      )
  
  I_obs_scaled = (observations.data() - np.mean(observations.data()))/np.std(observations.data())
  i_I_obs_sel = (flex.abs(I_obs_scaled) < sigma_max)
  
  observations = observations.customized_copy(indices=observations.indices().select(i_I_obs_sel), 
      data=observations.data().select(i_I_obs_sel), 
      sigmas=observations.sigmas().select(i_I_obs_sel),
      anomalous_flag=anomalous_flag,
      crystal_symmetry=miller_set.crystal_symmetry()
      )
  
  observations_asu = observations.customized_copy(
      anomalous_flag=anomalous_flag,
      crystal_symmetry=miller_set.crystal_symmetry()
      ).map_to_asu()
  
  observations_original = observations.customized_copy(
      anomalous_flag=anomalous_flag,
      crystal_symmetry=miller_set.crystal_symmetry()
      )
  
  if flag_polar:
    from cctbx import sgtbx
    cb_op = sgtbx.change_of_basis_op('k,h,-l')
    observations_rev = observations_asu.change_basis(cb_op).map_to_asu()
  else:
    observations_rev = observations_asu
  
  return observations_original, observations_asu, observations_rev
  
  
def determine_polar(observations_original, observations_asu, observations_rev, miller_array_main):
  
  """
  determine polarity based on input data
  """
  
  matches_asu = miller.match_multi_indices(
              miller_indices_unique=miller_array_main.indices(),
              miller_indices=observations_asu.indices())
  
  I_ref_asu = flex.double([miller_array_main.data()[pair[0]] for pair in matches_asu.pairs()])
  I_obs_asu = flex.double([observations_asu.data()[pair[1]] for pair in matches_asu.pairs()])
  miller_indices_ori_asu = flex.miller_index((observations_original.indices()[pair[1]] for pair in matches_asu.pairs()))
  
  corr_raw_asu, slope_raw_asu = get_overall_correlation(I_ref_asu, I_obs_asu)
  
  matches_rev = miller.match_multi_indices(
              miller_indices_unique=miller_array_main.indices(),
              miller_indices=observations_rev.indices())
  
  I_ref_rev = flex.double([miller_array_main.data()[pair[0]] for pair in matches_rev.pairs()])
  I_obs_rev = flex.double([observations_rev.data()[pair[1]] for pair in matches_rev.pairs()])
  miller_indices_ori_rev = flex.miller_index((observations_original.indices()[pair[1]] for pair in matches_rev.pairs()))
  
  corr_raw_rev, slope_raw_rev = get_overall_correlation(I_ref_rev, I_obs_rev)
  
  if corr_raw_asu >= corr_raw_rev:
    polar_hkl = 'h,k,l'
    I_ref = I_ref_asu
    I_obs = I_obs_asu
    miller_indices_ori = miller_indices_ori_asu
  else:
    polar_hkl = 'k,h,-l'
    I_ref = I_ref_rev
    I_obs = I_obs_rev
    miller_indices_ori = miller_indices_ori_rev
  
  if len(I_ref) != len(miller_indices_ori):
    print "Error"
    exit()        
    
  return polar_hkl, I_ref, I_obs, miller_indices_ori, corr_raw_asu, corr_raw_rev
  
    
def justsum_mproc(frame_no, miller_array_main):
  """
  Main module to support multi-processor
  """
  
  #1. Organize input
  pickle_filename = frame_files[frame_no]
  trial_results = pickle.load(open(pickle_filename,"rb"))
  observations = trial_results["observations"][0]
  wavelength = trial_results["wavelength"]
        
  observations_original, observations_asu, observations_rev = organize_input(observations, anomalous_flag=target_anomalous_flag)
  
  #2. Determine polarity
  polar_hkl, I_ref, I_obs, miller_indices_ori, corr_raw_asu, corr_raw_rev = determine_polar(observations_original,
            observations_asu, observations_rev, miller_array_main)
  #reset polarity if flag_polar is off
  if flag_polar==False:
    polar_hkl = 'h,k,l'
  
  #3 output mtz file
  if polar_hkl == 'h,k,l':
    observations_sel = observations_asu
  elif polar_hkl == 'k,h,-l':
    observations_sel = observations_rev
  
  return frame_no, (corr_raw_asu, corr_raw_rev), observations_sel.indices(), observations_sel.data(), observations_sel.sigmas(), observations_sel, wavelength
  
  
if __name__=="__main__":

  from iotbx import reflection_file_reader
  from cctbx.array_family import flex
  from cctbx import miller
  from cctbx import crystal
  from cctbx.crystal import symmetry
  from scitbx.matrix import sqr, col
  
  from libtbx.easy_mp import pool_map, get_processes
  from cctbx.crystal_orientation import crystal_orientation
  from scitbx.lstbx import normal_eqns
  
  import matplotlib
  import matplotlib.cm as cm
  import matplotlib.mlab as mlab
  import matplotlib.pyplot as plt
  
  import sys
  import numpy as np
  import os,cPickle as pickle,math
  from datetime import date, datetime, time, timedelta
  
  from mod_polar import polar_manager
  from mod_partiality import partiality_handler
  from mod_util import intensities_scaler
  
  
  #capture starting time
  time_global_start=datetime.now()
  
  #Extracting input
  frame_start = 0
  frame_end = 0
  flag_plot = False
  file_name_in_mtz=''
  pickle_dir = '' 
  d_min = 1.5
  d_max = 45.0
  
  #target_unit_cell = (90.309,90.309,45.204,90,90,120) #L614-myoglobin
  #target_space_group = 'P6'
  
  target_unit_cell = (111.7610,111.7610,103.7700,90,90,90) #hydrogenase
  target_space_group = 'P 42 21 2'
  
  target_anomalous_flag = False
  
  flag_polar = False
  
  n_bins = 25
  sigma_max = 15
  
  for i in range(len(sys.argv)):
    pair=sys.argv[i].split('=')
    if pair[0]=='frame_start':
      frame_start=int(pair[1])
    elif pair[0]=='frame_end':
      frame_end=int(pair[1])
    elif pair[0]=='flag_plot':
      if int(pair[1])==1:
        flag_plot=True
    elif pair[0]=='flag_polar':
      if int(pair[1])==1:
        flag_polar=True
    elif pair[0]=='hklin':
      file_name_in_mtz=pair[1]
    elif pair[0]=='pickle_dir':
      pickle_dir=pair[1]
  
  
  if file_name_in_mtz == '':
    print "Please input mtz file with reference intensity"
    exit() 
  
  reflection_file = reflection_file_reader.any_reflection_file(file_name_in_mtz)
  miller_arrays=reflection_file.as_miller_arrays()
    
  miller_array_main = miller_arrays[0].as_intensity_array()
  
  data_subset = 0
  frame_files = get_observations(pickle_dir,data_subset)
  
  #start comparing the main set with frame data  
  frames_rand_600 = [601, 206, 110, 126, 639, 127, 173, 157, 85, 755, 454, 666, 309, 8, 79, 470, 596, 349, 328, 398, 117, 21, 735, 275, 274, 269, 331, 31, 475, 706, 211, 96, 611, 27, 131, 401, 139, 386, 529, 503, 189, 739, 148, 197, 367, 115, 548, 44, 378, 68, 301, 746, 594, 425, 748, 421, 74, 150, 230, 365, 58, 343, 749, 551, 315, 527, 403, 714, 67, 685, 582, 324, 518, 737, 215, 94, 678, 554, 125, 307, 743, 341, 210, 575, 465, 207, 663, 311, 730, 478, 636, 604, 580, 588, 124, 118, 648, 542, 734, 251, 75, 335, 72, 722, 708, 412, 121, 752, 513, 631, 101, 628, 40, 208, 95, 238, 457, 497, 395, 340, 719, 672, 618, 664, 281, 234, 107, 221, 741, 371, 106, 178, 214, 419, 84, 721, 279, 359, 193, 477, 709, 338, 615, 337, 396, 257, 506, 376, 621, 140, 637, 135, 538, 18, 612, 388, 501, 696, 120, 267, 154, 510, 668, 130, 199, 617, 574, 498, 176, 32, 317, 543, 599, 183, 268, 517, 509, 41, 583, 261, 545, 533, 35, 536, 24, 52, 92, 288, 56, 291, 416, 345, 305, 681, 306, 531, 348, 587, 521, 23, 30, 640, 249, 559, 97, 105, 190, 170, 687, 213, 679, 468, 123, 229, 143, 689, 225, 564, 622, 186, 12, 488, 6, 645, 377, 99, 423, 557, 270, 66, 81, 720, 222, 586, 287, 731, 476, 276, 308, 39, 330, 153, 727, 353, 534, 188, 233, 415, 619, 297, 304, 657, 654, 491, 525, 607, 203, 93, 363, 453, 576, 350, 420, 428, 332, 103, 644, 443, 321, 405, 361, 409, 298, 446, 10, 625, 201, 90, 602, 695, 490, 397, 255, 283, 296, 598, 732, 155, 556, 129, 7, 541, 733, 187, 241, 313, 633, 247, 567, 224, 616, 750, 146, 593, 603, 171, 228, 200, 333, 62, 511, 526, 684, 217, 407, 63, 391, 370, 310, 385, 192, 196, 469, 191, 411, 499, 392, 145, 523, 53, 480, 578, 152, 87, 711, 561, 537, 319, 245, 373, 610, 356, 572, 656, 751, 212, 694, 439, 65, 109, 48, 232, 220, 756, 627, 717, 165, 673, 216, 51, 198, 515, 180, 738, 591, 325, 500, 486, 33, 390, 461, 723, 318, 484, 209, 690, 263, 753, 384, 162, 539, 445, 659, 289, 248, 375, 69, 573, 563, 86, 204, 718, 579, 16, 464, 364, 623, 458, 316, 111, 462, 605, 505, 379, 449, 312, 514, 344, 662, 597, 83, 34, 435, 467, 650, 512, 447, 431, 177, 437, 293, 665, 707, 436, 400, 256, 430, 638, 374, 406, 37, 705, 9, 286, 380, 629, 358, 485, 132, 677, 277, 273, 102, 11, 546, 459, 726, 688, 566, 669, 745, 555, 272, 682, 674, 167, 181, 660, 71, 243, 496, 667, 360, 119, 236, 166, 581, 492, 701, 329, 3, 142, 704, 164, 540, 362, 661, 352, 322, 553, 265, 133, 620, 404, 495, 36, 552, 136, 502, 100, 487, 339, 520, 595, 426, 549, 179, 565, 2, 530, 262, 455, 231, 692, 354, 740, 300, 614, 336, 342, 634, 28, 299, 266, 494, 25, 466, 471, 144, 560, 302, 147, 651, 239, 149, 393, 161, 252, 609, 463, 0, 122, 218, 278, 429, 544, 314, 691, 571, 728, 590, 703, 175, 481, 700, 646, 402, 13, 452, 710, 70, 647, 451, 295, 713, 182, 351, 235, 712, 432, 116, 73, 55, 168, 729, 444, 473, 725, 205, 577, 26, 742, 346, 250, 50, 643, 716, 528, 670, 414, 76, 697, 474, 163, 641, 202, 17, 693, 642, 570, 440, 38, 355, 88, 456, 483, 1, 137, 158, 357]
  frames_rand_500 = [197, 171, 386, 722, 262, 605, 232, 168, 218, 411, 465, 115, 189, 199, 149, 23, 74, 336, 227, 549, 80, 554, 616, 120, 69, 166, 34, 648, 148, 304, 19, 573, 331, 277, 476, 93, 555, 206, 349, 582, 238, 423, 195, 569, 361, 285, 44, 729, 250, 180, 702, 45, 228, 172, 556, 50, 529, 389, 397, 577, 742, 612, 695, 301, 693, 313, 439, 16, 9, 42, 587, 723, 418, 715, 281, 273, 181, 750, 491, 665, 407, 570, 71, 39, 66, 483, 378, 161, 2, 367, 257, 632, 362, 255, 205, 221, 52, 210, 581, 128, 328, 13, 91, 173, 366, 7, 650, 388, 442, 260, 134, 615, 299, 167, 745, 437, 504, 398, 17, 358, 421, 698, 237, 175, 329, 716, 263, 89, 138, 55, 744, 652, 626, 725, 403, 429, 219, 466, 563, 705, 505, 586, 385, 163, 0, 113, 347, 506, 27, 333, 508, 749, 267, 539, 369, 372, 436, 392, 62, 733, 422, 408, 119, 86, 717, 321, 131, 396, 588, 54, 747, 662, 169, 535, 395, 21, 654, 584, 160, 190, 226, 122, 48, 461, 462, 420, 215, 492, 542, 496, 157, 415, 709, 534, 346, 72, 673, 309, 236, 356, 81, 441, 520, 217, 244, 158, 435, 320, 664, 110, 318, 562, 286, 552, 337, 499, 521, 595, 679, 593, 129, 324, 532, 453, 613, 70, 714, 8, 51, 12, 630, 414, 572, 116, 88, 302, 526, 193, 323, 295, 332, 264, 391, 604, 656, 342, 541, 335, 220, 92, 671, 127, 699, 609, 401, 374, 85, 59, 114, 710, 658, 150, 620, 482, 256, 338, 444, 143, 211, 732, 390, 603, 410, 474, 509, 674, 224, 457, 512, 458, 696, 371, 359, 564, 511, 99, 100, 433, 176, 104, 225, 708, 274, 399, 678, 419, 641, 351, 617, 405, 712, 463, 480, 249, 24, 667, 344, 472, 688, 459, 740, 208, 565, 545, 18, 248, 375, 501, 121, 478, 684, 354, 156, 475, 132, 265, 746, 523, 721, 469, 82, 680, 515, 293, 184, 682, 498, 687, 597, 619, 676, 298, 558, 43, 179, 287, 638, 578, 571, 254, 306, 752, 334, 443, 209, 636, 41, 625, 677, 607, 522, 488, 345, 290, 592, 83, 460, 440, 311, 452, 4, 449, 6, 213, 479, 487, 272, 657, 754, 112, 627, 424, 240, 279, 446, 694, 634, 292, 103, 538, 317, 513, 229, 484, 343, 697, 748, 235, 152, 659, 700, 655, 640, 326, 327, 622, 270, 65, 222, 95, 713, 413, 29, 187, 271, 25, 427, 106, 383, 690, 186, 536, 477, 296, 245, 58, 322, 691, 448, 735, 628, 38, 192, 464, 266, 580, 11, 633, 574, 288, 579, 1, 417, 547, 96, 212, 683, 97, 319, 243, 233, 404, 32, 364, 49, 406, 559, 494, 387, 182, 412, 382, 330, 473, 315, 447, 644, 204, 637, 610, 60, 365, 284, 471, 631, 78, 743, 608, 507, 170, 734, 566, 207, 239, 276, 719, 454, 258, 701, 376, 553, 75, 497, 548, 305, 133, 350, 668, 111, 142]
  frames_rand_400 = [470, 358, 667, 525, 3, 135, 648, 487, 291, 163, 294, 456, 308, 547, 725, 744, 373, 504, 72, 202, 146, 395, 727, 478, 465, 149, 552, 642, 534, 302, 15, 79, 123, 722, 29, 120, 177, 659, 2, 117, 688, 588, 28, 239, 541, 74, 211, 706, 371, 480, 571, 214, 274, 334, 133, 527, 208, 570, 340, 285, 339, 696, 157, 76, 186, 535, 500, 36, 513, 195, 320, 522, 108, 277, 741, 502, 0, 400, 425, 347, 17, 687, 631, 596, 112, 325, 4, 417, 636, 514, 249, 393, 284, 394, 621, 309, 196, 440, 179, 336, 140, 412, 305, 210, 429, 83, 349, 628, 656, 702, 675, 299, 322, 293, 226, 714, 353, 380, 273, 606, 341, 609, 677, 127, 413, 310, 415, 620, 689, 198, 574, 13, 497, 281, 217, 9, 136, 116, 430, 578, 558, 93, 383, 212, 60, 296, 185, 591, 418, 594, 323, 665, 171, 311, 61, 420, 20, 654, 498, 19, 361, 102, 147, 443, 427, 370, 242, 372, 450, 331, 730, 472, 272, 406, 247, 192, 449, 388, 652, 582, 483, 95, 703, 271, 754, 585, 584, 348, 197, 477, 318, 614, 22, 344, 118, 635, 188, 597, 555, 141, 356, 719, 668, 62, 8, 564, 99, 230, 661, 419, 170, 156, 14, 486, 382, 319, 550, 106, 520, 343, 11, 84, 701, 669, 252, 342, 73, 231, 266, 645, 262, 748, 617, 401, 664, 581, 364, 391, 121, 579, 595, 92, 52, 484, 379, 625, 144, 662, 228, 189, 115, 260, 224, 329, 332, 265, 390, 75, 304, 583, 16, 267, 35, 86, 543, 65, 182, 561, 407, 473, 641, 283, 559, 542, 694, 126, 404, 81, 474, 312, 301, 21, 464, 59, 363, 138, 180, 238, 203, 37, 492, 653, 148, 316, 655, 168, 330, 630, 257, 376, 261, 598, 264, 152, 32, 438, 629, 612, 399, 516, 119, 109, 649, 300, 10, 351, 604, 158, 155, 751, 124, 241, 684, 481, 282, 670, 85, 632, 357, 276, 503, 55, 433, 471, 286, 392, 577, 398, 306, 377, 566, 243, 485, 166, 77, 436, 47, 488, 269, 40, 409, 250, 567, 683, 508, 125, 698, 650, 460, 154, 445, 172, 536, 717, 206, 66, 592, 178, 619, 63, 110, 651, 104, 573, 544, 657, 569, 338, 187, 368, 697, 169, 350, 280, 452, 236, 410, 494, 275, 512, 718, 215, 221, 663, 422, 461, 24, 164, 451, 386]
  frames_rand_300 = [333, 657, 301, 569, 615, 382, 431, 441, 681, 372, 585, 493, 26, 707, 34, 217, 665, 503, 71, 499, 641, 716, 7, 35, 572, 284, 389, 222, 69, 515, 255, 498, 617, 464, 309, 586, 95, 473, 256, 298, 704, 497, 531, 332, 663, 310, 258, 410, 176, 330, 598, 165, 172, 728, 737, 366, 738, 579, 627, 521, 496, 510, 191, 257, 697, 352, 89, 594, 710, 400, 299, 388, 166, 216, 146, 705, 524, 727, 731, 596, 302, 136, 505, 546, 474, 109, 49, 187, 613, 346, 167, 642, 329, 490, 163, 262, 659, 580, 612, 267, 654, 159, 353, 750, 104, 199, 147, 206, 513, 359, 534, 135, 294, 412, 650, 518, 435, 574, 155, 670, 549, 351, 368, 362, 687, 203, 322, 633, 434, 288, 637, 387, 0, 508, 219, 406, 143, 511, 658, 386, 403, 655, 646, 476, 630, 190, 620, 233, 489, 138, 25, 601, 588, 224, 673, 100, 470, 421, 200, 597, 592, 94, 364, 375, 277, 55, 720, 462, 268, 300, 67, 632, 577, 317, 265, 339, 547, 113, 291, 148, 456, 381, 440, 622, 125, 444, 162, 28, 392, 246, 269, 15, 161, 422, 43, 286, 207, 528, 266, 551, 618, 573, 92, 272, 279, 137, 719, 128, 215, 370, 599, 563, 754, 231, 66, 605, 756, 278, 221, 409, 70, 564, 251, 407, 686, 374, 214, 643, 211, 469, 150, 81, 331, 280, 261, 52, 31, 726, 484, 97, 316, 608, 755, 417, 12, 164, 512, 338, 10, 180, 327, 205, 693, 576, 624, 408, 86, 411, 90, 629, 47, 53, 21, 416, 631, 393, 648, 369, 354, 319, 717, 725, 106, 471, 249, 14, 23, 606, 404, 313, 223, 702, 154, 123, 142, 36, 457, 56, 243, 275, 562, 706, 649, 58, 355, 396, 436, 538, 625, 380]
  frames_rand_200 = [347, 3, 526, 726, 172, 692, 436, 670, 701, 426, 279, 663, 137, 225, 390, 363, 110, 266, 150, 485, 508, 582, 600, 129, 407, 349, 56, 194, 441, 4, 562, 11, 496, 629, 381, 323, 410, 377, 333, 133, 34, 156, 59, 525, 199, 254, 527, 164, 309, 189, 606, 31, 229, 166, 592, 68, 265, 191, 248, 669, 675, 651, 424, 307, 78, 661, 736, 596, 601, 152, 487, 41, 70, 158, 745, 107, 329, 169, 126, 117, 210, 171, 585, 432, 359, 188, 438, 306, 462, 178, 258, 239, 51, 687, 623, 691, 27, 589, 184, 212, 529, 547, 520, 690, 135, 443, 626, 397, 396, 356, 47, 299, 593, 574, 476, 345, 69, 519, 219, 157, 186, 657, 18, 567, 537, 557, 565, 654, 76, 221, 83, 8, 754, 170, 46, 319, 656, 492, 388, 200, 636, 447, 569, 278, 5, 676, 723, 82, 461, 403, 486, 247, 314, 311, 634, 357, 44, 344, 452, 737, 700, 584, 560, 287, 732, 594, 375, 48, 506, 729, 650, 206, 142, 187, 748, 120, 538, 330, 342, 62, 331, 246, 662, 42, 412, 680, 98, 89, 23, 553, 118, 431, 433, 613, 36, 6, 75, 459, 481, 37]
  frames_rand_100 = [212, 380, 396, 445, 71, 455, 182, 553, 492, 732, 612, 218, 421, 260, 65, 77, 412, 626, 700, 590, 89, 441, 57, 171, 125, 654, 354, 448, 487, 64, 216, 68, 72, 728, 243, 78, 4, 280, 321, 390, 564, 210, 646, 620, 417, 1, 431, 495, 381, 499, 283, 551, 289, 184, 488, 504, 82, 9, 144, 652, 531, 284, 336, 708, 484, 247, 369, 598, 98, 534, 599, 497, 438, 587, 521, 165, 116, 115, 674, 465, 563, 738, 242, 450, 391, 663, 337, 503, 415, 404, 45, 332, 343, 684, 588, 399, 146, 265, 407, 481]
  frames_rand_50 = [363, 215, 253, 172, 198, 455, 705, 179, 463, 184, 755, 181, 218, 99, 251, 175, 653, 207, 398, 242, 306, 83, 356, 104, 252, 257, 689, 560, 646, 634, 304, 688, 290, 607, 494, 671, 16, 206, 6, 650, 403, 722, 447, 580, 695, 381, 338, 350, 738, 698]
  frames_rand_25 = [127, 258, 386, 138, 26, 131, 320, 349, 376, 728, 750, 465, 556, 662, 741, 133, 296, 648, 336, 398, 18, 10, 56, 743, 254]
  frames_rand_10 = [615, 162, 275, 316, 467, 623, 293, 255, 357, 482]
  frames = frames_rand_100
  #frames = range(frame_start, frame_end)
    
  def justsum_mproc_wrapper(arg):
    return justsum_mproc(arg, miller_array_main)
  
  justsum_result = pool_map(
        args=frames,
        func=justsum_mproc_wrapper,
        processes=None)
  
  
  cc_correct_polar = flex.double()
  cc_wrong_polar = flex.double()
  miller_indices_all = flex.miller_index()
  I_refined_all = flex.double()
  sigI_refined_all = flex.double()
  cn_accept_frame = 0
  observations_all = []
  wavelength_all = flex.double()
  for frame_no, cc_polar, miller_indices_refined, I_obs_refined, sig_I_obs_refined, observations_sel, wavelength  in justsum_result:
    print frame_no, cc_polar[0], cc_polar[1]
    cn_accept_frame += 1
    observations_all.append(observations_sel)
    wavelength_all.append(wavelength)
    if cc_polar[0] > cc_polar[1]:
      cc_correct_polar.append(cc_polar[0])
      cc_wrong_polar.append(cc_polar[1])
    else:
      cc_correct_polar.append(cc_polar[1])
      cc_wrong_polar.append(cc_polar[0])
    
    for index_new, I_new, sigI_new in zip(miller_indices_refined, I_obs_refined, sig_I_obs_refined):
      miller_indices_all.append(index_new)
      I_refined_all.append(I_new)
      sigI_refined_all.append(sigI_new)
        
  
  miller_array_all = miller_array_main.customized_copy(indices=miller_indices_all,
        data=I_refined_all,
        sigmas=sigI_refined_all)
  
  
  
  #scale and output to mtz files
  inten_scaler = intensities_scaler()
  
  """
  inten_scaler.output_mtz_file_mosflm(target_unit_cell, 
        target_space_group, 
        target_anomalous_flag,
        observations_all,
        wavelength_all,
        miller_array_main,
        'justsum_757_res1p5')
  """
  
  miller_arrays_merge, cc_merge_sets, slope_merge_sets, multiplicities = inten_scaler.output_mtz_files(target_unit_cell, 
        target_space_group, 
        target_anomalous_flag,
        miller_indices_all, 
        I_refined_all, 
        sigI_refined_all,
        flex.double([1]*len(I_refined_all)),
        miller_array_main,
        'justsum_757')
  
   
  miller_array_merge = miller_arrays_merge[0]
  binner_merge = miller_array_merge.setup_binner(n_bins=n_bins)
  binner_merge_indices = binner_merge.bin_indices()
  completeness_merge = miller_array_merge.completeness(use_binning=True)
  
  print 'Completeness (%) Multiplicities (mean, median, std., max, min) and correlation (mean, median, max)'
  bin_multiplicities = []
  for i in range(1,n_bins+1):
    i_binner = (binner_merge_indices == i)
    multiplicities_sel = multiplicities.select(i_binner)
    bin_multiplicities.append((np.mean(multiplicities_sel), np.median(multiplicities_sel), np.std(multiplicities_sel), flex.max(multiplicities_sel), flex.min(multiplicities_sel)))
    
    #calculate C.C.iso per bin
    matches_bin = miller.match_multi_indices(
              miller_indices_unique=miller_array_main.indices(),
              miller_indices=miller_arrays_merge[0].indices().select(i_binner))
  
    I_ref_bin = flex.double([miller_array_main.data()[pair[0]] for pair in matches_bin.pairs()])
    I_obs_bin_mean = flex.double([miller_arrays_merge[0].data().select(i_binner)[pair[1]] for pair in matches_bin.pairs()])
    I_obs_bin_median = flex.double([miller_arrays_merge[1].data().select(i_binner)[pair[1]] for pair in matches_bin.pairs()])
    I_obs_bin_max = flex.double([miller_arrays_merge[2].data().select(i_binner)[pair[1]] for pair in matches_bin.pairs()])
    
    corr_bin_mean, slope_bin_mean = get_overall_correlation(I_obs_bin_mean, I_ref_bin)
    corr_bin_median, slope_bin_median = get_overall_correlation(I_obs_bin_median, I_ref_bin)
    corr_bin_max, slope_bin_max = get_overall_correlation(I_obs_bin_max, I_ref_bin)
    
    completeness = completeness_merge.data[i]
          
    print 'bin ', i, ': %2.4f - %2.4f' %binner_merge.bin_d_range(i), \
      '%3.2f'%(completeness*100), '%4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f'%(np.mean(multiplicities_sel), np.median(multiplicities_sel), np.std(multiplicities_sel), flex.max(multiplicities_sel), flex.min(multiplicities_sel), corr_bin_mean, corr_bin_median, corr_bin_max)
    
  print 'Summary'  
  print 'C.C. correct polar mean=%6.5f median=%6.5f sigma=%6.5f' %(np.mean(cc_correct_polar), np.median(cc_correct_polar), np.std(cc_correct_polar))
  print 'C.C. wrong polar mean=%6.5f median=%6.5f sigma=%6.5f' %(np.mean(cc_wrong_polar), np.median(cc_wrong_polar), np.std(cc_wrong_polar))
  print 'No. of accepted frames:', cn_accept_frame
  print 'No. of all reflections', len(I_refined_all)
  print 'No. of reflections after merge', len(miller_arrays_merge[0].indices())
  print 'Final-merge mean (iso) C.C=%6.5f Slope=%6.5f' %(cc_merge_sets[0], slope_merge_sets[0])
  print 'Final-merge median (iso) C.C=%6.5f Slope=%6.5f' %(cc_merge_sets[1], slope_merge_sets[1])
  print 'Final-merge max (iso) C.C=%6.5f Slope=%6.5f' %(cc_merge_sets[2], slope_merge_sets[2])
      
  """
  plt.scatter(I_ref_all,I_obs_all,s=[10]*len(I_ref_all), marker='o', c='b')
  plt.title('CC Final='+str(cc_all))
  plt.xlabel('Reference intensity')
  plt.ylabel('Observed intensity')
  plt.show()
  """
